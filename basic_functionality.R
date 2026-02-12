library(data.table)
library(tidyverse)
library(glue)

source('~/src/coverage/coordinate_conversion.r')

# Convert genomic (hg38) to cDNA for PAH intronic variants
# PAH is on chr12 MINUS strand

construct_label <- function(location, pmid, total, label_type) {
  if (label_type == 'country (size, study)') {
    label = glue("{location} (N={total}, {pmid})")
  } else if (label_type == 'country, (size)') {
    label = glue("{location} (N={total})")
  } else if (label_type == 'country, (study)') {
    label = glue("{location} ({pmid})")
  }
  label
}

if (file.exists('~/pah/genotype_tables/')) {
  gt_file = '~/pah/genotype_tables/combined.csv'
  bb_file = '~/disease_genetics/versions/v1.0/pku.csv'
  cpku_file = '~/disease_genetics/versions/v1.0/cpku.csv'
  annots_file = '~/disease_genetics/v1.0/annots.csv'
} else {
  gt_file = '~/data/combined.csv'
  bb_file = '~/data/v1.0/pku.csv'
  cpku_file = '~/data/v1.0/cpku.csv'
  annots_file = '~/data/v1.0/annots.csv'
}
# run once when updating carrier variant data
if (F) {
  source('~/src/coverage/protein_annotations.R')
  dfx = df_pku2 %>% select(chrm, pos, ref, alt, gene) %>% distinct()
  annots = get_all_variant_annotations(dfx, skip_api=F) %>% mutate
  annots %>% fwrite(annots_file)
}



df_pku2 = fread(bb_file) %>% mutate(study_type='het carriers') %>% filter(gene == 'PAH')
annots = fread(annots_file)
nc = annots %>% filter(grepl('chr', label))
cod = annots %>% filter(!grepl('chr', label))
annots = annots %>% separate(label, into=c(''))



df_pku = fread(gt_file) %>% mutate(study_type='patients') %>%
  mutate(variant=orig_commonize_variant(aa_codes, variant)) %>%
  select(variant, num, pmid, location, study_type)
#df_pku_ivs = df_pku %>% filter(grepl('^IVS', variant)) %>%
  #mutate(ivs=sub('ntt', 'nt', standardized_variant)) %>%
  #mutate(ivs=sub('[ACGTacg].*', '', ivs)) %>%
  #mutate(ivs=sub('nt\\+', '+', ivs)) %>%
  #mutate(ivs=sub('nt-', '-', ivs)) %>%
  #rowwise() %>%  mutate(cdna=ivs_to_cdna(ivs)) %>% ungroup() %>%
  #mutate(standardized_variant=cdna) %>% select(-ivs, -cdna)
#df_pku = df_pku %>% filter(!grepl('^IVS', standardized_variant)) %>%
  #bind_rows(df_pku_ivs)
df = bind_rows(df_pku,
               NULL)

studies =  df %>%
  filter(variant != 'unk', variant != 'Unknown') %>%
  group_by(location, pmid, study_type) %>% summarize(total=sum(num)) %>% ungroup()

variant_lists = list(small=c('R408W', 'c.1066-11G>A', 'IVS10-11G>A', 'P281L'),
                     larger=c('R408W', 'c.1066-11G>A', 'IVS10-11G>A', 'P281L', 'R243X', 'R261Q')
)
label_type = 'country (size, study)'
chosen_studies = studies %>% select(location, pmid) # will be user-defined

df = df %>%
    filter(variant != 'unk', variant != 'Unknown') %>%
    group_by(location, pmid) %>% mutate(total=sum(num)) %>% ungroup() %>%
    inner_join(chosen_studies, by=c('location'='location', 'pmid'='pmid'))

coverages = NULL
for (lname in names(variant_lists)) {
  vlist = variant_lists[[lname]]
  covs =  df %>%
    mutate(frac=num / total, chosen=variant %in% vlist) %>%
    group_by(location, pmid, study_type) %>% summarize(frac_covered=sum(as.numeric(chosen) * frac),
                                                       total=max(total)) %>% ungroup() %>%
    mutate(pop_covered=frac_covered**2 + 2 * frac_covered * (1 - frac_covered)) %>%
    rowwise() %>% mutate(label=construct_label(location, pmid, total, label_type)) %>% ungroup() %>%
    mutate(variant_set=lname, set_size=length(vlist))
  coverages = bind_rows(coverages, covs)
}

levs = coverages %>% group_by(variant_set) %>% summarize(v=mean(pop_covered)) %>% arrange(-v) %>% .$variant_set
coverages$variant_set = factor(coverages$variant_set, levels=levs)

coverages %>%
  arrange(as.numeric(variant_set)) %>%
  ggplot(aes(x=label, y=pop_covered, fill=variant_set)) +
  geom_bar(stat='identity', alpha=.5, position='identity') +
  theme_minimal() +
  coord_flip() +
  labs(x=label_type, y='estimated fraction of patients treated by selected variants')

