
aa_codes <- tribble(
  ~three, ~one, ~name,
  "Ala", "A", "Alanine",
  "Arg", "R", "Arginine",
  "Asn", "N", "Asparagine",
  "Asp", "D", "Aspartic acid",
  "Cys", "C", "Cysteine",
  "Gln", "Q", "Glutamine",
  "Glu", "E", "Glutamic acid",
  "Gly", "G", "Glycine",
  "His", "H", "Histidine",
  "Ile", "I", "Isoleucine",
  "Leu", "L", "Leucine",
  "Lys", "K", "Lysine",
  "Met", "M", "Methionine",
  "Phe", "F", "Phenylalanine",
  "Pro", "P", "Proline",
  "Ser", "S", "Serine",
  "Thr", "T", "Threonine",
  "Trp", "W", "Tryptophan",
  "Tyr", "Y", "Tyrosine",
  "Val", "V", "Valine",
  "Sec", "U", "Selenocysteine",
  "Pyl", "O", "Pyrrolysine",
  "Asx", "B", "Asp or Asn",
  "Glx", "Z", "Glu or Gln",
  "Ter", "X", "Stop"
)

orig_commonize_variant <- function(aa_codes, variant, unknown='unk') {
  cvar = sub('^p.', '', variant)
  cvar = gsub('\\[', '', cvar)
  cvar = gsub('\\]', '', cvar)
  cvar = gsub('\\(', '', cvar)
  cvar = gsub('\\)', '', cvar)
  cvar = gsub(' ', '', cvar)
  cvar = sub('EX', 'Ex', cvar)
  cvar = sub('→', '>', cvar)
  cvar = sub('=', 'X', cvar)
  cvar = sub('\\*', 'X', cvar)
  cvar = sub('Ivs', 'IVS', cvar)
  for (i in seq(1, nrow(aa_codes))) {
    three = aa_codes$three[i]
    one = aa_codes$one[i]
    cvar = gsub(three, one, cvar)
  }
  cvar = case_when(cvar == '?' ~ unknown,
                   cvar == '-' ~ unknown,
                   cvar == '' ~ unknown,
                   cvar == ' ' ~ unknown,
                   cvar == 'Notdetected' ~ unknown,
                   cvar == 'unk' ~ unknown,
                   cvar == 'n.i.' ~ unknown,
                   cvar == 'c.?' ~ unknown,
                   cvar == 'others' ~ unknown,
                   is.na(cvar) ~ unknown,
                   cvar == 'nomutationdetected' ~ unknown,
                   T ~ cvar
                   )
  cvar = case_when(grepl('^E[0-9]+del', cvar) ~ sub('^E', 'Ex', cvar), # annoying overlap between exon and Glu
                   T ~ cvar)
  cvar = ifelse(grepl('^[0-9]', cvar),  paste0('c.', cvar), cvar)
  cvar = ifelse(grepl('Y204C', cvar),  'Ex6-96A>G', cvar)

  # screw these old naming systems
  cvar = case_when(grepl('IVS4nt5', cvar) ~ sub('IVS4nt5', 'IVS4+5', cvar),
                   grepl('IVS2nt5', cvar) ~ sub('IVS2nt5', 'IVS2+5', cvar),
                   grepl('IVS1nt5', cvar) ~ sub('IVS1nt5', 'IVS1+5', cvar),
                   grepl('IVS8nt1', cvar) ~ sub('IVS8nt1', 'IVS8+1', cvar),
                   grepl('IVS8nt\\+1', cvar) ~ sub('IVS8nt\\+1', 'IVS8+1', cvar),
                   grepl('IVS8nt-7', cvar) ~ sub('IVS8nt-7', 'IVS8-7', cvar),
                   grepl('IVS8nt1', cvar) ~ sub('IVS8nt1', 'IVS8+1', cvar),
                   grepl('IVS8nt7', cvar) ~ sub('IVS8nt7', 'IVS8-7', cvar),
                   grepl('IVS7nt1', cvar) ~ sub('IVS7nt1', 'IVS7+1', cvar),
                   grepl('Ivs4nt', cvar) ~ sub('Ivs4nt', 'IVS4+', cvar),
                   T ~ cvar
                         )
  cvar = case_when(!grepl('^IVS', cvar) ~ cvar,
                       grepl('^IVS2\\+', cvar) ~ sub('IVS2', 'c.168', cvar),
                       grepl('^IVS2-', cvar) ~ sub('IVS2', 'c.169', cvar),
                       grepl('^IVS3\\+', cvar) ~ sub('IVS3', 'c.352', cvar),
                       grepl('^IVS3-', cvar) ~ sub('IVS3', 'c.353', cvar),
                       grepl('^IVS4\\+', cvar) ~ sub('IVS4', 'c.441', cvar),
                       grepl('^IVS4-', cvar) ~ sub('IVS4', 'c.442', cvar),
                       grepl('^IVS5\\+', cvar) ~ sub('IVS5', 'c.509', cvar),
                       grepl('^IVS5-', cvar) ~ sub('IVS5', 'c.510', cvar),
                       grepl('^IVS6\\+', cvar) ~ sub('IVS6', 'c.706', cvar),
                       grepl('^IVS6-', cvar) ~ sub('IVS6', 'c.707', cvar),
                       grepl('^IVS7\\+', cvar) ~ sub('IVS7', 'c.842', cvar),
                       grepl('^IVS7-', cvar) ~ sub('IVS7', 'c.843', cvar),
                       grepl('^IVS8\\+', cvar) ~ sub('IVS8', 'c.912', cvar),
                       grepl('^IVS8-', cvar) ~ sub('IVS8', 'c.913', cvar),
                       grepl('^IVS9\\+', cvar) ~ sub('IVS9', 'c.969', cvar),
                       grepl('^IVS9-', cvar) ~ sub('IVS9', 'c.970', cvar),
                       grepl('^IVS10\\+', cvar) ~ sub('IVS10', 'c.1065', cvar),
                       grepl('^IVS10-', cvar) ~ sub('IVS10', 'c.1066', cvar),
                       grepl('^IVS11\\+', cvar) ~ sub('IVS11', 'c.1199', cvar),
                       grepl('^IVS11-', cvar) ~ sub('IVS11', 'c.1200', cvar),
                       grepl('^IVS12\\+', cvar) ~ sub('IVS12', 'c.1315', cvar),
                       grepl('^IVS12-', cvar) ~ sub('IVS12', 'c.1316', cvar),
                       grepl('^IVS1\\+', cvar) ~ sub('IVS1', 'c.60', cvar),
                       grepl('^IVS1-', cvar) ~ sub('IVS1', 'c.61', cvar),
                       T ~ cvar )
  cvar = case_when(!grepl('^IVS', cvar) ~ cvar,
                   cvar == 'IVS10' ~ 'c.1066-11G>A',
                   cvar == 'IVS12' ~ 'c.1315+1G>A',
                   cvar == 'IVS12nt1' ~ 'c.1315+1G>A',
                   T ~ cvar)
  cvar
}

comminze_ivs <- function(ivar) {
  if (!grepl('^IVS', ivar)) {
    return(ivar)
  }
  ivar = case_when(!grepl('^IVS', ivar) ~ ivar,
                       grepl('^IVS2\\+', ivar) ~ sub('IVS2', 'c.168', ivar),
                       grepl('^IVS2-', ivar) ~ sub('IVS2', 'c.169', ivar),
                       grepl('^IVS3\\+', ivar) ~ sub('IVS3', 'c.352', ivar),
                       grepl('^IVS3-', ivar) ~ sub('IVS3', 'c.353', ivar),
                       grepl('^IVS4\\+', ivar) ~ sub('IVS4', 'c.441', ivar),
                       grepl('^IVS4-', ivar) ~ sub('IVS4', 'c.442', ivar),
                       grepl('^IVS5\\+', ivar) ~ sub('IVS5', 'c.509', ivar),
                       grepl('^IVS5-', ivar) ~ sub('IVS5', 'c.510', ivar),
                       grepl('^IVS6\\+', ivar) ~ sub('IVS6', 'c.706', ivar),
                       grepl('^IVS6-', ivar) ~ sub('IVS6', 'c.707', ivar),
                       grepl('^IVS7\\+', ivar) ~ sub('IVS7', 'c.842', ivar),
                       grepl('^IVS7-', ivar) ~ sub('IVS7', 'c.843', ivar),
                       grepl('^IVS8\\+', ivar) ~ sub('IVS8', 'c.912', ivar),
                       grepl('^IVS8-', ivar) ~ sub('IVS8', 'c.913', ivar),
                       grepl('^IVS9\\+', ivar) ~ sub('IVS9', 'c.969', ivar),
                       grepl('^IVS9-', ivar) ~ sub('IVS9', 'c.970', ivar),
                       grepl('^IVS10\\+', ivar) ~ sub('IVS10', 'c.1065', ivar),
                       grepl('^IVS10-', ivar) ~ sub('IVS10', 'c.1066', ivar),
                       grepl('^IVS11\\+', ivar) ~ sub('IVS11', 'c.1199', ivar),
                       grepl('^IVS11-', ivar) ~ sub('IVS11', 'c.1200', ivar),
                       grepl('^IVS12\\+', ivar) ~ sub('IVS12', 'c.1315', ivar),
                       grepl('^IVS12-', ivar) ~ sub('IVS12', 'c.1316', ivar),
                       grepl('^IVS1\\+', ivar) ~ sub('IVS1', 'c.60', ivar),
                       grepl('^IVS1-', ivar) ~ sub('IVS1', 'c.61', ivar),
                       T ~ ivar )
  ivar = case_when(!grepl('^IVS', ivar) ~ ivar,
                   ivar == 'IVS10' ~ 'c.1066-11G>A',
                   ivar == 'IVS12' ~ 'c.1315+1G>A',
                   ivar == 'IVS12nt1' ~ 'c.1315+1G>A',
                   T ~ ivar)
  ivar
}

commonize_variant_annots <- function(amino_acids, protein_position) {
  if (is.na(amino_acids)) {
    return(NULL)
  } else {
    amino_acids = sub('\\*', 'X', amino_acids)
    if ((protein_position == 204) & (amino_acids == 'Y/C')) {
      cvar = 'Ex6-96A>G'
    } else {
      amino_acids = str_split(amino_acids, '/') %>% unlist
      orig = amino_acids[1]
      new = amino_acids[2]
      if (is.na(new)) {
        # splice site
        new = orig
      } else {
        new = sub(paste0('^', orig), '', new)
      }
      cvar = paste0(orig, protein_position, new)
    }
    return(cvar)
  }
}

ivs_to_cdna <- function(ivs_notation) {
  # PAH exon end positions (last nucleotide before each intron)
  exon_ends <- c(
    IVS1 = 117,
    IVS2 = 238,
    IVS3 = 353,
    IVS4 = 442,
    IVS5 = 509,
    IVS6 = 611,
    IVS7 = 707,
    IVS8 = 842,
    IVS9 = 969,
    IVS10 = 1066,
    IVS11 = 1197,
    IVS12 = 1315
  )

  # Parse IVS notation
  # Pattern: IVS<number>[+/-]<offset>[nucleotide change]
  pattern <- "^IVS(\\d+)([+-]\\d+)(.*)$"

  if (!grepl(pattern, ivs_notation)) {
    print("Invalid IVS notation. Expected format: IVS<n>+/-<offset> or IVS<n>+/-<offset><change>")
    print(ivs_notation)
    return(ivs_notation)
    stop("Invalid IVS notation. Expected format: IVS<n>+/-<offset> or IVS<n>+/-<offset><change>")

  }

  # Extract components
  ivs_num <- as.numeric(sub(pattern, "\\1", ivs_notation))
  offset <- sub(pattern, "\\2", ivs_notation)
  change <- sub(pattern, "\\3", ivs_notation)

  # Look up exon end position
  ivs_key <- paste0("IVS", ivs_num)
  if (!ivs_key %in% names(exon_ends)) {
    stop(paste("Invalid IVS number:", ivs_num, ". PAH has 12 introns (IVS1-IVS12)"))
  }

  exon_end <- exon_ends[ivs_key]

  # Construct c. notation
  cdna_notation <- paste0("c.", exon_end, offset, change)

  return(cdna_notation)
}


genomic_to_cdna <- function(introns, pos, ref = NULL, alt = NULL, strand='-') {

  positive = strand == '+'
  if (positive) {
    intronic = introns %>% filter(pos >= genomic_start, pos <= genomic_end)
  } else {
    intronic = introns %>% filter(pos <= genomic_start, pos >= genomic_end)
  }
  #print(introns %>% filter(pos >= genomic_start))
  #print(introns %>% filter(pos <= genomic_end))
  if (nrow(intronic) == 0) {
    return("not_intronic")
  }
  if (nrow(intronic) > 1) {
    return("multi_intronic")
  }
  # genomic_start/genomic_end are the first/last intronic positions,
  # so offset from exon boundary is distance + 1
  dstart = abs(intronic$genomic_start[1] - pos) + 1
  dend = abs(intronic$genomic_end[1] - pos) + 1
  # genomic_start is donor side (cdna_start, +offset), genomic_end is acceptor side (cdna_end, -offset)
  cpos = ifelse(dstart < dend, intronic$cdna_start, intronic$cdna_end)
  relative = ifelse(dstart < dend, paste0('+', dstart), paste0('-', dend))
  if (!positive && !is.null(ref) && !is.null(alt)) {
    ref = chartr("ATCG", "TAGC", toupper(ref))
    alt = chartr("ATCG", "TAGC", toupper(alt))
  }
  cdna = paste0('c.', cpos, relative, ref, '>', alt)
  return(cdna)
}
