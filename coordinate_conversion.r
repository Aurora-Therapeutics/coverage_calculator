
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

commonize_variant <- function(aa_codes, ivs_conv, amino_acids, protein_position, cdna_position, cdna_end) {
  if (is.na(amino_acids)) {
    istarts = as.numeric(sub('c.', '', ivs_conv$Upstream_Exon_End))
    gt = cdna_position >= istarts
    lt = cdna_position <= istarts
  } else {
    amino_acids = sub('\\*', 'X', amino_acids)
    if ((protein_position == 204) & (amino_acids == 'Y/C')) {
      cvar = 'Ex6-96A>G'
    } else {
      amino_acids = str_split(amino_acids, sep='/') %>% unlist
      orig = amino_acids[1]
      new = amino_acids[2]
      if (is.na(new)) {
        # splice site
        new = old
      } else {
        new = sub(paste0('^', orig), '', new)
      }
      cvar = paste0(orig, protein_position, new)
      return(cvar)
    }
  }

  cvar = ifelse(grepl('^[0-9]', cvar),  paste0('c.', cvar), cvar)

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


genomic_to_cdna_pah <- function(pos, ref = NULL, alt = NULL) {

  # Correct PAH exon boundaries based on your IVS mapping
  exons <- data.frame(
    exon = 1:13,
    cdna_end = c(60, 168, 352, 441, 509, 706, 842, 912, 969, 1065, 1199, 1315, 1356),
    cdna_start = c(1, 61, 169, 353, 442, 510, 707, 843, 913, 970, 1066, 1200, 1316),
    # Genomic coordinates (hg38, chr12, minus strand)
    start = c(102958335, 102946866, 102945242, 102939792, 102936175,
              102925359, 102914637, 102903826, 102894474, 102891850,
              102883932, 102879321, 102836982),
    end = c(102958422, 102946986, 102945395, 102939901, 102936241,
            102925470, 102914778, 102904000, 102894570, 102891947,
            102884063, 102879438, 102837108)
  )

  # Calculate genomic size of each exon
  exons$genomic_size <- exons$end - exons$start + 1
  exons$cdna_size <- exons$cdna_end - exons$cdna_start + 1

  # Find which exon (if any)
  in_exon <- which(pos >= exons$start & pos <= exons$end)

  if (length(in_exon) > 0) {
    return(list(cdna = NA, ivs = NA, location = "exonic"))
  }

  # Calculate distances to all splice boundaries
  dist_to_ends <- abs(pos - exons$end)
  dist_to_starts <- abs(exons$start)

  # Find closest boundary
  closest_end <- which.min(dist_to_ends)
  closest_start <- which.min(dist_to_starts)

  if (dist_to_ends[closest_end] < dist_to_starts[closest_start]) {
    # Closest to exon end (donor, + notation)
    intron <- closest_end
    cdna_base <- exons$cdna_end[intron]
    offset <- pos - exons$end[intron]
    offset_str <- paste0("+", offset)
  } else {
    # Closest to exon start (acceptor, - notation)
    intron <- closest_start - 1
    cdna_base <- exons$cdna_start[closest_start]
    offset <- exons$start[closest_start] - pos
    offset_str <- paste0("-", offset)
  }

  # Variant
  if (!is.null(ref) && !is.null(alt)) {
    ref_comp <- chartr("ATCG", "TAGC", toupper(ref))
    alt_comp <- chartr("ATCG", "TAGC", toupper(alt))
    cdna <- paste0("c.", cdna_base, offset_str, ref_comp, ">", alt_comp)
    ivs <- paste0("IVS", intron, offset_str, ref_comp, ">", alt_comp)
  } else {
    cdna <- paste0("c.", cdna_base, offset_str)
    ivs <- paste0("IVS", intron, offset_str)
  }

  return(list(cdna = cdna, ivs = ivs, intron = intron))
}
