library(httr)
library(jsonlite)

get_end_pos <- function(pos, ref, alt) {
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)

  if (ref_len == alt_len) {
    # SNV or MNV
    end_pos <- as.numeric(pos) + ref_len - 1
  } else if (ref_len > alt_len) {
    # Deletion
    end_pos <- as.numeric(pos) + ref_len - 1
  } else {
    # Insertion
    end_pos <- as.numeric(pos)
  }
  end_pos
}

query_vep_rest <- function(chr, pos, end_pos, alt, ref) {
  variant_string <- paste0(chr, ":", pos, "-", end_pos, "/", alt, "/", ref)
  url <- paste0("https://rest.ensembl.org/vep/human/region/",
                variant_string,
                "?canonical=1&mane_select=1")
  response <- GET(url,
                  add_headers(`Content-Type` = "application/json"),
                  timeout(30))

  if (status_code(response) == 200) {
    result <- fromJSON(content(response, "text", encoding='UTF-8'))
  } else {
    result = NULL
  }
  result
}

get_best_transcript <- function(tc) {
  # tc is the transcript_consequences data frame

  # 1. MANE Select (gold standard - one per gene)
  mane_select <- tc[!is.na(tc$mane_select), ]
  if (nrow(mane_select) > 0)
    return(mane_select[1, ])

  # 2. MANE Plus Clinical (clinical significance)
  mane_plus <- tc[!is.na(tc$mane_plus_clinical), ]
  if (nrow(mane_plus) > 0)
    return(mane_plus[1, ])

  # 3. Canonical transcript (Ensembl's representative)
  canonical <- tc[!is.na(tc$canonical) & tc$canonical == 1, ]
  if (nrow(canonical) > 0)
    return(canonical[1, ])

  # 4. Protein coding with amino acid change
  protein_coding <- tc[tc$biotype == "protein_coding" & !is.na(tc$amino_acids), ]
  if (nrow(protein_coding) > 0)
    return(protein_coding[1, ])

  # 5. Any protein coding
  any_protein <- tc[tc$biotype == "protein_coding", ]
  if (nrow(any_protein) > 0)
    return(any_protein[1, ])

  # 6. Fallback to first transcript
  return(tc[1, ])
}

get_protein_effect_vep <- function(chr, pos, ref, alt, gene) {
  # Calculate end position for deletions
  empty = data.frame(hgvsp = NA, consequence = NA, gene_symbol = NA)
  end_pos = get_end_pos(pos, ref, alt)
  #variant_string <- paste0(chr, ":", pos, "-", end_pos, "/", ref, "/", alt)
  result = query_vep_rest(chr, pos, end_pos, alt, ref)
  if (is.null(result))
    return(empty)
  if (length(result) == 0)
    return(empty)
  if (is.null(result$transcript_consequences))
    return(empty)
  tc <- result$transcript_consequences[[1]]
  tc = tc[tc$gene_symbol == gene,]
  #print(tc)
  if (nrow(tc) == 0)
    return(empty)
  best = get_best_transcript(tc)
  return(data.frame(
    gene_symbol = ifelse(!is.null(best$gene_symbol), best$gene_symbol, NA),
    amino_acids = ifelse(!is.null(best$amino_acids), best$amino_acids, NA),
    protein_position = ifelse(!is.null(best$protein_start), best$protein_start, NA),
    cdna_position = ifelse(!is.null(best$cdna_start), best$cdna_start, NA),
    cdna_end = ifelse(!is.null(best$cdna_end), best$cdna_end, NA),
    consequence = paste(best$consequence_terms[[1]], collapse = ",")
    ))
}

get_all_variant_annotations <- function(variants, skip_api=F) {
  variants = variants %>% select(gene, chrm, pos, ref, alt) %>% distinct()
  if (skip_api) {
    variants_annotated = variants %>% mutate(label=paste0(gene, ':', chrm, ':', pos, ':', ref, ':', alt))
    variants_annotated$label = factor(variants_annotated$label, levels=variants_annotated$label)
  } else {
    annotations <- lapply(1:nrow(variants), function(i) {
      Sys.sleep(0.15)
      get_protein_effect_vep(sub('chr', '', variants$chrm[i]), variants$pos[i],
                         variants$ref[i], variants$alt[i], variants$gene[i])
    })
    variants_annotated <- cbind(variants, bind_rows(annotations))
    variants_annotated = variants_annotated %>%
      #mutate(label=ifelse(is.na(amino_acids), paste0(allele_num, ", ", chrm, ':', pos, ':', ref, ':', alt),
      #                    paste0(allele_num, ', ', gene_symbol, ':', protein_position, ':', amino_acids))) %>%
      mutate(label=ifelse(is.na(amino_acids), paste0(gene, ':', chrm, ':', pos, ':', ref, ':', alt),
                          paste0(gene_symbol, ':', protein_position, ':', amino_acids)))
    variants_annotated$label = factor(variants_annotated$label, levels=unique(variants_annotated$label))
  }
  variants_annotated
}
