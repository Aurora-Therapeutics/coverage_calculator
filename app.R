library(data.table)
library(tidyverse)
library(glue)
library(shiny)
library(jcolors)
# library(ggpattern)

source('~/src/coverage/coordinate_conversion.r')

# --- Data loading (same logic as basic_functionality.R) ---

carriers = 'hetz carriers'
patientx = 'obs. patients'

load_patient_data <- function(gt_file) {

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

  df_pku = fread(gt_file) %>% mutate(study_type=patientx) %>%
    mutate(variant=orig_commonize_variant(aa_codes, variant)) %>%
    select(variant, num, pmid, location, study_type)
  df_pku
}

load_carrier_data <- function(bb_file, annots_file) {

  introns = tribble(~intron, ~cdna_start, ~cdna_end, ~genomic_start, ~genomic_end,
                  1, 60, 61,       102917070, 102912899,
                  2, 168, 169,     102912790, 102894919,
                  3, 352, 353,     102894734, 102877551,
                  4, 441, 442,     102877461, 102866664,
                  5, 509, 510,     102866595, 102855333,
                  6, 706, 707,     102855135, 102852951,
                  7, 842, 843,     102852814, 102851757,
                  8, 912, 913,     102851686, 102846952,
                  9, 969, 970,     102846894, 102844432,
                  10, 1065, 1066,  102844335, 102843780,
                  11, 1199, 1200,  102843645, 102840516,
                  12, 1315, 1316,  102840399, 102839219,
                  )

  df_pku2 = fread(bb_file) %>%  filter(gene == 'PAH')
  population_sizes = tribble(~gnomad_population, ~total_people, ~source, ~pmid,
                             'Admixed American', 30019, 'gnomad',  38057664,
                             'African/African American', 37545, 'gnomad',  38057664,
                             'Ashkenazi Jewish', 14804, 'gnomad',  38057664,
                             'East Asian', 22448, 'gnomad',  38057664,
                             'Non-Finnish European', 590031, 'gnomad',  38057664,
                             'Middle Eastern', 3031, 'gnomad',  38057664,
                             'South Asian', 45546, 'gnomad',  38057664,
                             'Finnish', 32026, 'gnomad',  38057664,
                             'UK BioBank', 981098/2, 'UKBB',  40770095,
                             'TopMed v10', 150899, 'TopMed',  33568819,
                             'All Of Us', 414840, 'AOU',  38374255,
                             'Japan', 61000, 'BBJ', 28189464
  )
  missing_sizes = df_pku2 %>% anti_join(population_sizes, by='gnomad_population') %>% select(gnomad_population) %>% distinct()
  if (nrow(missing_sizes) > 0) {
    print('population dataset sizes are missing:')
    print(missing_sizes)
    stop('fix this')
  }
  df_pku2 = df_pku2 %>%
    inner_join(population_sizes, by='gnomad_population') %>%
    mutate(num=round(AF * total_people * 2)) %>%
    mutate(location=case_when(gnomad_population == 'Japan' & source == 'BBJ' ~ 'Japan: BBJ',
                              gnomad_population == 'All Of Us' & source == 'AOU' ~ 'USA: All Of Us',
                              gnomad_population == 'TopMed v10' & source == 'TopMed' ~ 'USA: TopMedv10',
                              gnomad_population == 'UK BioBank' & source == 'UKBB' ~ 'UK: UKBB',
                              T ~ paste('gnomad:', gnomad_population) ),
           study_type=carriers) %>%
    select(chrm, pos, ref, alt, num, pmid, location, study_type)
  # run once when updating carrier variant data
  if (F) {
    source('~/src/coverage/protein_annotations.R')
    dfx = df_pku2 %>% select(chrm, pos, ref, alt, gene) %>% distinct()
    annots = get_all_variant_annotations(dfx, skip_api=F) %>% mutate
    annots %>% fwrite(annots_file)
  }
  annots = fread(annots_file)
  nc = annots %>% filter(is.na(protein_position)) %>%
    rowwise() %>% mutate(cvar=genomic_to_cdna(introns, pos, ref, alt, strand='-')) %>% ungroup() %>%
    filter(gene == 'PAH') %>%  select(chrm, pos, ref, alt, cvar)
  cod = annots %>% filter(!is.na(protein_position)) %>%
    rowwise() %>% mutate(cvar=commonize_variant_annots(amino_acids, protein_position)) %>% ungroup() %>%
    filter(gene == 'PAH') %>%  select(chrm, pos, ref, alt, cvar)
  annotated = bind_rows(nc, cod)
  difficult = tribble(~chrm, ~pos, ~ref, ~alt, ~cvar,
                        'chr12', 102839172, 'GCTTTA', 'G', 'c.1357_*2del',
                        'chr12', 102840395, 'CTTACTG', 'C', 'c.1314_1315+4del',
                        )
  annotated = bind_rows(nc, cod, difficult)
  missing_annots = df_pku2 %>% anti_join(annotated, by=c('pos', 'ref', 'alt'))
  if (nrow(missing_annots) > 0) {
    print('population dataset sizes are missing:')
    print(missing_sizes)
    stop('fix this')
  }

  df_pku2 = df_pku2 %>% inner_join(annotated, by=c('pos', 'ref', 'alt')) %>%
    rename(variant=cvar) %>%
    select(variant, num, pmid, location, study_type)

}

construct_label <- function(location, pmid, total, label_type) {
  if (label_type == 'country (number of alleles, study PMID)') {
    label = glue("{location} (N={total}, {pmid})")
  } else if (label_type == 'country (number of alleles)') {
    label = glue("{location} (N={total})")
  } else if (label_type == 'country (study PMID)') {
    label = glue("{location} ({pmid})")
  }
  label
}

if (file.exists('~/pah/genotype_tables/')) {
  gt_file = '~/pah/genotype_tables/combined.csv'
  bb_file = '~/disease_genetics/versions/v1.0/pku.csv'
  cpku_file = '~/disease_genetics/versions/v1.0/cpku.csv'
  annots_file = '~/disease_genetics/versions/v1.0/annots.csv'
} else {
  gt_file = '~/data/combined.csv'
  bb_file = '~/data/v1.0/pku.csv'
  cpku_file = '~/data/v1.0/cpku.csv'
  annots_file = '~/data/v1.0/annots.csv'
}

df_pku2 = load_carrier_data(bb_file, annots_file)
df_pku = load_patient_data(gt_file)

df_pku$pmid = as.character(df_pku$pmid)
df_pku2$pmid = as.character(df_pku2$pmid)
df_raw = bind_rows(df_pku, df_pku2)

studies = df_raw %>%
  filter(variant != 'unk', variant != 'Unknown') %>%
  group_by(location, pmid, study_type) %>% summarize(total=sum(num)) %>% ungroup()

variants = df_raw %>%
  filter(variant != 'unk', variant != 'Unknown') %>%
  group_by(location, pmid) %>% mutate(total=sum(num)) %>% ungroup() %>%
  mutate(frac=num/total) %>%
  group_by(variant) %>% summarize(val=sum(frac)) %>%
  arrange(-val)
all_variants = variants$variant

label_types = c('country (number of alleles, study PMID)', 'country (number of alleles)', 'country (study PMID)')
linetypes = c("solid", "dotted")
names(linetypes) = c(carriers, patientx)

# Pre-built variant lists
default_variant_lists = list(
  `Aurora 3` = c('R408W', 'c.1066-11G>A', 'P281L'),
  `Aurora 6` = c('R408W', 'c.1066-11G>A', 'P281L', 'R243X', 'R261Q', 'R111X')
)

# --- UI ---

starting_values = c("Australia | 24368688",
                    "China | 16256386",
                    "France | 26666653",
                    "Iran | 37525467",
                    "Norway | 8831077",
                    "Lithuania | 12655550",
                    "Portugal | 33465300",
                    "Japan: BBJ | 28189464",
                    "UK: UKBB | 40770095",
                    "USA: All Of Us | 38374255",
                    "USA: TopMedv10 | 33568819",
                    "gnomad: Ashkenazi Jewish | 38057664",
                    "gnomad: Non-Finnish European | 38057664"
                    )
starting_studies = c(glue('Australia (N=222, {patientx})'),
                     glue('China (N=349, {patientx})'),
                     glue('France (N=727, {patientx})'),
                     glue('Iran (N=1515, {patientx})'),
                     glue('Japan: BBJ (N=815, {carriers})'),
                     glue('Norway (N=108, {patientx})'),
                     glue('Lithuania (N=175, {patientx})'),
                     glue('Portugal (N=446, {patientx})'),
                     glue('UK: UKBB (N=15022, {carriers})'),
                     glue('USA: All Of Us (N=12240, {carriers})'),
                     glue('USA: TopMedv10 (N=3982, {carriers})'),
                     glue('gnomad: Ashkenazi Jewish (N=1585, {carriers})'),
                     glue('gnomad: Non-Finnish European (N=19213, {carriers})')
                     )
# Build named choices: display name -> internal "location | pmid" value
study_values = paste(studies$location, studies$pmid, sep = " | ")
study_labels = paste0(studies$location, " (N=", studies$total, ", ", studies$study_type, ")")
study_choices = setNames(study_values, study_labels)
# Map starting_studies display names to their internal values
#starting_values = study_choices[starting_studies]
#starting_values = starting_values[!is.na(starting_values)]

ui <- fluidPage(
  titlePanel("PKU Coverage Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Variant Lists"),
      textInput("new_list_name", "New list name", placeholder = "e.g. my_list"),
      selectizeInput("new_list_variants", "Variants for new list",
                     choices = all_variants, multiple = TRUE,
                     options = list(placeholder = "Select variants...")),
      actionButton("add_list", "Add List"),
      hr(),
      h5("Active variant lists:"),
      uiOutput("variant_lists_display"),
      hr(),
      selectizeInput("chosen_studies", "Studies to include",
                     choices = setNames(
                       paste(studies$location, studies$pmid, sep = " | "),
                       paste0(studies$location, " (N=", studies$total, ", ", studies$study_type, ")")
                     ),
                     selected = starting_values,
                     multiple = TRUE,
                     options = list(placeholder = "Click to select studies...")),
      hr(),
      selectInput("label_type", "Label type",
                  choices = label_types, selected = label_types[1]),
      sliderInput("text_size", "Text size",
                  min = 50, max = 200, value = 100, step = 10, post = "%")),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Plot",
                 downloadButton("download_pdf", "Download PDF"),
                 #plotOutput("coverage_plot", height = "600px")
                 uiOutput("coverage_plot_ui")
                 ),
        tabPanel("Data",
                 downloadButton("download_csv", "Download CSV"),
                 br(), br(),
                 tableOutput("coverage_table"))
      )
    )
  )
)

# --- Server ---

server <- function(input, output, session) {

  variant_lists <- reactiveVal(default_variant_lists)

  observeEvent(input$add_list, {
    req(input$new_list_name, input$new_list_variants)
    name <- trimws(input$new_list_name)
    if (nchar(name) == 0) return()
    current <- variant_lists()
    current[[name]] <- input$new_list_variants
    variant_lists(current)
    updateTextInput(session, "new_list_name", value = "")
    updateSelectizeInput(session, "new_list_variants", selected = character(0))
  })

  # Observer for dynamic remove buttons
  observe({
    vl <- variant_lists()
    lapply(names(vl), function(n) {
      btn_id <- paste0("rm_", make.names(n))
      observeEvent(input[[btn_id]], {
        current <- variant_lists()
        current[[n]] <- NULL
        variant_lists(current)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  output$variant_lists_display <- renderUI({
    vl <- variant_lists()
    if (length(vl) == 0) return(em("None - add a list above"))
    tags$div(
      lapply(names(vl), function(n) {
        btn_id <- paste0("rm_", make.names(n))
        tags$div(
          style = "margin-bottom: 6px;",
          actionButton(btn_id, "X", style = "padding: 2px 6px; margin-right: 6px;"),
          strong(n), ": ", paste(vl[[n]], collapse = ", ")
        )
      })
    )
  })

  coverages_data <- reactive({
    req(length(input$chosen_studies) > 0)

    chosen <- do.call(rbind, strsplit(input$chosen_studies, " \\| "))
    chosen_studies <- data.frame(location = chosen[, 1], pmid = chosen[, 2],
                                  stringsAsFactors = FALSE)

    vl <- variant_lists()
    req(length(vl) > 0)
    label_type <- input$label_type

    df <- df_raw %>%
      filter(variant != 'unk', variant != 'Unknown') %>%
      group_by(location, pmid) %>% mutate(total = sum(num)) %>% ungroup() %>%
      mutate(pmid=as.character(pmid)) %>%
      inner_join(chosen_studies, by = c('location', 'pmid'))

    req(nrow(df) > 0)

    coverages <- NULL
    for (lname in names(vl)) {
      vlist <- vl[[lname]]
      covs <- df %>%
        mutate(frac = num / total, chosen = variant %in% vlist) %>%
        group_by(location, pmid, study_type) %>%
        summarize(frac_covered = sum(as.numeric(chosen) * frac),
                  total = max(total), .groups = 'drop') %>%
        mutate(pop_covered = frac_covered**2 + 2 * frac_covered * (1 - frac_covered)) %>%
        rowwise() %>%
        mutate(label = construct_label(location, pmid, total, label_type)) %>%
        ungroup() %>%
        mutate(variant_set = lname, set_size = length(vlist))
      coverages <- bind_rows(coverages, covs)
    }

    levs <- coverages %>%
      group_by(variant_set) %>%
      summarize(v = mean(pop_covered)) %>%
      arrange(-v) %>%
      .$variant_set
    coverages$variant_set <- factor(coverages$variant_set, levels = levs)

    coverages
  })

  make_plot <- function(coverages) {
    base_size <- 18 * input$text_size / 100
    coverages %>%
      arrange(as.numeric(variant_set)) %>%
      ggplot(aes(x = label, y = pop_covered, fill = variant_set, linetype = study_type)) +
      geom_bar(stat = 'identity', alpha = .5, position = 'identity',
               color = 'grey30', linewidth = 0.8,
               key_glyph = "path") +
      theme_minimal(base_size = base_size) +
      coord_flip() +
      scale_fill_jcolors('pal6') +
      scale_linetype_manual(values = linetypes) +
      labs(x = input$label_type,
           y = 'estimated fraction of patients treated by selected variants',
           fill = 'Variant set', linetype = 'Study type')
  }

  output$coverage_plot <- renderPlot({
    make_plot(coverages_data())
  })

  output$coverage_plot_ui <- renderUI({
    n_studies <- length(input$chosen_studies)
    plot_height <- max(400, n_studies * 40)
    plotOutput("coverage_plot", height = paste0(plot_height, "px"))
  })

  output$download_pdf <- downloadHandler(
    filename = function() { "coverage_plot.pdf" },
    content = function(file) {
      p <- make_plot(coverages_data())
      ggsave(file, plot = p, device = "pdf", width = 10, height = 7)
    }
  )

  output$coverage_table <- renderTable({
    coverages_data()
  })

  output$download_csv <- downloadHandler(
    filename = function() { "coverage_data.csv" },
    content = function(file) {
      write.csv(coverages_data(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
