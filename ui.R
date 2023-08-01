popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'

shinyUI(
  fluidPage(
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
        color: #ff0000;
        font-weight: bold;
      }
    "))
    ),
    useShinyjs(),
    extendShinyjs(text = jscode, functions = c("closeWindow")),
    actionButton("close", "",icon = icon("window-close")),
    tags$a(href="javascript:history.go(0)", 
           popify(tags$i(class="fa fa-refresh"),
                  title = "Reload", 
                  content = "Click here to restart the session",
                  placement = "right")),
    tags$style(
      tags$head(
        includeScript("navAppend.js")
      ),
      tags$head(
        includeScript("navAppend2.js")
      ),
      type = 'text/css',
      HTML(
        ".container-fluid > .nav > li >
                        a[data-value='Title'] {font-size: 20px}",
        ".navbar{margin:3px;}",
        ".navbar-default {background-color: #EAF3EE !important;}",
        ".well { background-color: #EAF3EE !important;}"
      )
    ),
    tags$style(".popover{
            max-width: 500px;
          }"),

    navbarPage(
      footer=p(hr(),p("Need help? Create an issue on", a("Github", href = "https://github.com/Kan-E/RNAseqChef/issues"), 
                      "or", a("contact us", href = "maito:kaneto@kumamoto-u.ac.jp"),".",align="center",width=4)
      ),
      "",
      id='navBar',
      tabPanel("EpigenomeChef" ,value='Title', icon = icon("utensils"),
               fluidRow(
                 column(12,
                        br(),br(),
                        h1(strong("EpigenomeChef"),align="center"),br(),
                        p("EpigenomeChef is a web-based application for highlighting the differential epigenome-transcriptome interaction.",
                          align="center"),br(),br(),style={'background-color:mintcream;font-size: 16px;'},
                 ),
                 column(12,br(),
                        h3("Pair-wise DAR detects and visualizes differentially accessible regions"),br(),
                        img(src="pair-wise_DAR.png", width = 1200,height = 1600),br(),br(),hr(),
                        h3("Venn diagram detects and visualizes the overlap between DARs from multiple datasets"),br(),
                        img(src="Venn.png", width = 1200,height = 700),br(),br(),hr(), 
                        h3("Clustering identifies similar samples and DNA binding patterns by clustering methods"),br(),
                        img(src="Clustering.png", width = 1000,height = 900),br(),hr(),
                        img(src="combined.png", width = 1200,height = 600),br(),hr(),
                 )
               )
      ),
      # pair-wise -------------------------------------
      tabPanel("Pair-wise DAR",
               sidebarLayout(
                 # sidebar---------------------------------
                 sidebarPanel(
                   radioButtons('data_file_type','Input:',
                                c('BigWig files'="Row1",
                                  'Bam files'="Row2"
                                ),selected = "Row1"),
                   conditionalPanel(condition="input.data_file_type=='Row1'",
                                    fileInput("file1",
                                              strong(
                                                span("Select BigWig files"),
                                                span(icon("info-circle"), id = "icon1", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("bw","BigWig"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon1", "BigWig files (bw, BigWig):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                                            strong("Do not use it for anything else"),".<br>Short names are recommended for file names."), 
                                              placement = "right",options = list(container = "body"))
                                    
                   ),
                   conditionalPanel(condition="input.data_file_type=='Row2'",
                                    fileInput("file_bam",
                                              strong(
                                                span("Select Bam files"),
                                                span(icon("info-circle"), id = "icon_bam", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("bam"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon1", "Bam files (bam):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                                            strong("Do not use it for anything else"),".<br>Short names are recommended for file names.<br>"), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type=='Row2'",
                                    radioButtons('Pair_or_single','Sequence type:',
                                                 c('Paired-end'="Paired-end",
                                                   'Single-end'="Single-end"),selected = "Paired-end"),         
                   ),
                   radioButtons('Genomic_region','Genomic region:',
                                c('Genome-wide'="Genome-wide",
                                  'Promoter'="Promoter"),selected = "Genome-wide"),
                   conditionalPanel(condition="input.Genomic_region=='Promoter'",
                                    fluidRow(
                                      column(5, numericInput("upstream", "upstream", value = 500, min = 0)),
                                      column(5, numericInput("downstream", "downstream", value = 500, min = 0))
                                    )
                   ),
                   conditionalPanel(condition="input.Genomic_region=='Genome-wide'",
                                    fileInput("peak_call_file1",
                                              strong(
                                                span("Select peak call files"),
                                                span(icon("info-circle"), id = "icon2", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("bed","narrowPeak"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    
                                    bsPopover("icon2", "peak call files (bed, narrowPeak):", 
                                              content=paste("The first column is chromosome (e.g. chr3, chrY) name.<br>",
                                                            "The second column is start position on the chromosome.<br>",
                                                            "The third column is end position on the chromosome.<br>",
                                                            strong("Multiple peak call files are required for the filteration.<br>"),
                                                            img(src="pair_peakcallfilter.png", width = 200,height = 400)),
                                              placement = "right",options = list(container = "body"))
                   ),
                   fluidRow(
                     column(12, selectInput("Species", "Species", species_list, selected = "not selected"))),
                   fluidRow(
                     column(4, numericInput("fc", "Fold Change", min   = 1, max   = NA, value = 2)),
                     column(4, numericInput("fdr", "FDR", min   = 0, max   = 1, value = 0.1, step = 0.005)),
                     column(4, numericInput("basemean", "Basemean", min   = 0, max   = NA, value = 0))
                   ),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "pair_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("pair_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("pair_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("pair_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Clustering:"), "height = 3.5, width = 9<br>", 
                                           strong("Volcano-plot:"), "height = 4, width = 4 <br>",
                                           strong("Heatmap:"), "height = 4, width = 3 <br>",
                                           strong("Trackplot:"), "height = 4, width = 7 <br>",
                                           strong("Peak distribution:"), "height = 6, width = 10 <br>",
                                           strong("Peak pattern heatmap:"), "height = 6, width = 6 <br>",
                                           strong("Peak pattern line plot:"), "height = 5, width = 5 <br>",
                                           strong("Dotplot (GREAT):"), "height = 5, width = 6 <br>",
                                           strong("Region gene associations plot (GREAT):"), "height = 5, width = 12 <br>",
                                           strong("Motif dot plot:"), "height = 6, width = 7 <br>",
                                           strong("Boxplot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("Barplot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("KS-plot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("Dotplot (with RNAseq):"), "height = 5, width = 8 <br>",
                                           strong("combined heatmap:"), "height = 6, width = 10 <br>"),trigger = "click"), 
                   actionButton("goButton", "example data (hg19)"),
                   tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             )),
                   fluidRow(column(7),
                   column(3, downloadButton("pair_report", "Download summary"),
                          tags$head(tags$style("#pair_report{color: red;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                                    tags$style("
          body {
            padding: 0 !important;
          }"
                                    )))
                   )
                 ), #sidebarPanel
                 
                 # Main Panel -------------------------------------
                 mainPanel(
                   tabsetPanel(
                     type = "tabs",
                     tabPanel("Input Data",
                              bsCollapse(id="input_collapse_panel",open="bw_files_panel",multiple = TRUE,
                                         bsCollapsePanel(title="uploaded bigwig or bam files:",
                                                         value="bw_files_panel",
                                                         dataTableOutput('input_bw_files') 
                                         ),
                                         bsCollapsePanel(title="uploaded peak call files:",
                                                         value="peak_call_files_panel",
                                                         downloadButton("download_filtered_peakcall_bed", "Download filtered_marged (.bed)"),
                                                         dataTableOutput('input_peak_call_files')
                                         ),
                                         bsCollapsePanel(title="Count data:",
                                                         value="raw_count_panel",
                                                         column(4, textOutput("Spe"),
                                                                tags$head(tags$style("#Spe{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                                         downloadButton("download_raw_count_table", "Download count table"),
                                                         dataTableOutput('raw_count_table')
                                         )
                              )
                     ),
                     tabPanel("Result overview",
                              fluidRow(
                                column(4, downloadButton("download_pair_PCA", "Download clustering analysis"))
                              ),
                              plotOutput("PCA"),
                              fluidRow(
                                column(4, downloadButton("download_pair_volcano", "Download volcano plot"))
                              ),
                              fluidRow(
                                column(6, htmlOutput("volcano_x")),
                                column(6, htmlOutput("volcano_y"))),
                              plotOutput("volcano1",brush = "plot1_brush"),
                              bsCollapse(id="input_collapse_pair_DEG",open="DEG_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="Trackplot_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_trackplot", "Download trackplot")),
                                                           column(8, htmlOutput("trackplot_additional"))
                                                         ),
                                                         fluidRow(
                                                           column(8, htmlOutput("igv_uprange")),
                                                           column(4, htmlOutput("igv_ylim"))),
                                                         column(4, textOutput("Spe_track"),
                                                                tags$head(tags$style("#Spe_track{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                                         plotOutput("trackplot_goi")
                                         ),
                                         bsCollapsePanel(title="DAR result:",
                                                         value="DEG_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_DEG_result", "Download DEG result"))
                                                         ),
                                                         DTOutput("DEG_result")
                                         ),
                                         bsCollapsePanel(title="Up_peaks:",
                                                         value="DEG_up_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_DEG_up_result", "Download Up regions (.bed)"))
                                                         ),
                                                         DTOutput("DEG_up_result")
                                         ),
                                         bsCollapsePanel(title="Down_peaks:",
                                                         value="DEG_down_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_DEG_down_result", "Download Down regions (.bed)"))
                                                         ),
                                                         DTOutput("DEG_down_result")
                                         ),
                                         bsCollapsePanel(title="Normalized_Count_matrix:",
                                                         value="norm_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_norm_count", "Download normalized count matrix"))
                                                         ),
                                                         dataTableOutput("Normalized_Count_matrix")
                                         ),
                                         bsCollapsePanel(title="PCA:",
                                                         value="PCA_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_PCA_table", "Download PCA table"))
                                                         ),
                                                         dataTableOutput("pair_PCA_data")
                                         )
                              )),
                     tabPanel("Peak distribution",
                              downloadButton("download_input_peak_distribution", "Download Up DAR distribution"),
                              column(4, textOutput("Spe_dist_promoter"),
                                     tags$head(tags$style("#Spe_dist_promoter{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                              textOutput("Spe_dist"),
                                     tags$head(tags$style("#Spe_dist{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("input_peak_distribution"),
                              downloadButton("download_deg_peak_distribution", "Download Down DAR distribution"),
                              plotOutput("deg_peak_distribution"),
                     ),
                     tabPanel("Peak pattern",
                              fluidRow(
                                column(4, downloadButton("download_peak_pattern_up_heatmap", "Download peak heatmap"),
                                       htmlOutput("peak_pattern_up_bw"),
                                       htmlOutput("peak_pattern_up_heat_range")),
                                column(4, downloadButton("download_peak_pattern_up_line", "Download peak aligned distribution"),
                                       htmlOutput("peak_pattern_up_additional")),
                              ),
                              fluidRow(
                                column(6,  plotOutput("peak_pattern_up_heatmap")),
                                column(6,  plotOutput("peak_pattern_up_line"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_peak_pattern_down_heatmap", "Download peak heatmap"),
                                       htmlOutput("peak_pattern_down_heat_range")),
                                column(4, downloadButton("download_peak_pattern_down_line", "Download peak aligned distribution"))
                              ),
                              fluidRow(
                                column(6,  plotOutput("peak_pattern_down_heatmap")),
                                column(6,  plotOutput("peak_pattern_down_line"))
                              ),
                     ),
                     tabPanel("GREAT",
                              column(4, textOutput("Spe_great_promoter"),
                                     tags$head(tags$style("#Spe_great_promoter{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                              fluidRow(
                                column(4, textOutput("Spe_GREAT"),
                                       tags$head(tags$style("#Spe_GREAT{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(4, htmlOutput("Gene_set")),
                                column(4, downloadButton("download_pair_enrichment", "Download"))
                              ),
                              plotOutput("enrichment1"),
                              bsCollapse(id="input_collapse_pair_enrich",open="ORA_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Enrichment result:",
                                                         value="ORA_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_enrichment_table", "Download enrichment result"))
                                                         ),
                                                         dataTableOutput('pair_enrichment_result')
                                         )
                              ),
                              fluidRow(
                                column(4, htmlOutput("whichGenes"), 
                                       htmlOutput("whichGeneSet")),
                                column(4, downloadButton("download_region_gene_associations_plot", "Download plot"))
                              ),
                              plotOutput('region_gene_associations_plot'),
                              fluidRow(
                                column(4, downloadButton("download_region_gene_associations", "Download table data"))
                              ),
                              DTOutput('region_gene_associations')
                     ),
                     tabPanel("HOMER",
                              fluidRow(
                                column(3, downloadButton("download_motif_plot", "Download motif plot")),
                                column(3, downloadButton("download_homer_report", "Download homer report"),
                                       tags$head(tags$style("#download_homer_report{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),)
                              ),
                              textOutput("Spe_motif"),
                                     tags$head(tags$style("#Spe_motif{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              htmlOutput("homer_unknown"),
                              fluidRow(
                                column(4, htmlOutput("homer_size")),
                                column(4, htmlOutput("homer_size2"))
                              ),
                              htmlOutput("homer_bg"),
                              fluidRow(
                                column(4, actionButton("motifButton", "Start"),
                                       tags$head(tags$style("#motifButton{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                )
                              ),
                              htmlOutput("homer_showCategory"),
                              plotOutput("motif_plot"),
                              bsCollapse(id="Promoter_motif_collapse_panel",open="motif_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Known motif",
                                                         value="motif_result_table",
                                                         downloadButton("download_motif_table", "Download motif enrichment result"),
                                                         DTOutput('motif_result')
                                         ),
                                         bsCollapsePanel(title= "de_novo_motif",
                                                         value="Promoter_motif_region_panel",
                                                         downloadButton("download_denovo_motif_table", "Download de novo motif"),
                                                         dataTableOutput("denovo_motif_result")
                                         )
                              )
                     ),
                     tabPanel("with RNA-seq",
                              bsCollapse(id="RNAseqresult_pair_panel",open="RNAseqresult_pair_table",
                              bsCollapsePanel(title= p(span("RNA-seq result"),span(icon("info-circle"), id = "RNAseqresult_pair", 
                                                                               options = list(template = popoverTempate))),
                                              bsPopover("RNAseqresult_pair", "RNA-seq result file:", 
                                                        content=paste("You can use a pair-wise RNA-seq DEG result file as input.<br>",
                                                                      "First column must be gene name (Gene symbol or ENSEMBL ID).<br>", 
                                                                      "The file must contain", strong("log2FoldChange"), "and", strong("padj"), "columns.<br>",
                                                                      "<br>", 
                                                                      img(src="input_format_volcano.png", width = 500,height = 230)), 
                                                        placement = "right",options = list(container = "body")),
                                              value="RNAseqresult_pair_table",
                                              fluidRow(
                                                column(8, htmlOutput("pairRNAseqresult"))
                                              ),
                                              dataTableOutput('pair_DEG_result')
                              )),
                              htmlOutput("peak_distance"),
                              fluidRow(
                                column(4, downloadButton("download_pairintbox","Download box plot"))
                              ),
                              textOutput("Spe_int"),
                              tags$head(tags$style("#Spe_int{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              textOutput("RNAseq_boxplot_error"),
                              tags$head(tags$style("#RNAseq_boxplot_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("int_boxplot"),
                              fluidRow(
                                column(4, htmlOutput("DEG_fc")),
                                column(4, htmlOutput("DEG_fdr"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_pairintbar","Download bar plot")),
                                column(4, downloadButton("download_pairKSplot","Download ks plot"))
                              ),
                              fluidRow(
                                column(6,plotOutput("int_bar")),
                                column(6,plotOutput('ks_plot'))
                              ),
                              bsCollapse(id="input_collapse_pair_RP",open="RP_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="int_Trackplot_panel",
                                                         fluidRow(
                                                           column(6, downloadButton("download_pair_int_trackplot", "Download trackplot")),
                                                           column(6, htmlOutput("int_trackplot_additional"))
                                                         ),
                                                         fluidRow(
                                                           column(8, htmlOutput("int_igv_uprange")),
                                                           column(4, htmlOutput("int_igv_ylim"))),
                                                         plotOutput("int_trackplot_goi")
                                         ),
                                         bsCollapsePanel(title="Result table:",
                                                         value="RP_panel",
                                                         htmlOutput("RNAseqGroup"),
                                                         htmlOutput("ChIPseqGroup"),
                                                         fluidRow(
                                                           column(4, downloadButton("download_RP_table","Download summary table")),
                                                           column(4, downloadButton("download_selected_RP_table","Download selected table"))
                                                         ),
                                                         fluidRow(
                                                           column(6, downloadButton("download_selected_int_peak","Download selected peak file (bed)"))
                                                         ),
                                                         DTOutput('RP_table')
                                         ),
                                         bsCollapsePanel(title="Functional enrichment analysis:",
                                                         value="function_panel",
                                                         fluidRow(
                                                           column(6, htmlOutput("intGroup")),
                                                           column(6, htmlOutput("intGeneset"))
                                                         ),
                                                         downloadButton("download_pair_int_enrichment", "Download"),
                                                         plotOutput('int_enrichment1'),
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_int_enrichment_table", "Download enrichment result"))
                                                         ),
                                                         dataTableOutput("int_enrichment_result")
                                         )
                                         
                              )
                     ),
                     tabPanel("Combined heatmap",
                              textOutput("Spe_rnaseq2"),
                                     tags$head(tags$style("#Spe_rnaseq2{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(6, htmlOutput("rnaseq_DEGs")),
                                column(6, htmlOutput("rnaseq_count"))
                              ),
                              bsCollapse(id="z-score_count",open="z-score_multiple_count_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Uploaded DEG result files",
                                                         value="Uploaded_DEGs",
                                                         dataTableOutput('rnaseq_DEGs_output')
                                         ),
                                         bsCollapsePanel(title="z-score of uploaded count files",
                                                         value="z-score_multiple_count_panel",
                                                         dataTableOutput('rnaseq_count_output')
                                         )
                              ),
                              fluidRow(
                                column(4, htmlOutput("integrated_bw1")),
                                column(4, htmlOutput("integrated_bw2")),
                                column(4, htmlOutput("integrated_bw3"))
                              ),
                              fluidRow(
                                column(4, actionButton("integrated_heatmapButton", "Start"),
                                       tags$head(tags$style("#integrated_heatmapButton{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                ),
                                column(4, downloadButton("download_integrated_heatmap","Download heatmap"))
                              ),
                              div(
                                plotOutput('integrated_heatmap', height = "100%"),
                                style = "height: calc(75vh  - 75px)"
                              ),
                     )
                   )
                 ) # main panel
               ) #sidebarLayout
      ), #tabPanel
      tabPanel("Venn diagram",
               sidebarLayout(
                 # venn diagram analysis---------------------------------
                 sidebarPanel(
                   fileInput("peak_call_file_venn1",
                             strong(
                               span("Select bed files"),
                               span(icon("info-circle"), id = "icon_venn1", 
                                    options = list(template = popoverTempate))
                             ),
                             accept = c("bed","narrowPeak"),
                             multiple = TRUE,
                             width = "80%"),
                   bsPopover("icon_venn1", "Bed files (bed, narrowPeak):", 
                             content=paste(strong("The maximum number of uploads is five."),
                                           "The first column is chromosome (e.g. chr3, chrY) name.<br>",
                                           "The second column is start position on the chromosome.<br>",
                                           "The third column is end position on the chromosome.<br>",
                                           strong("Short names are recommended for file names.")), 
                             placement = "right",options = list(container = "body")),
                   fileInput("file_venn1",
                             strong(
                               span("Select BigWig files"),
                               span(icon("info-circle"), id = "icon_venn2", 
                                    options = list(template = popoverTempate))
                             ),
                             accept = c("bw","BigWig"),
                             multiple = TRUE,
                             width = "80%"),
                   bsPopover("icon_venn2", "BigWig files (bw, BigWig):", 
                             content=paste(strong("Short names are recommended for file names.")), 
                             placement = "right",options = list(container = "body")),
                   fluidRow(
                     column(12, selectInput("Species_venn", "Species", species_list, selected = "not selected")),
                   ),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "venn_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("venn_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("venn_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("venn_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Venn diagram:"), "height = 3, width = 3<br>", 
                                           strong("Peak distribution:"), "height = 4.5, width = 6 <br>",
                                           strong("Peak pattern heatmap:"), "height = 6, width = 6 <br>",
                                           strong("Peak pattern line plot:"), "height = 5, width = 5 <br>",
                                           strong("Track plot:"), "height = 4, width = 7 <br>",
                                           strong("Dotplot (GREAT):"), "height = 6, width = 8 <br>",
                                           strong("Region gene associations plot (GREAT):"), "height = 5, width = 12 <br>",
                                           strong("Motif Dotplot:"), "height = 6, width = 8 <br>",
                                           strong("Barplot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("Boxplot (with RNAseq):"), "height = 5, width = 6 <br>",
                                           strong("KS-plot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("Dotplot (with RNAseq):"), "height = 5, width = 8 <br>",
                                           strong("combined heatmap:"), "height = 6, width = 10 <br>"),trigger = "click"), 
                   actionButton("goButton_venn", "example data (hg19)"),
                   tags$head(tags$style("#goButton_venn{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             )
                   ),
                   fluidRow(column(7),
                            column(3, downloadButton("venn_report", "Download summary"),
                                   tags$head(tags$style("#venn_report{color: red;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                                             tags$style("
          body {
            padding: 0 !important;
          }"
                                             )))
                   )#sidebarPanel
                 ),
                 
                 # Main Panel -------------------------------------
                 mainPanel(
                   tabsetPanel(
                     type = "tabs",
                     tabPanel("Venn diagram",
                              fluidRow(
                                column(4, downloadButton("download_vennplot", "Download venn diagram"))
                              ),
                              plotOutput("venn"),
                              bsCollapse(id="Lineplot of intersections",multiple = TRUE,
                                         bsCollapsePanel(title="Lineplot of intersections (comparison between intersections):",
                                                         value="vennLineplot_intersection_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_venn_intersection_line", "Download line plot"))
                                                         ),
                                                         plotOutput("venn_batch_bed_line")
                                         ),
                                         bsCollapsePanel(title="Lineplot of intersections (comparison between bigwig):",
                                                         value="vennLineplot_bigwig_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_venn_bigwig_line", "Download line plot"))
                                                         ),
                                                         plotOutput("venn_batch_bigwig_line")
                                         )
                              )
                     ),
                     tabPanel("Peak pattern",
                              htmlOutput("select_file2"),
                              fluidRow(
                                column(4, downloadButton("download_peak_pattern_venn_heatmap", "Download peak heatmap"),
                                       htmlOutput("peak_pattern_venn_heat_range")),
                                column(4, downloadButton("download_peak_pattern_venn_line", "Download peak aligned distribution")
                              )),
                              fluidRow(
                                column(6,  plotOutput("peak_pattern_venn_heatmap")),
                                column(6,  plotOutput("peak_pattern_venn_line"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_venn_peak_distribution", "download peak distribution"))
                              ),
                              textOutput("Spe_venn_distribution"),
                              tags$head(tags$style("#Spe_venn_distribution{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("venn_peak_distribution"),
                              bsCollapse(id="int_result_collapse_panel",open="selected_intersect_annotation_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="Trackplot_venn_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_venn_trackplot", "Download trackplot"))
                                                         ),
                                                         fluidRow(
                                                           column(8, htmlOutput("igv_venn_uprange")),
                                                           column(4, htmlOutput("igv_venn_ylim"))),
                                                         plotOutput("trackplot_venn_goi")
                                         ),
                                         bsCollapsePanel(title="selected intersect data:",
                                                         value="selected_intersect_annotation_panel",
                                                         downloadButton("download_selected_intersect_annotation_table", "Download selected intersection table (.txt)"),
                                                         downloadButton("download_selected_intersect_annotation_table_bed", "Download selected intersection table (.bed)"),
                                                         DTOutput("selected_intersect_annotation")
                                         )
                              )
                     ),
                     tabPanel("GREAT",
                              fluidRow(
                                column(4, htmlOutput("venn_whichGroup2")),
                                column(4, htmlOutput("Gene_set_venn")),
                                column(4, downloadButton("download_venn_enrichment", "Download"))
                              ),
                              plotOutput("enrichment_venn"),
                              bsCollapse(id="input_collapse_venn_enrich",open="ORA_venn_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Enrichment result:",
                                                         value="ORA_venn_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_venn_enrichment_table", "Download enrichment result"))
                                                         ),
                                                         dataTableOutput('venn_enrichment_result')
                                         )
                              ),
                              fluidRow(
                                column(6, htmlOutput("whichGenes_venn"), 
                                       htmlOutput("whichGeneSet_venn")),
                                column(4, downloadButton("download_venn_region_gene_associations_plot", "Download gene associations plot"))
                              ),
                              plotOutput('region_gene_associations_venn_plot'),
                              fluidRow(
                                column(4, downloadButton("download_venn_region_gene_associations", "Download table data"))
                              ),
                              DTOutput('region_gene_venn_associations')
                     ),
                     tabPanel("HOMER",
                              fluidRow(
                                column(3, downloadButton("download_motif_venn_plot", "Download motif plot")),
                                column(3, downloadButton("download_homer_report_venn", "Download homer report"),
                                       tags$head(tags$style("#download_homer_report_venn{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),)
                              ),
                              fluidRow(
                                column(8, htmlOutput("venn_whichGroup1"))
                              ),
                              textOutput("Spe_motif_venn"),
                              tags$head(tags$style("#Spe_motif_venn{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, htmlOutput("homer_unknown_venn")),
                                column(4, htmlOutput("homer_size_venn")),
                                column(4, htmlOutput("homer_size2_venn"))
                              ),
                              fluidRow(
                                column(4, htmlOutput("homer_bg_venn")),
                                column(4, htmlOutput("homer_bg2_venn"))
                              ),
                              fluidRow(
                                column(4, actionButton("motifButton_venn", "Start"),
                                       tags$head(tags$style("#motifButton_venn{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                )
                              ),
                              
                              htmlOutput("homer_showCategory_venn"),
                              plotOutput("motif_venn_plot"),
                              bsCollapse(id="Promoter_motif_venn_collapse_panel",open="motif_venn_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Known motif",
                                                         value="motif_venn_result_table",
                                                         downloadButton("download_motif_venn_table", "Download known motif result"),
                                                         DTOutput('motif_venn_result')
                                         ),
                                         bsCollapsePanel(title= "de_novo_motif",
                                                         value="Promoter_motif_region_panel",
                                                         downloadButton("download_denovo_motif_venn_table", "Download motif region"),
                                                         dataTableOutput("denovo_motif_venn_result")
                                         )
                              )
                     ),
                     tabPanel("with RNA-seq",
                              bsCollapse(id="RNAseqresult_venn_panel",open="RNAseqresult_venn_table",
                                         bsCollapsePanel(title= p(span("RNA-seq result"),span(icon("info-circle"), id = "RNAseqresult_venn", 
                                                                                              options = list(template = popoverTempate))),
                                                         bsPopover("RNAseqresult_venn", "RNA-seq result file:", 
                                                                   content=paste("You can use a pair-wise RNA-seq DEG result file as input.<br>",
                                                                                 "First column must be gene name (Gene symbol or ENSEMBL ID).<br>", 
                                                                                 "The file must contain", strong("log2FoldChange"), "and", strong("padj"), "columns.<br>",
                                                                                 "<br>", 
                                                                                 img(src="input_format_volcano.png", width = 500,height = 230)), 
                                                                   placement = "right",options = list(container = "body")),
                                                         value="RNAseqresult_venn_table",
                              fluidRow(
                                column(8, htmlOutput("vennRNAseqresult"))
                              ),
                              dataTableOutput('venn_DEG_result')
                                         )),
                              htmlOutput("venn_select_RNA"),
                              tags$head(tags$style("#venn_select_RNA{color: black;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, htmlOutput("peak_distance_venn")),
                                column(4, downloadButton("download_vennintbox","Download box plot"))
                              ),
                              textOutput("Spe_int_venn"),
                              tags$head(tags$style("#Spe_int_venn{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              textOutput("RNAseq_boxplot_venn_error"),
                              tags$head(tags$style("#RNAseq_boxplot_venn_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("int_box_venn"),
                              fluidRow(
                                column(4, htmlOutput("DEG_fc_venn")),
                                column(4, htmlOutput("DEG_fdr_venn"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_vennintbar","Download bar plot")),
                                column(4, downloadButton("download_vennKSplot","Download ks plot"))
                              ),
                              fluidRow(
                                column(6,plotOutput("int_bar_venn")),
                                column(6,plotOutput('ks_plot_venn'))
                              ),
                              bsCollapse(id="input_collapse_venn_RP",open="RP_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="int_Trackplot_panel",
                                                         fluidRow(
                                                           column(6, downloadButton("download_venn_int_trackplot", "Download trackplot")),
                                                           column(6, htmlOutput("int_trackplot_additional_venn"))
                                                         ),
                                                         fluidRow(
                                                           column(8, htmlOutput("int_igv_uprange_venn")),
                                                           column(4, htmlOutput("int_igv_ylim_venn"))),
                                                         plotOutput("int_trackplot_goi_venn")
                                         ),
                                         bsCollapsePanel(title="Result table:",
                                                         value="RP_panel",
                                                         htmlOutput("RNAseqGroup_venn"),
                                                         htmlOutput("ChIPseqGroup_venn"),
                                                         fluidRow(
                                                           column(4, downloadButton("download_RP_venn_table","Download summary table")),
                                                           column(4, downloadButton("download_selected_RP_venn_table","Download selected table"))
                                                         ),
                                                         fluidRow(
                                                           column(6, downloadButton("download_selected_int_peak_venn","Download selected peak file (bed)"))
                                                         ),
                                                         DTOutput('RP_table_venn')
                                         ),
                                         bsCollapsePanel(title="Functional enrichment analysis:",
                                                         value="function_panel",
                                                         fluidRow(
                                                           column(6, htmlOutput("intGroup_venn")),
                                                           column(6, htmlOutput("intGeneset_venn"))
                                                         ),
                                                         downloadButton("download_venn_int_enrichment", "Download"),
                                                         plotOutput('int_enrichment1_venn'),
                                                         fluidRow(
                                                           column(4, downloadButton("download_venn_int_enrichment_venn_table", "Download enrichment result"))
                                                         ),
                                                         dataTableOutput("int_enrichment_result_venn")
                                         )
                                         
                              )
                     ),
                     tabPanel("Combined heatmap",
                              textOutput("Spe_rnaseq2_venn"),
                              tags$head(tags$style("#Spe_rnaseq2_venn{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(6, htmlOutput("venn_heatmap_group"),
                                       tags$head(tags$style("#venn_heatmap_group{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),)
                              ),
                              fluidRow(
                                column(6, htmlOutput("rnaseq_DEGs_venn")),
                                column(6, htmlOutput("rnaseq_count_venn"))
                              ),
                              bsCollapse(id="z-score_count_venn",open="z-score_multiple_count_venn_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Uploaded DEG result files",
                                                         value="Uploaded_DEGs_venn",
                                                         dataTableOutput('rnaseq_DEGs_output_venn')
                                         ),
                                         bsCollapsePanel(title="z-score of uploaded count files",
                                                         value="z-score_multiple_count_venn_panel",
                                                         dataTableOutput('rnaseq_count_output_venn')
                                         )
                              ),
                              fluidRow(
                                column(4, htmlOutput("integrated_bw2_venn")),
                                column(4, htmlOutput("integrated_bw3_venn")),
                                column(4, htmlOutput("integrated_bw4_venn"))
                              ),
                              fluidRow(
                                column(4, actionButton("integrated_heatmapButton_venn", "Start"),
                                       tags$head(tags$style("#integrated_heatmapButton_venn{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                ),
                                column(4, downloadButton("download_integrated_heatmap_venn","Download heatmap"))
                              ),
                              div(
                                plotOutput('integrated_heatmap_venn', height = "100%"),
                                style = "height: calc(75vh  - 75px)"
                              ),
                     )
                   )
                 ) #sidebarLayout
               ) #tabPanel
      ),
      # Clustering -------------------------------------
      tabPanel("Clustering",
               sidebarLayout(
                 # Clustering---------------------------------
                 sidebarPanel(
                   radioButtons('data_file_type_clustering','Input:',
                                c('BigWig files'="Row1"
                                ),selected = "Row1"),
                   conditionalPanel(condition="input.data_file_type_clustering=='Row1'",
                                    fileInput("file1_clustering",
                                              strong(
                                                span("Select BigWig files"),
                                                span(icon("info-circle"), id = "icon1_clustering", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("bw","BigWig"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon1_clustering", "BigWig files (bw, BigWig):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                                            strong("Do not use it for anything else"),".<br>",
                                                            strong("Short names are recommended for file names.")), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type_clustering=='Row2'",
                                    fileInput("file_bam_clustering",
                                              strong(
                                                span("Select Bam files"),
                                                span(icon("info-circle"), id = "icon_bam_clustering", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("bam"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon_bam_clustering", "Bam files (bam):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                                            strong("Do not use it for anything else"),".<br>",
                                                            strong("Short names are recommended for file names.")), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type_clustering=='Row2'",
                                    radioButtons('Pair_or_single_clustering','Sequence type:',
                                                 c('Paired-end'="Paired-end",
                                                   'Single-end'="Single-end"),selected = "Paired-end"),         
                   ),
                   radioButtons('Genomic_region_clustering','Genomic region:',
                                c('Genome-wide'="Genome-wide",
                                  'Promoter'="Promoter"),selected = "Genome-wide"),
                   conditionalPanel(condition="input.Genomic_region_clustering=='Promoter'",
                                    fluidRow(
                                      column(5, numericInput("upstream_clustering", "upstream", value = 500, min = 0)),
                                      column(5, numericInput("downstream_clustering", "downstream", value = 500, min = 0))
                                    ),
                                    fileInput("genelist_file1_clustering",
                                              strong(
                                                span("Option: Select a gene list file for gene extraction"),
                                                span(icon("info-circle"), id = "icon_genelist_clustering", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt","csv","xlsx"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon_genelist_clustering", "gene list file (txt, csv, xlsx):", 
                                              content=paste("The first column is", strong("gene name"), ".<br>", 
                                                            "The second and subsequent columns do not affect the analysis.<br>", 
                                                            "File names do not use `.` other than the extension.<br><br>", 
                                                            img(src="genelist_input.png", width = 100,height = 300)), 
                                              placement = "right",options = list(container = "body"))
                   ),
                   conditionalPanel(condition="input.Genomic_region_clustering=='Genome-wide'",
                                    fileInput("peak_call_file1_clustering",
                                              strong(
                                                span("Select bed files"),
                                                span(icon("info-circle"), id = "icon2_clustering", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("bed","narrowPeak"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon2_clustering", "Bed files (bed, narrowPeak):", 
                                              content=paste("The first column is chromosome (e.g. chr3, chrY) name.<br>",
                                                            "The second column is start position on the chromosome.<br>",
                                                            "The third column is end position on the chromosome.<br>"), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(
                     column(12, selectInput("Species_clustering", "Species", species_list, selected = "not selected"))),
                   fluidRow(
                     column(4, numericInput("fc_clustering", "Fold Change", min   = 1, max   = NA, value = 2)),
                     column(4, numericInput("basemean_clustering", "Basemean", min   = 0, max   = NA, value = 0))
                   ),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "clustering_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("clustering_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("clustering_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("clustering_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>", 
                                           strong("Clustering:"), "height = 3.5, width = 9<br>", 
                                           strong("Correlation plot:"), "height = 5, width = 7<br>", 
                                           strong("Heatmap:"), "height = 10, width = 7<br>", 
                                           strong("Trackplot:"), "height = 4, width = 7 <br>",
                                           strong("Peak distribution:"), "height = 6, width = 10 <br>",
                                           strong("Peak pattern heatmap:"), "height = 6, width = 6 <br>",
                                           strong("Peak pattern line plot:"), "height = 5, width = 5 <br>",
                                           strong("Dotplot (GREAT):"), "height = 5, width = 6 <br>",
                                           strong("Region gene associations plot (GREAT):"), "height = 5, width = 12 <br>",
                                           strong("Motif dot plot:"), "height = 6, width = 7 <br>",
                                           strong("Boxplot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("Barplot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("KS-plot (with RNAseq):"), "height = 5, width = 7 <br>",
                                           strong("Dotplot (with RNAseq):"), "height = 5, width = 8 <br>"),trigger = "click"), 
                   actionButton("goButton_clustering", "example data (hg19)"),
                   tags$head(tags$style("#goButton_clustering{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             )
                   ) #sidebarPanel
                 ),
                 
                 # Main Panel -------------------------------------
                 mainPanel(
                   tabsetPanel(
                     type = "tabs",
                     tabPanel("Input Data",
                              bsCollapse(id="input_collapse_panel_clustering",open="bw_files_panel",multiple = TRUE,
                                         bsCollapsePanel(title="uploaded bigwig files:",
                                                         value="bw_files_panel",
                                                         DTOutput('input_bw_files_clustering') 
                                         ),
                                         bsCollapsePanel(title="uploaded peak call files:",
                                                         value="peak_call_files_panel",
                                                         dataTableOutput('input_peak_call_files_clustering')
                                         ),
                                         bsCollapsePanel(title="Raw count data:",
                                                         value="raw_count_panel",
                                                         column(4, textOutput("Spe_clustering"),
                                                                tags$head(tags$style("#Spe_clustering{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                                         downloadButton("download_raw_count_clustering_table", "Download raw count table"),
                                                         dataTableOutput('raw_count_table_clustering')
                                         )
                              )
                     ),
                     tabPanel("Result overview",
                              fluidRow(
                                column(4, downloadButton("download_clustering_PCA", "Download clustering analysis"))
                              ),
                              textOutput("clustering_pca_error"),
                              tags$head(tags$style("#clustering_pca_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("PCA_clustering"),
                              downloadButton("download_clustering_corrplot", "Download correlation plot"),
                              div(
                                plotOutput("correlationplot", height = "100%"),
                                style = "height: calc(75vh  - 75px)"
                              ), 
                              bsCollapse(id="input_collapse_clustering_DEG",open="PCA_panel",multiple = TRUE,
                                         bsCollapsePanel(title="PCA:",
                                                         value="PCA_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_clustering_PCA_table", "Download PCA table"))
                                                         ),
                                                         dataTableOutput("clustering_PCA_data")
                                         )
                              )),
                     tabPanel("k-means clustering",
                              fluidRow(
                                column(4, 
                                       htmlOutput("selectFC"),
                                       textOutput("filtered_region"),
                                       tags$head(tags$style("#filtered_region{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }")),
                                       htmlOutput("clustering_kmeans_num"),
                                       htmlOutput("kmeans_cv"),
                                       downloadButton("download_clustering_kmeans_heatmap", "Download heatmap"),
                                       actionButton("kmeans_start", "Start"),
                                       tags$head(tags$style("#kmeans_start{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))),
                                column(8, plotOutput("clustering_kmeans_heatmap"))
                              ),
                              htmlOutput("clustering_select_kmean"),
                              tags$head(tags$style("#clustering_select_kmean{color: black;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, downloadButton("download_peak_pattern_kmeans_heatmap", "Download peak heatmap"),
                                       htmlOutput("peak_pattern_kmeans_bw"),
                                       htmlOutput("peak_pattern_kmeans_heat_range")),
                                column(5, downloadButton("download_peak_pattern_kmeans_line", "Download peak aligned distribution"),
                                       htmlOutput("peak_pattern_kmeans_additional")),
                              ),
                              fluidRow(
                                column(6,  plotOutput("peak_pattern_kmeans_heatmap")),
                                column(6,  plotOutput("peak_pattern_kmeans_line"))
                              ),
                              bsCollapse(id="clustering_kmeans_collapse_panel",open="clustering_kmeans_extract_count",multiple = TRUE,
                                         bsCollapsePanel(title="selected cluster data:",
                                                         value="clustering_kmeans_extract_count",
                                                         downloadButton("download_clustering_kmeans_extract_count", "Download selected cluster"),
                                                         downloadButton("download_clustering_kmeans_extract_count_bed","Download selected cluster(bed)"),
                                                         DTOutput("clustering_kmeans_extract_table")
                                         ),
                                         bsCollapsePanel(title= "Trackplot:",
                                                         value="kmeans_track_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_clustering_trackplot", "Download trackplot")),
                                                           column(8, htmlOutput("trackplot_additional_clustering"))
                                                         ),
                                                         fluidRow(
                                                           column(8, htmlOutput("igv_uprange_clustering")),
                                                           column(4, htmlOutput("igv_ylim_clustering"))),
                                                         column(4, textOutput("Spe_clustering_track"),
                                                                tags$head(tags$style("#Spe_clustering_track{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                                         div(
                                                           plotOutput("trackplot_goi_clustering", height = "100%"),
                                                           style = "height: calc(75vh  - 75px)"
                                                         )      
                                         )
                              )
                     ),
                     tabPanel("with RNA-seq",
                              htmlOutput("warning_kmeans_RNA"),
                              tags$head(tags$style("#warning_kmeans_RNA{color: red;
                                 font-size: 25px;
            font-style: bold;
            }")),
                              bsCollapse(id="RNAseqresult_clustering_panel",open="RNAseqresult_clustering_table",
                                         bsCollapsePanel(title= p(span("RNA-seq result"),span(icon("info-circle"), id = "RNAseqresult_clustering", 
                                                                                              options = list(template = popoverTempate))),
                                                         bsPopover("RNAseqresult_clustering", "RNA-seq result file:", 
                                                                   content=paste("You can use a pair-wise RNA-seq DEG result file as input.<br>",
                                                                                 "First column must be gene name (Gene symbol or ENSEMBL ID).<br>", 
                                                                                 "The file must contain", strong("log2FoldChange"), "and", strong("padj"), "columns.<br>",
                                                                                 "<br>", 
                                                                                 img(src="input_format_volcano.png", width = 500,height = 230)), 
                                                                   placement = "right",options = list(container = "body")),
                                                         value="RNAseqresult_clustering_table",
                              fluidRow(
                                column(8, htmlOutput("clusteringRNAseqresult"))
                              ),
                              dataTableOutput('clustering_DEG_result')
                                         )),
                              htmlOutput("clustering_select_kmean_RNA"),
                              tags$head(tags$style("#clustering_select_kmean_RNA{color: black;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, htmlOutput("peak_distance_clustering")),
                                column(4, downloadButton("download_clusteringintbox","Download box plot"))
                              ),
                              textOutput("Spe_int_clustering"),
                              tags$head(tags$style("#Spe_int_clustering{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              textOutput("RNAseq_boxplot_clustering_error"),
                              tags$head(tags$style("#RNAseq_boxplot_clustering_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("int_box_clustering"),
                              fluidRow(
                                column(4, htmlOutput("DEG_fc_clustering")),
                                column(4, htmlOutput("DEG_fdr_clustering"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_clusteringintbar","Download bar plot")),
                                column(4, downloadButton("download_clusteringKSplot","Download ks plot"))
                              ),
                              fluidRow(
                                column(6,plotOutput("int_bar_clustering")),
                                column(6,plotOutput('ks_plot_clustering'))
                              ),
                              bsCollapse(id="input_collapse_clustering_RP",open="RP_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="int_Trackplot_panel",
                                                         fluidRow(
                                                           column(6, downloadButton("download_clustering_int_trackplot", "Download trackplot")),
                                                           column(6, htmlOutput("int_trackplot_additional_clustering"))
                                                         ),
                                                         fluidRow(
                                                           column(8, htmlOutput("int_igv_uprange_clustering")),
                                                           column(4, htmlOutput("int_igv_ylim_clustering"))),
                                                         plotOutput("int_trackplot_goi_clustering")
                                         ),
                                         bsCollapsePanel(title="Result table:",
                                                         value="RP_panel",
                                                         htmlOutput("RNAseqGroup_clustering"),
                                                         htmlOutput("ChIPseqGroup_clustering"),
                                                         fluidRow(
                                                           column(4, downloadButton("download_RP_clustering_table","Download summary table")),
                                                           column(4, downloadButton("download_selected_RP_clustering_table","Download selected table"))
                                                         ),
                                                         fluidRow(
                                                           column(6, downloadButton("download_selected_int_peak_clustering","Download selected peak file (bed)"))
                                                         ),
                                                         DTOutput('RP_table_clustering')
                                         ),
                                         bsCollapsePanel(title="Functional enrichment analysis:",
                                                         value="function_panel",
                                                         fluidRow(
                                                           column(6, htmlOutput("intGroup_clustering")),
                                                           column(6, htmlOutput("intGeneset_clustering"))
                                                         ),
                                                         downloadButton("download_clustering_int_enrichment", "Download"),
                                                         plotOutput('int_enrichment1_clustering'),
                                                         fluidRow(
                                                           column(4, downloadButton("download_clustering_int_enrichment_clustering_table", "Download enrichment result"))
                                                         ),
                                                         dataTableOutput("int_enrichment_result_clustering")
                                         )
                                         
                              )
                     )
                   )
                 ),

               ) #sidebarLayout
      ),
      # enrichment viewer -------------------------------------
      tabPanel("Enrichment viewer",
               sidebarLayout(
                 # enrichment viewer---------------------------------
                 sidebarPanel(
                   fileInput("enrich_data_file",
                             strong(
                               span("Select bed files (bed, narrowPeak)"),
                               span(icon("info-circle"), id = "icon_enrich1", 
                                    options = list(template = popoverTempate))
                             ),
                             accept = c("bed", "narrowPeak"),
                             multiple = FALSE,
                             width = "80%"),
                   bsPopover("icon_enrich1", "Uploaded file format (bed, narrowPeak): ", 
                             content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                           strong("Do not use it for anything else"),".<br>",
                                           strong("Short names are recommended for file names."),".<br>",
                                           "File names do not use `.` other than the extension.<br><br>"), 
                             placement = "right",options = list(container = "body")),
                   fluidRow(
                     column(12, selectInput("Species_enrich", "Species", species_list, selected = "not selected"))),
                   sliderInput("enrich_showCategory", "Most significant pathways",
                               min = 1, max = 20, value = 5,step = 1),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "enrich_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("enrich_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("enrich_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("enrich_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",
                                           "Dotplot:", "height = 5, width = 6.5 <br>", 
                                           "cnet plot:","height = 6, width = 6 <br><br>",
                                           "combined heatmap:", "height = 6, width = 10 <br>"), trigger = "click"), 
                   actionButton("goButton_enrich", "example data (hg19)"),
                   tags$head(tags$style("#goButton_enrich{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             )
                   ) #sidebarPanel
                 ),
                 
                 # Main Panel -------------------------------------
                 mainPanel(
                   tabsetPanel(
                     type = "tabs",
                     tabPanel("Input list",
                              dataTableOutput('enrichment_input')
                     ),
                     tabPanel("Peak distribution",
                              downloadButton("download_input_peak_distribution_enrich", "Download distribution"),
                              textOutput("Spe_dist_enrich"),
                              tags$head(tags$style("#Spe_dist_enrich{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("input_peak_distribution_enrich"),
                              bsCollapsePanel(title="Annotation:",
                                              value="annotation_enrich_panel",
                                              fluidRow(
                                                column(4, htmlOutput("annotation_select_enrich"))
                                              ),
                              DTOutput("enrich_annotation")
                              )
                     ),
                     tabPanel("GREAT",
                              fluidRow(
                                column(4, htmlOutput("Gene_set_enrich")),
                                column(4, downloadButton("download_enrich_enrichment", "Download"))
                              ),
                              plotOutput("enrichment_enrich"),
                              bsCollapse(id="input_collapse_enrich_enrich",open="ORA_enrich_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Enrichment result:",
                                                         value="ORA_enrich_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_enrich_enrichment_table", "Download enrichment result"))
                                                         ),
                                                         dataTableOutput('enrich_enrichment_result')
                                         )
                              ),
                              fluidRow(
                                column(6, htmlOutput("whichGenes_enrich"), 
                                       htmlOutput("whichGeneSet_enrich")),
                                column(4, downloadButton("download_enrich_region_gene_associations_plot", "Download gene associations plot"))
                              ),
                              plotOutput('region_gene_associations_enrich_plot'),
                              fluidRow(
                                column(4, downloadButton("download_enrich_region_gene_associations", "Download table data"))
                              ),
                              DTOutput('region_gene_enrich_associations')
                     ),
                     tabPanel("HOMER",
                              fluidRow(
                                column(3, downloadButton("download_motif_enrich_plot", "Download motif plot")),
                                column(3, downloadButton("download_homer_report_enrich", "Download homer report"),
                                       tags$head(tags$style("#download_homer_report_enrich{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),)
                              ),
                              textOutput("Spe_motif_enrich"),
                              tags$head(tags$style("#Spe_motif_enrich{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, htmlOutput("homer_unknown_enrich")),
                                column(4, htmlOutput("homer_size_enrich")),
                                column(4, htmlOutput("homer_size2_enrich"))
                              ),
                              fluidRow(
                                column(4, htmlOutput("homer_bg_enrich")),
                                column(4, htmlOutput("homer_bg2_enrich"))
                              ),
                              fluidRow(
                                column(4, actionButton("motifButton_enrich", "Start"),
                                       tags$head(tags$style("#motifButton_enrich{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                )
                              ),
                              htmlOutput("homer_fdr_enrich"),
                              plotOutput("motif_enrich_plot"),
                              bsCollapse(id="Promoter_motif_enrich_collapse_panel",open="motif_enrich_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Known motif",
                                                         value="motif_enrich_result_table",
                                                         downloadButton("download_motif_enrich_table", "Download known motif result"),
                                                         DTOutput('motif_enrich_result')
                                         ),
                                         bsCollapsePanel(title= "de_novo_motif",
                                                         value="Promoter_motif_region_panel",
                                                         downloadButton("download_denovo_motif_enrich_table", "Download de novo motif"),
                                                         dataTableOutput("denovo_motif_enrich_result")
                                         )
                              )
                     ),
                     tabPanel("Combined heatmap",
                              textOutput("Spe_rnaseq2_enrich"),
                              tags$head(tags$style("#Spe_rnaseq2_enrich{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(6, htmlOutput("integrated_bw1_enrich"))
                              ),
                              fluidRow(
                                column(6, htmlOutput("rnaseq_DEGs_enrich")),
                                column(6, htmlOutput("rnaseq_count_enrich"))
                              ),
                              bsCollapse(id="z-score_count_enrich",open="z-score_multiple_count_enrich_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Uploaded DEG result files",
                                                         value="Uploaded_DEGs_enrich",
                                                         dataTableOutput('rnaseq_DEGs_output_enrich')
                                         ),
                                         bsCollapsePanel(title="z-score of uploaded count files",
                                                         value="z-score_multiple_count_enrich_panel",
                                                         dataTableOutput('rnaseq_count_output_enrich')
                                         )
                              ),
                              fluidRow(
                                column(4, htmlOutput("integrated_bw2_enrich")),
                                column(4, htmlOutput("integrated_bw3_enrich")),
                                column(4, htmlOutput("integrated_bw4_enrich"))
                              ),
                              fluidRow(
                                column(4, actionButton("integrated_heatmapButton_enrich", "Start"),
                                       tags$head(tags$style("#integrated_heatmapButton_enrich{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                ),
                                column(4, downloadButton("download_integrated_heatmap_enrich","Download heatmap"))
                              ),
                              div(
                                plotOutput('integrated_heatmap_enrich', height = "100%"),
                                style = "height: calc(75vh  - 75px)"
                              ),
                     )
                   )
                 )
               ) #sidebarLayout
      ),
      #Instruction--------------------------
      navbarMenu("More",
                 tabPanel("Bedtools",
                          sidebarLayout(
                            # Bedtools---------------------------------
                            sidebarPanel(
                              radioButtons('data_file_type_bed','Function:',
                                           choiceNames = list(
                                             tagList(
                                               tags$span("Merge intervals"),
                                               tags$span(icon("info-circle"), id = "icon_bedtool1", style = "color: black;")
                                             ), 
                                             tagList(
                                               tags$span("Intersect"),
                                               tags$span(icon("info-circle"), id = "icon_bedtool2", style = "color: black;")
                                             ),
                                             tagList(
                                               tags$span("Subtract"),
                                               tags$span(icon("info-circle"), id = "icon_bedtool3", style = "color: black;")
                                             )
                                           ),
                                           choiceValues = c("type1","type3","type2"),
                                           selected="type1"),
                              bsPopover("icon_bedtool1", "Merge interval:", 
                                        content=paste("This function can fill the gap sequence.<br>", 
                                                      img(src="merge_interval.png", width = 400,height = 180))), 
                              bsPopover("icon_bedtool3", "Subtract:", 
                                        content=paste("This function searches for features in B that overlap A.<br>", 
                                                      "If an overlapping feature is found in B, the overlapping portion is removed from A and the remaining portion of A is reported.",
                                                      img(src="subtract.png", width = 400,height = 180))), 
                              bsPopover("icon_bedtool2", "Intersect:", 
                                        content=paste("By ",strong("default"),", if an overlap is found, this function reports the shared interval between the two overlapping features.<br>", 
                                                      "The ",strong("wa"), " option reports the original entry in A for each overlap.<br>",
                                                      "The ",strong("v"), " option only reports those entries in A that have no overlap in B.<br>",
                                                      img(src="intersect.png", width = 400,height = 180))), 
                              fileInput("bed_data_file",
                                        strong(
                                          span("A: Select bed files (bed, narrowPeak)",style="color: red")
                                        ),
                                        accept = c("bed", "narrowPeak"),
                                        multiple = FALSE,
                                        width = "80%"),
                              conditionalPanel(condition="input.data_file_type_bed=='type1'",
                                               htmlOutput("bed_merge_dist")
                              ),
                              conditionalPanel(condition="input.data_file_type_bed=='type2' || input.data_file_type_bed=='type3'",
                                               fileInput("bed_data_file2",
                                                         strong(
                                                           span("B :Select bed files  (bed, narrowPeak)",style="color: blue")
                                                         ),
                                                         accept = c("bed", "narrowPeak"),
                                                         multiple = FALSE,
                                                         width = "80%"),
                              ),
                              conditionalPanel(condition="input.data_file_type_bed=='type3'",
                                               radioButtons('intersect_bed','Mode:',
                                                            c('default'="default",
                                                              'wa'="wa",
                                                              'v'="exclude"
                                                            ),selected = "default"),           
                              ),
                              actionButton("goButton_bed", "example data (hg19)"),
                              tags$head(tags$style("#goButton_bed{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                                        tags$style("
          body {
            padding: 0 !important;
          }"
                                        )
                              ) 
                            ),
                            
                            # Main Panel -------------------------------------
                            mainPanel(
                              textOutput("bed_file1_warning"),
                              tags$head(tags$style("#bed_file1_warning{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              textOutput("bed_file2_warning"),
                              tags$head(tags$style("#bed_file2_warning{color: blue;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              dataTableOutput('bed_input'),
                              downloadButton("download_bed", "Download processed bed files"),
                              tags$head(tags$style("#download_bed{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                            )
                          ) #sidebarLayout
                 ),
                 tabPanel("Reference",
                          fluidRow(
                            column(12,
                                   h2("Reference:"),
                                   "- Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert and Barbara Borges (2021). shiny: Web Application Framework for R. R package version 1.7.1. https://CRAN.R-project.org/package=shiny",br(),
                                   "- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)",br(),
                                   "- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro",br(),
                                   "- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141",br(),
                                   "- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609",br(),
                                   "- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.", br(),
                                   "- Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. 'Benchmark and integration of resources for the estimation of human transcription factor activities.' Genome Research. 2019. DOI: 10.1101/gr.240663.118.", br(),
                                   "- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi",br(),
                                   "- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.",br(),
                                   "- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.",br(),
                                   "- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.",br(),
                                   "- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.", br(),
                                   "- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr",br(),
                                   "- Mitjavila-Ventura A (2022). _plotmics: Visualization of Omics And Sequencing Data In R_. https://amitjavilaventura.github.io/plotmics/index.html,
  https://github.com/amitjavilaventura/plotmics.",br(),
                                   "- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr",br(),
                                   "- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr",br(),
                                   "- Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118",br(),
                                   "- Team BC, Maintainer BP (2019). _TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package for TxDb object(s)_. R package version 3.10.0.",br(),
                                   "- Carlson M, Maintainer BP (2015). _TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package for TxDb object(s)_. R package version 3.2.2.",br(),
                                   "- Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing
  and microarray studies. Nucleic Acids Research 43(7), e47.",br(),
                                   "- Liao Y, Smyth GK and Shi W (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.
  Nucleic Acids Research 47(8), e47.",br(),
                                   "- Morgan M, Pagès H, Obenchain V, Hayden N (2022). _Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import_. R package version
  2.12.0, <https://bioconductor.org/packages/Rsamtools>.",br(),
                                   "- Gu Z (2022). _rGREAT: GREAT Analysis - Functional Enrichment on Genomic Regions_. https://github.com/jokergoo/rGREAT,
  http://great.stanford.edu/public/html/.",br(),
                                   "- - Shang, G.-D., Xu, Z.-G., Wan, M.-C., Wang, F.-X. & Wang, J.-W.  FindIT2: an R/Bioconductor package to identify influential transcription factor and
  targets based on multi-omics data.  BMC Genomics 23, 272 (2022)",br(),
                                   "- Wagih O (2017). _ggseqlogo: A 'ggplot2' Extension for Drawing Publication-Ready Sequence Logos_. R package version 0.1,
  <https://CRAN.R-project.org/package=ggseqlogo>.",br(),
                                   "- Amezquita R (2022). _marge: API for HOMER in R for Genomic Analysis using Tidy Conventions_. R package version 0.0.4.9999.",br(),
                                   "- Zeileis A, Hornik K, Murrell P (2009). Escaping RGBland: Selecting Colors for Statistical Graphics. _Computational Statistics & Data Analysis_,
  *53*(9), 3259-3270. doi:10.1016/j.csda.2008.11.033 <https://doi.org/10.1016/j.csda.2008.11.033>.",br(),
                                   "- Kassambara A (2022). _ggcorrplot: Visualization of a Correlation Matrix using 'ggplot2'_. R package version 0.1.4,
  <https://CRAN.R-project.org/package=ggcorrplot>.",br(),
                                   "- Zheng H (2021). _bedtorch: R package for fast BED-file manipulation_. R package version 0.1.12.12, <https://github.com/haizi-zh/bedtorch>.",br(),
                                   "- Hahne F, Ivanek R. Visualizing Genomic Data Using Gviz and Bioconductor. Methods Mol Biol. 1418:335-51 (2016).",br()
                            )
                          )
                 )
      )
    )
  ))
