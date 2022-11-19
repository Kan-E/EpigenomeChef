popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'

shinyUI(
  fluidPage(
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
                      "or", a("contact us", href = "kaneto@kumamoto-u.ac.jp"),".",align="center",width=4)
      ),
      "",
      id='navBar',
      tabPanel("EpigenomeChef" ,value='Title', icon = icon("utensils"),
               fluidRow(
                 column(12,
                        br(),br(),
                        h1(strong("EpigenomeChef"),align="center"),br(),
                        p("EpigenomeChef is a web-based application for automated, systematic, and integrated epigenetic analysis.",
                          align="center"),br(),br(),style={'background-color:mintcream;font-size: 16px;'},
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
                                  'Option1: '="Row2",
                                  'Option2: '="Row11"
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
                                                            strong("Do not use it for anything else"),".<br><br>"), 
                                              placement = "right",options = list(container = "body")),
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
                                              content=paste("File names must be the same as bigwig files.<br><br>"), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(
                     column(12, selectInput("Species", "Species", species_list, selected = "not selected"))),
                   selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH"),
                   fluidRow(
                     column(4, numericInput("fc", "Fold Change", min   = 1, max   = NA, value = 5)),
                     column(4, numericInput("fdr", "FDR", min   = 0, max   = 1, value = 0.01, step = 0.005)),
                     column(4, numericInput("basemean", "Basemean", min   = 0, max   = NA, value = 0))
                   ),
                   
                   fluidRow(
                     column(5, numericInput("pair_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("pair_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("pair_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Clustering:"), "height = 3.5, width = 9<br>", strong("MA-plot:"), "height = 4, width = 7 <br>",strong("Volcano plot:"), "height = 5, width = 5 <br>",pdfSize_for_GOI,
                                           strong("Enrichment analysis:"), "height = 10, width = 12 <br>"),trigger = "click"), 
                   actionButton("goButton", "example data (hg19)"),
                   tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             ))
                 ), #sidebarPanel
                 
                 # Main Panel -------------------------------------
                 mainPanel(
                   tabsetPanel(
                     type = "tabs",
                     tabPanel("Input Data",
                              bsCollapse(id="input_collapse_panel",open="bw_files_panel",multiple = TRUE,
                                         bsCollapsePanel(title="uploaded bigwig files:",
                                                         value="bw_files_panel",
                                                         dataTableOutput('input_bw_files') 
                                         ),
                                         bsCollapsePanel(title="uploaded peak call files:",
                                                         value="peak_call_files_panel",
                                                         dataTableOutput('input_peak_call_files')
                                         ),
                                         bsCollapsePanel(title="Raw count data:",
                                                         value="raw_count_panel",
                                                         column(4, textOutput("Spe"),
                                                                tags$head(tags$style("#Spe{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                                         downloadButton("download_raw_count_table", "Download raw count table"),
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
                                column(4, downloadButton("download_pair_volcano", "Download volcano plot and heatmap"))
                              ),
                              fluidRow(
                                column(6, htmlOutput("volcano_x")),
                                column(6, htmlOutput("volcano_y"))),
                              fluidRow(
                                column(8, plotOutput("volcano1")),
                                column(4, plotOutput("GOIheatmap"))
                              ),
                              bsCollapse(id="input_collapse_pair_DEG",open="DEG_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="Trackplot_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_trackplot", "Download trackplot")),
                                                           column(4, htmlOutput("trackplot_additional"))
                                                         ),
                                                         fluidRow(
                                                           column(4, htmlOutput("igv_uprange")),
                                                           column(4, htmlOutput("igv_downrange")),
                                                           column(4, htmlOutput("igv_ylim"))),
                                                         plotOutput("trackplot_goi")
                                         ),
                                         bsCollapsePanel(title="DESeq2 result:",
                                                         value="DEG_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_DEG_result", "Download DEG result"))
                                                         ),
                                                         DTOutput("DEG_result")
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
                              downloadButton("download_input_peak_distribution", "Download input peak distribution"),
                              plotOutput("input_peak_distribution"),
                              downloadButton("download_deg_peak_distribution", "Download DAR peak distribution"),
                              plotOutput("deg_peak_distribution"),
                     ),
                     tabPanel("Motif analysis",
                              fluidRow(
                                column(3, downloadButton("download_motif_plot", "Download motif plot"))
                              ),
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
                              plotOutput("motif_plot"),
                              bsCollapse(id="Promoter_motif_collapse_panel",open="motif_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Motif table:",
                                                         value="motif_result_table",
                                                         downloadButton("download_motif_table", "Download motif enrichment result"),
                                                         DTOutput('motif_result')
                                         ),
                                         bsCollapsePanel(title= p(span("Motif region"),span(icon("info-circle"), id = "icon_promoter_motif_region", 
                                                                                            options = list(template = popoverTempate))),
                                                         bsPopover("icon_promoter_motif_region", "Motif region:", 
                                                                   content=paste("Please select genes in", strong("k-means clustering result"),".<br><br>",
                                                                                 img(src="enrich_motif.png", width = 450,height = 600)), 
                                                                   placement = "right",options = list(container = "body")),
                                                         value="Promoter_motif_region_panel",
                                                         downloadButton("download_promoter_motif_region", "Download motif region"),
                                                         dataTableOutput("promoter_motif_region_table")
                                         )
                              )
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
                     tabPanel("Functional analysis",
                              fluidRow(
                                column(4, textOutput("Spe1"),
                                       tags$head(tags$style("#Spe1{color: red;
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
                                column(4, downloadButton("download_region_gene_associations_plot", "Download table data"))
                              ),
                              plotOutput('region_gene_associations_plot'),
                              fluidRow(
                                column(4, downloadButton("download_region_gene_associations", "Download table data"))
                              ),
                              DTOutput('region_gene_associations')
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
                               span("Select peak call files"),
                               span(icon("info-circle"), id = "icon_venn1", 
                                    options = list(template = popoverTempate))
                             ),
                             accept = c("bed","narrowPeak"),
                             multiple = TRUE,
                             width = "80%"),
                   bsPopover("icon_venn1", "peak call files (bed, narrowPeak):", 
                             content=paste(""), 
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
                             content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                           strong("Do not use it for anything else"),".<br><br>"), 
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
                                           "Default size: <br>",strong("Venn diagram:"), "height = 3.5, width = 9<br>", strong("Integrated heatmap:"), "height = 8, width = 8 <br>",
                                           strong("Enrichment analysis:"), "height = 6, width = 8 <br>",strong("cnet plot:"), "height = 6, width = 6 <br>"),trigger = "click"), 
                   actionButton("goButton_venn", "example data"),
                   tags$head(tags$style("#goButton{color: black;
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
                     tabPanel("Venn diagram",
                              fluidRow(
                                column(4, downloadButton("download_vennplot", "Download venn diagram"))
                              ),
                              plotOutput("venn"),
                              downloadButton("download_venn_result", "Download venn result"),
                              dataTableOutput("venn_result")
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
                              plotOutput("venn_peak_distribution"),
                              bsCollapse(id="int_result_collapse_panel",open="selected_intersect_annotation_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Trackplot:",
                                                         value="Trackplot_venn_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_venn_trackplot", "Download trackplot"))
                                                         ),
                                                         fluidRow(
                                                           column(4, htmlOutput("igv_venn_uprange")),
                                                           column(4, htmlOutput("igv_venn_downrange")),
                                                           column(4, htmlOutput("igv_venn_ylim"))),
                                                         plotOutput("trackplot_venn_goi")
                                         ),
                                         bsCollapsePanel(title="selected intersect data:",
                                                         value="selected_intersect_annotation_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_selected_intersect_annotation_table", "Download selected intersection table (.txt)")),
                                                           column(4, downloadButton("download_selected_intersect_annotation_table_bed", "Download selected intersection table (bed format)"))
                                                         ),
                                                         
                                                         DTOutput("selected_intersect_annotation")
                                         )
                              )
                     ),
                     tabPanel("Motif analysis",
                              fluidRow(
                                column(8, htmlOutput("venn_whichGroup1")),
                                column(4, downloadButton("download_motif_venn_plot", "Download motif plot"))
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
                              plotOutput("motif_venn_plot"),
                              bsCollapse(id="Promoter_motif_venn_collapse_panel",open="motif_venn_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Motif table:",
                                                         value="motif_venn_result_table",
                                                         downloadButton("download_motif_venn_table", "Download motif enrichment result"),
                                                         DTOutput('motif_venn_result')
                                         ),
                                         bsCollapsePanel(title= p(span("Motif region"),span(icon("info-circle"), id = "icon_promoter_motif_venn_region", 
                                                                                            options = list(template = popoverTempate))),
                                                         bsPopover("icon_promoter_motif_venn_region", "Motif region:", 
                                                                   content=paste("Please select genes in", strong("k-means clustering result"),".<br><br>",
                                                                                 img(src="enrich_motif.png", width = 450,height = 600)), 
                                                                   placement = "right",options = list(container = "body")),
                                                         value="Promoter_motif_region_venn_panel",
                                                         downloadButton("download_promoter_motif_venn_region", "Download motif region"),
                                                         dataTableOutput("promoter_motif_region_venn_table")
                                         )
                              )
                     ),
                     tabPanel("Functional analysis",
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
                     )
                   )
                 ) #sidebarLayout
               ) #tabPanel
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
                             content=paste("File names do not use `.` other than the extension.<br><br>"), 
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
                                           "Default size: <br>","Dotplot:", "height = 5, width = 6.5 <br>", "cnet plot:","height = 6, width = 6 <br><br>"), trigger = "click"), 
                   actionButton("goButton_enrich", "example data (mouse)"),
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
                     tabPanel("Functional analysis",
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
                     tabPanel("Motif analysis",
                              fluidRow(
                                column(4, downloadButton("download_motif_enrich_plot", "Download motif plot"))
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
                              plotOutput("motif_enrich_plot"),
                              bsCollapse(id="Promoter_motif_enrich_collapse_panel",open="motif_enrich_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Motif table:",
                                                         value="motif_enrich_result_table",
                                                         downloadButton("download_motif_enrich_table", "Download motif enrichment result"),
                                                         DTOutput('motif_enrich_result')
                                         ),
                                         bsCollapsePanel(title= p(span("Motif region"),span(icon("info-circle"), id = "icon_promoter_motif_enrich_region", 
                                                                                            options = list(template = popoverTempate))),
                                                         bsPopover("icon_promoter_motif_enrich_region", "Motif region:", 
                                                                   content=paste("Please select genes in", strong("k-means clustering result"),".<br><br>",
                                                                                 img(src="enrich_motif.png", width = 450,height = 600)), 
                                                                   placement = "right",options = list(container = "body")),
                                                         value="Promoter_motif_region_enrich_panel",
                                                         downloadButton("download_promoter_motif_enrich_region", "Download motif region"),
                                                         dataTableOutput("promoter_motif_region_enrich_table")
                                         )
                              )
                     ),
                   )
                 )
               ) #sidebarLayout
      ),
      #Instruction--------------------------
      navbarMenu("More",
                 tabPanel("Reference",
                          fluidRow(
                            column(12,
                                   h2("Reference:"),
                                   "- Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert and Barbara Borges (2021). shiny: Web Application Framework for R. R package version 1.7.1. https://CRAN.R-project.org/package=shiny",br(),
                                   "- Ning Leng and Christina Kendziorski (2020). EBSeq: An R package for gene and isoform
  differential expression analysis of RNA-seq data. R package version 1.30.0.",br(),
                                   "- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)",br(),
                                   "- Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential
  expression analysis of digital gene expression data. Bioinformatics 26, 139-140",br(),
                                   "- Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg and Wolfgang Huber (2016): Data-driven hypothesis
  weighting increases detection power in genome-scale multiple testing. Nature Methods 13:577,
  doi: 10.1038/nmeth.3885",br(),
                                   "- John D. Storey, Andrew J. Bass, Alan Dabney and David Robinson (2021). qvalue: Q-value
  estimation for false discovery rate control. R package version 2.26.0.
  http://github.com/jdstorey/qvalue",br(),
                                   "- Pantano L (2022). DEGreport: Report of DEG analysis. R package version 1.32.0, http://lpantano.github.io/DEGreport", br(),
                                   "- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro",br(),
                                   "- Konopka T (2022). _umap: Uniform Manifold Approximation and Projection_. R package version 0.2.8.0, <https://CRAN.R-project.org/package=umap>.", br(),
                                   "- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141",br(),
                                   "- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609",br(),
                                   "- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R
  package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.", br(),
                                   "- Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. 'Benchmark and integration of resources for the estimation of human
  transcription factor activities.' Genome Research. 2019. DOI: 10.1101/gr.240663.118.", br(),
                                   "- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi",br(),
                                   "- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Rn.eg.db: Genome wide annotation for Rat. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Xl.eg.db: Genome wide annotation for Xenopus. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Dm.eg.db: Genome wide annotation for Fly. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Ce.eg.db: Genome wide annotation for Worm. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2022). org.Bt.eg.db: Genome wide annotation for Bovine. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Cf.eg.db: Genome wide annotation for Canine. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Dr.eg.db: Genome wide annotation for Zebrafish. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Gg.eg.db: Genome wide annotation for Chicken. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Mmu.eg.db: Genome wide annotation for Rhesus. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Pt.eg.db: Genome wide annotation for Chimp. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Sc.sgd.db: Genome wide annotation for Yeast. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Ss.eg.db: Genome wide annotation for Pig. R package version 3.15.0.",br(),
                                   "- Morgan M, Shepherd L (2022). AnnotationHub: Client to access AnnotationHub resources. R package version 3.4.0.",br(),
                                   "- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.",br(),
                                   "- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.",br(),
                                   "- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.", br(),
                                   "- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr",br(),
                                   "- Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn",br(),
                                   "- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr",br(),
                                   "- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr",br(),
                                   "- Machlab D, Burger L, Soneson C, Rijli FM, Schübeler D, Stadler MB. monaLisa: an R/Bioconductor package for identifying regulatory motifs. Bioinformatics (2022).",br(),
                                   "- Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118",br(),
                                   "- Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel:
  Bioconductor facilities for parallel evaluation_. R package version 1.30.3,
  <https://github.com/Bioconductor/BiocParallel>.",br(),
                                   "- Morgan M, Obenchain V, Hester J, Pagès H (2022). _SummarizedExperiment:
  SummarizedExperiment container_. R package version 1.26.1,
  <https://bioconductor.org/packages/SummarizedExperiment>.",br(),
                                   "- Baranasic D (2020). _JASPAR2020: Data package for JASPAR database (version 2020)_. R
  package version 0.99.10, <http://jaspar.genereg.net/>.",br(),
                                   "- Team BC, Maintainer BP (2019). _TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package
  for TxDb object(s)_. R package version 3.10.0.",br(),
                                   "- Team TBD (2021). _BSgenome.Mmusculus.UCSC.mm10: Full genome sequences for Mus musculus
  (UCSC version mm10, based on GRCm38.p6)_. R package version 1.4.3.",br(),
                                   "- Carlson M, Maintainer BP (2015). _TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package
  for TxDb object(s)_. R package version 3.2.2.",br(),
                                   "- Team TBD (2020). _BSgenome.Hsapiens.UCSC.hg19: Full genome sequences for Homo sapiens
  (UCSC version hg19, based on GRCh37.p13)_. R package version 1.4.3.",br(),
                                   "Tan, G., and Lenhard, B. (2016). TFBSTools: an R/bioconductor package for transcription factor
  binding site analysis. Bioinformatics 32, 1555-1556.",br(),
                            )
                          )
                 )
      )
    )
  ))
