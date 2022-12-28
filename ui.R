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
                 ),
                 column(12,br(),
                        h3("Pair-wise DAR detects and visualizes differentially accessible regions"),br(),
                        img(src="pair-wise_DAR.png", width = 1200,height = 1600),br(),br(),hr(),
                        h3("Venn diagram detects and visualizes the overlap between DARs from multiple datasets"),br(),
                        img(src="Venn.png", width = 1200,height = 700),br(),br(),hr(), 
                        h3("Clustering identifies similar samples and DNA binding patterns by clustering methods"),br(),
                        img(src="Clustering.png", width = 1000,height = 900),br(),hr(),
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
                                                            strong("Do not use it for anything else"),".<br><br>"), 
                                              placement = "right",options = list(container = "body")),
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
                                                            strong("Do not use it for anything else"),".<br><br>"), 
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
                                              content=paste("File names must be the same as bigwig files.<br><br>"), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(
                     column(12, selectInput("Species", "Species", species_list, selected = "not selected"))),
                   conditionalPanel(condition="input.data_file_type=='Row2'",
                   selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")
                   ),
                   fluidRow(
                     column(4, numericInput("fc", "Fold Change", min   = 1, max   = NA, value = 2)),
                     column(4, numericInput("fdr", "FDR", min   = 0, max   = 1, value = 0.1, step = 0.005)),
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
                                         bsCollapsePanel(title="uploaded bigwig or bam files:",
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
                                column(4, downloadButton("download_pair_volcano", "Download volcano plot")),
                                column(4, downloadButton("download_pair_heatmap", "Download heatmap"))
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
                                         bsCollapsePanel(title="DESeq2 result:",
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
                              fluidRow(
                                column(8, htmlOutput("pairRNAseqresult"))
                              ),
                              dataTableOutput('pair_DEG_result'),
                              htmlOutput("peak_distance"),
                              fluidRow(
                                column(4, downloadButton("download_pairintbox","Download box plot"))
                              ),
                              textOutput("Spe_int"),
                              tags$head(tags$style("#Spe_int{color: red;
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
                     tabPanel("with RNA-seq",
                              fluidRow(
                                column(8, htmlOutput("vennRNAseqresult"))
                              ),
                              dataTableOutput('venn_DEG_result'),
                              htmlOutput("venn_select_RNA"),
                              tags$head(tags$style("#venn_select_RNA{color: black;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, htmlOutput("peak_distance_venn")),
                                column(4, downloadButton("download_vennintbox","Download box plot"))
                              ),
                              plotOutput("int_box_venn"),
                              fluidRow(
                                column(4, htmlOutput("DEG_fc_venn")),
                                column(4, htmlOutput("DEG_fdr_venn"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_vennintbar","Download bar plot")),
                                column(4, downloadButton("download_vennKSplot","Download ks plot"))
                              ),
                              textOutput("Spe_int_venn"),
                              tags$head(tags$style("#Spe_int_venn{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
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
                                c('BigWig files'="Row1",
                                  'Bam files'="Row2"
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
                                                            strong("Do not use it for anything else"),".<br><br>"), 
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
                                    bsPopover("icon1_clustering", "Bam files (bam):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),"in file names.<br>", 
                                                            strong("Do not use it for anything else"),".<br><br>"), 
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
                                    )
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
                                    bsPopover("icon2_clustering", "peak call files (bed, narrowPeak):", 
                                              content=paste("File names must be the same as bigwig files.<br><br>"), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(
                     column(12, selectInput("Species_clustering", "Species", species_list, selected = "not selected"))),
                   fluidRow(
                     column(5, numericInput("clustering_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("clustering_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("clustering_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>","Dotplot:", "height = 5, width = 6.5 <br>", "cnet plot:","height = 6, width = 6 <br><br>"), trigger = "click"), 
                   actionButton("goButton_clustering", "example data (mouse)"),
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
                                         bsCollapsePanel(title="uploaded bigwig or bam files:",
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
                              plotOutput("PCA_clustering"),
                              fluidRow(
                                column(6, htmlOutput("clustering_umap_n"),
                                       downloadButton("download_clustering_umap", "Download umap"),
                                       textOutput("clustering_umap_error"),
                                       tags$head(tags$style("#clustering_umap_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")))
                              ),
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
                                column(4, htmlOutput("clustering_kmeans_num"),
                                       htmlOutput("kmeans_cv"),
                                       downloadButton("download_clustering_kmeans_heatmap", "Download heatmap")),
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
                                         bsCollapsePanel(title= p(span("Track plot"),span(icon("info-circle"), id = "icon_clustering_kmeans_boxplot", 
                                                                                       options = list(template = popoverTempate))),
                                                         bsPopover("icon_clustering_kmeans_boxplot", "Boxplot of GOI:", 
                                                                   content=paste("Please select genes in", strong("k-means clustering result"),".<br><br>",
                                                                                 img(src="clustering_GOIboxplot.png", width = 450,height = 640)), 
                                                                   placement = "right",options = list(container = "body")),
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
                              fluidRow(
                                column(8, htmlOutput("clusteringRNAseqresult"))
                              ),
                              dataTableOutput('clustering_DEG_result'),
                              htmlOutput("clustering_select_kmean_RNA"),
                              tags$head(tags$style("#clustering_select_kmean_RNA{color: black;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              fluidRow(
                                column(4, htmlOutput("peak_distance_clustering")),
                                column(4, downloadButton("download_clusteringintbox","Download box plot"))
                              ),
                              plotOutput("int_box_clustering"),
                              fluidRow(
                                column(4, htmlOutput("DEG_fc_clustering")),
                                column(4, htmlOutput("DEG_fdr_clustering"))
                              ),
                              fluidRow(
                                column(4, downloadButton("download_clusteringintbar","Download bar plot")),
                                column(4, downloadButton("download_clusteringKSplot","Download ks plot"))
                              ),
                              textOutput("Spe_int_clustering"),
                              tags$head(tags$style("#Spe_int_clustering{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
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
                              actionButton("goButton_bed", "example data (mouse)"),
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
