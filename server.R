popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=10000000*1024^2)
  observeEvent(input$close, {
    session$close()
    js$closeWindow()
    stopApp()
    cat(sprintf("Closing session %s\n", session$token))
    lapply(paste("package:",names(sessionInfo()$otherPkgs),sep=""),
           detach,character.only = TRUE,
           unload = TRUE)
  })
  
  
  output$Spe <- renderText({
    if(input$Species == "not selected" && input$Genomic_region == "Promoter") print("Please select 'Species'")
  })
  output$Spe_track <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_dist <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  observeEvent(input$goButton,({
    updateSelectInput(session,inputId = "Species","Species",species_list, selected = "Homo sapiens (hg19)")
  }))
  output$Spe_dist_promoter <- renderText({
    if(input$Genomic_region == "Promoter") print("In the case of 'Promoter' mode, this function is not available. This function is for 'Genome-wide' mode. ")
  })
  output$Spe_great_promoter <- renderText({
    if(input$Genomic_region == "Promoter") print("In the case of 'Promoter' mode, this function is not available. This function is for 'Genome-wide' mode. ")
  })
  output$Spe_motif <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_GREAT <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_int <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_rnaseq2 <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  
  
  
  
  observeEvent(pre_bw_count(),({
    updateSelectizeInput(session,inputId = "sample_order","Sample order:",
                         choices = colnames(pre_bw_count()),selected = colnames(pre_bw_count()))
  }))
  output$sample_order_pair_track <- renderUI({
    if(!is.null(pre_track_additional_files())){
      selectInput(inputId = "sample_order_pair_track","Sample order:",
                  choices =  names(pre_track_additional_files()),
                  selected = names(pre_track_additional_files()),multiple = TRUE)
    }
  })
  output$sample_order_pair_pattern <- renderUI({
    if(!is.null(pre_up_additional())){
      selectInput(inputId = "sample_order_pair_pattern","Sample order:",
                  choices =  names(pre_up_additional()),
                  selected = names(pre_up_additional()),multiple = TRUE)
    }
  })
  output$sample_order_pair_track_int <- renderUI({
    if(!is.null(pre_int_track_additional_files())){
      selectInput(inputId = "sample_order_pair_track_int","Sample order:",
                  choices =  names(pre_int_track_additional_files()),
                  selected = names(pre_int_track_additional_files()),multiple = TRUE)
    }
  })
  output$with_sample_order_pair_comb1 <- renderUI({
    if(!is.null(with_pre_integrated_additional1())){
      selectInput(inputId = "with_sample_order_pair_comb1","Sample order (blue):",
                  choices =  names(with_pre_integrated_additional1()),
                  selected = names(with_pre_integrated_additional1()),multiple = TRUE)
    }
  })
  output$with_sample_order_pair_comb2 <- renderUI({
    if(!is.null(with_pre_integrated_additional2())){
      selectInput(inputId = "with_sample_order_pair_comb2","Sample order (green):",
                  choices =  names(with_pre_integrated_additional2()),
                  selected = names(with_pre_integrated_additional2()),multiple = TRUE)
    }
  })
  output$with_sample_order_pair_comb3 <- renderUI({
    if(!is.null(with_pre_integrated_additional3())){
      selectInput(inputId = "with_sample_order_pair_comb3","Sample order (purple):",
                  choices =  names(with_pre_integrated_additional3()),
                  selected = names(with_pre_integrated_additional3()),multiple = TRUE)
    }
  })
  output$sample_order_pair_comb1 <- renderUI({
    if(!is.null(pre_integrated_additional1())){
      selectInput(inputId = "sample_order_pair_comb1","Sample order (blue):",
                  choices =  names(pre_integrated_additional1()),
                  selected = names(pre_integrated_additional1()),multiple = TRUE)
    }
  })
  output$sample_order_pair_comb2 <- renderUI({
    if(!is.null(pre_integrated_additional2())){
      selectInput(inputId = "sample_order_pair_comb2","Sample order (green):",
                  choices =  names(pre_integrated_additional2()),
                  selected = names(pre_integrated_additional2()),multiple = TRUE)
    }
  })
  output$sample_order_pair_comb3 <- renderUI({
    if(!is.null(pre_integrated_additional3())){
      selectInput(inputId = "sample_order_pair_comb3","Sample order (purple):",
                  choices =  names(pre_integrated_additional3()),
                  selected = names(pre_integrated_additional3()),multiple = TRUE)
    }
  })
  
  
  output$file1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("file1",label=NULL,choices = bw,multiple=T)
    }else{
      fileInput("file1",NULL,
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$file_bam <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bam <- c(str_subset(list,".bam"))
      selectInput("file_bam",label=NULL,choices = bam,multiple=T)
    }else{
      fileInput("file_bam",NULL,
                accept = c("bam"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$peak_call_file1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bed <- c(str_subset(list,".narrowPeak"),str_subset(list,".bed"))
      selectInput("peak_call_file1",label=NULL,choices = bed,multiple=T)
    }else{
      fileInput("peak_call_file1",NULL,
                accept = c("bed","narrowPeak"),
                multiple = TRUE,
                width = "80%")
    }
  })
  
  updateCounter_createCount <- reactiveValues(i = 0)
  
  observe({
    input$createcountButton
    isolate({
      updateCounter_createCount$i <- updateCounter_createCount$i + 1
    })
  })
  
  
  #Restart
  observeEvent(bws_count(), {
    isolate(updateCounter_createCount$i == 0)
    updateCounter_createCount <<- reactiveValues(i = 0)
  }) 
  observeEvent(bws(), {
    isolate(updateCounter_createCount$i == 0)
    updateCounter_createCount <<- reactiveValues(i = 0)
  }) 
  observeEvent(peak_call_files(), {
    isolate(updateCounter_createCount$i == 0)
    updateCounter_createCount <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$upstream, {
    isolate(updateCounter_createCount$i == 0)
    updateCounter_createCount <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$downstream, {
    isolate(updateCounter_createCount$i == 0)
    updateCounter_createCount <<- reactiveValues(i = 0)
  }) 
  # pair-wise ------------------------------------------------------------------------------
  org1 <- reactive({
    return(org(Species = input$Species))
  })
  txdb <- reactive({
    if(input$Species != "not selected"){
      return(txdb_function(Species = input$Species))
    }
  })
  promoter_region <- reactive({
    if(input$Genomic_region == "Promoter"){
      if(input$Species != "not selected"){
        return(promoter(txdb(),upstream = input$upstream, downstream = input$downstream,filter = input$pair_filter))
      }
    }else {
      if(input$data_file_type == "Row1_count"){
        data <- range_changer(bw_count())
        data2 <- with(data, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
        return(data2)
      }else{
        return(promoter(upstream = input$upstream, downstream = input$downstream,
                        input_type = "Genome-wide",files =peak_call_files(),filter = input$pair_filter))
      }
    }
  })
  gene_position <- reactive({
    if(input$Species != "not selected"){
      return(genes(txdb()))
    }
  })
  
  bws <- reactive({
    if(input$data_file_type == "Row1" || input$data_file_type == "Row1_count"){
      files<-NULL
      if(input$goButton > 0 ){
        files<-list()
        files[["A_1"]] <- "data/bigwig/A_1.BigWig"
        files[["A_2"]] <- "data/bigwig/A_2.BigWig"
        files[["B_1"]] <- "data/bigwig/B_1.BigWig"
        files[["B_2"]] <- "data/bigwig/B_2.BigWig"
      }
      if(!is.null(input$file1)){
        if(length(list.files("./Volume/")) > 0){
          files <- input$file1
          names(files) <- bigwig_name(input$file1)
        }else{
          files<-c()
          name<-c()
          for(nr in 1:length(input$file1[, 1])){
            file <- input$file1[[nr, 'datapath']]
            name <- c(name, bigwig_name(input$file1[nr,]$name))
            files <- c(files,file)
          }
          names(files)<-name
        }
      }
      return(files)
    }
  })
  bws_order<-reactive({
    return(bws_ordering(bws=bws(),sample_order=input$sample_order,additional=FALSE))
  })
  bws_count <- reactive({
    tmp <- input$file1_count$datapath
    if(is.null(input$file1_count) && input$goButton > 0 )  tmp = "data/bws_count.txt"
    if(!is.null(tmp)) return(read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = ""))
  })
  observeEvent(input$goButton,({
    if(input$data_file_type == "Row2"){
      updateRadioButtons(session,'Pair_or_single','Sequence type:',
                         c('Paired-end'="Paired-end",
                           'Single-end'="Single-end"),selected = "Single-end")
    }
  }))
  bam <- reactive({
    if(input$data_file_type == "Row2"){
      if(is.null(input$file_bam)){
        if(input$goButton > 0 ){
          df<-c()
          df["A_1"] <- "data/bam/A_1.bam"
          df["A_2"] <- "data/bam/A_2.bam"
          df["B_1"] <- "data/bam/B_1.bam"
          df["B_2"] <- "data/bam/B_2.bam"
          return(df)
        }
      }else{
        if(length(list.files("./Volume/")) > 0){
          files <- input$file_bam
          name <- gsub(".+\\/","",input$file_bam)
          names(files) <- gsub("\\..+$", "", name)
        }else{
          files<-c()
          name<-c()
          for(nr in 1:length(input$file_bam[, 1])){
            file <- input$file_bam[[nr, 'datapath']]
            name <- c(name, gsub("\\..+$", "", input$file_bam[nr,]$name))
            files <- c(files,file)
          }
          names(files)<-name
        }
        return(files)
      }
    }
  })
  
  peak_call_files <- reactive({
    if(input$Genomic_region == "Genome-wide" && input$data_file_type != "Row1_count"){
      if(is.null(input$peak_call_file1)){
        if(input$goButton > 0 ){
          df <- list()
          df[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
          df[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
          df[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
          df[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
          name <- c("A_1","A_2","B_1","B_2")
          files2 <- lapply(df, GetGRanges, simple = TRUE)
          names(files2)<-name
          return(files2)
        }
        return(NULL)
      }else{
        if(length(list.files("./Volume/")) > 0){
          files <- input$peak_call_file1
          files2 <- lapply(files, GetGRanges, simple = TRUE)
          names(files2) <- bed_name(input$peak_call_file1)
        }else{
          files<-c()
          name<-c()
          for(nr in 1:length(input$peak_call_file1[, 1])){
            file <- input$peak_call_file1[[nr, 'datapath']]
            name <- c(name, bed_name(input$peak_call_file1[nr,]$name))
            files <- c(files,file)
          }
          files2 <- lapply(files, GetGRanges, simple = TRUE)
          names(files2)<-name
        }
        return(files2)
      }
    }
  })
  
  bw_files <- reactive({
    return(BigWigFileList(bws()))
  })
  
  output$input_bw_files <- DT::renderDataTable({
    if(input$data_file_type == "Row1" || input$data_file_type == "Row1_count"){
      uploaded_files = names(bw_files())
      as.data.frame(uploaded_files)
    }else{
      uploaded_files = names(bam())
      as.data.frame(uploaded_files) 
    }
  })
  
  
  output$input_peak_call_files <- DT::renderDataTable({
    if(input$Genomic_region == "Genome-wide"){
      uploaded_files = names(peak_call_files())
      as.data.frame(uploaded_files)
    }
  })
  
  pre_bw_count <- reactive({
    if(input$data_file_type == "Row1_count" && !is.null(bws_count())){
      return(bws_count())
    }
    if(input$createcountButton > 0 && updateCounter_createCount$i > 0){
      if(input$data_file_type == "Row1" && !is.null(bw_files())){
        if(input$Genomic_region == "Promoter"){
          if(input$Species != "not selected"){
            return(Bigwig2count(bw = bw_files(),promoter_region(),
                                Species = input$Species,input_type =input$Genomic_region)) 
          }
        }else{
          if(!is.null(peak_call_files())){
            return(Bigwig2count(bw = bw_files(),promoter_region(),input_type =input$Genomic_region))
          }
        }
      }
      if(input$Species != "not selected" && input$data_file_type == "Row2" && 
         !is.null(bam())){
        withProgress(message = "converting bam to gene count data",{
          switch(input$Pair_or_single,
                 "Paired-end" = seq <- TRUE,
                 "Single-end" = seq <- FALSE)
          regionsToCount <- data.frame(GeneID = paste("ID", seqnames(promoter_region()), 
                                                      start(promoter_region()), end(promoter_region()), sep = "_"), Chr = seqnames(promoter_region()), 
                                       Start = start(promoter_region()), End = end(promoter_region()), Strand = strand(promoter_region()))
          a <- as.data.frame(promoter_region())
          Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
          
          if(input$Genomic_region == "Promoter"){
            regionsToCount <- data.frame(GeneID = as.data.frame(promoter_region())$gene_id, Chr = seqnames(promoter_region()), 
                                         Start = start(promoter_region()), End = end(promoter_region()), Strand = strand(promoter_region()))
            count <- featureCounts(bam(),annot.ext = regionsToCount,isPairedEnd = seq,
                                   countMultiMappingReads =F)
            count <- count$counts
            if(!str_detect(rownames(count)[1], "FBgn")){
              my.symbols <- rownames(count)
              or <- org(input$Species)
              gene_IDs<-AnnotationDbi::select(or,keys = my.symbols,
                                              keytype = "ENTREZID",
                                              columns = c("ENTREZID","SYMBOL"))
              colnames(gene_IDs) <- c("Row.names","SYMBOL")
              gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
              gene_IDs <- data.frame(SYMBOL = gene_IDs$SYMBOL, row.names = gene_IDs$Row.names)
              data2 <- merge(gene_IDs,count, by=0)
              rownames(data2) <- data2$SYMBOL
              count <- data2[,-1:-2]
            }
            colnames(count) <- names(bam())
            print(count)
            return(count)
          }else{
            if(!is.null(peak_call_files())){
              count <- featureCounts(bam(),annot.ext = regionsToCount,isPairedEnd = seq,
                                     countMultiMappingReads =F)$counts
              rownames(count) <- Row.name
              colnames(count) <- names(bam())
              return(count)
            }
          }
        })
      }
    }
  })
  bw_count <- reactive({
    count <- pre_bw_count()
    order <- input$sample_order
    if(!is.null(count)){
      if(dim(count)[2] == length(order)){
        data <- count[,order]
        return(data)
      }
    }
  })
  output$raw_count_table <- DT::renderDataTable({
    bw_count()
  })
  
  output$download_raw_count_table <- downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("bigwig2count-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("bigwig2count-genomeWide",".txt")
      }},
    content = function(file){write.table(bw_count(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_filtered_peakcall_bed <- downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("filtered_marged-promoter(", -input$upstream,"-",input$downstream,").bed")
      }else{
        paste0("filtered_marged-genomeWide",".bed")
      }},
    content = function(file){write.table(as.data.frame(promoter_region()), file, row.names = F,col.names=F, sep = "\t", quote = F)}
  )
  
  ##pair-wise deseq2--------
  collist_bw_pair <- reactive({
    count <- bw_count()
    collist <- gsub("\\_.+$", "", colnames(count))
    return(collist)
  })
  dds <- reactive({
    count <- bw_count()
    collist <- collist_bw_pair()
    if(length(unique(collist)) != 2) validate(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 2 conditions."))
    withProgress(message = "DESeq2",{
      group <- data.frame(con = factor(collist))
      dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
      dds$con <- factor(dds$con, levels = unique(collist))
      dds <- DESeq(dds)
      incProgress(1)
    })
    return(dds)
  })
  output$regression_mode <- renderUI({
    if(input$data_file_type == "Row1"){
      radioButtons("regression_mode","Regression",
                   c('least squares'=FALSE,
                     'robust'=TRUE
                   ), selected = TRUE)
    }
  })
  dds_limma <- reactive({
    withProgress(message = "Detecting differential accessible region",{
      count <- log(bw_count() + 1,2)
      if(!is.null(input$regression_mode)){
        collist <- collist_bw_pair()
        print(collist)
        if(length(unique(collist)) != 2) validate(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 2 conditions."))
        collist <- gsub(" ", ".", collist)
        collist <- gsub("\\_.+$", "", colnames(count))
        colnames(count) <- gsub("\\-","\\_",colnames(count))
        collist <- gsub("\\-","\\_",collist)
        
        collist<- factor(collist, levels = unique(collist),ordered = TRUE)
        eset = new("ExpressionSet", exprs=as.matrix(count))
        design <- model.matrix(~0+collist)
        colnames(design) <- factor(unique(collist),levels = unique(collist))
        fit <- lmFit(eset, design)
        comparisons <-  paste(unique(collist)[1],"-",unique(collist)[2],sep="")
        cont.matrix <- makeContrasts(contrasts=comparisons, levels=design)
        fit <- contrasts.fit(fit, cont.matrix)
        fit2 <- try(eBayes(fit,trend = TRUE,robust = input$regression_mode))
        if(class(fit2) == "try-error") fit2 <- eBayes(fit)
        result =topTable(fit2, number = 1e12)
        colnames(result) <- c("log2FoldChange","baseMean","t","pval","padj","B")
        return(result)
      }
      incProgress(1)
    })
  })
  
  
  deg_result <- reactive({
    if(is.null(bw_count())){
      return(NULL)
    }else{
      if(input$data_file_type == "Row2"){
        count <- bw_count()
        collist <- collist_bw_pair()
        contrast <- c("con", unique(collist))
        dds <- dds()
        res <- results(dds,  contrast = contrast)
        res <- as.data.frame(res)
      }
      if(input$data_file_type == "Row1"){
        res <- dds_limma()
      }
      if(input$data_file_type == "Row1_count"){
        if(input$count_file_type == "Norm") {
          res <- dds_limma()
        }else{
          count <- bw_count()
          collist <- collist_bw_pair()
          dds <- dds()
          contrast <- c("con", unique(collist))
          res <- results(dds,  contrast = contrast)
          res <- as.data.frame(res)
        }
      }
      return(res)
    }
  })
  
  deg_norm_count <- reactive({
    if(is.null(bw_count())){
      return(NULL)
    }else{
      count <- bw_count()
      if(input$data_file_type == "Row2"){
        collist <- collist_bw_pair()
        group <- data.frame(con = factor(collist))
        dds <- dds()
        contrast <- c("con", unique(collist))
        normalized_counts <- counts(dds, normalized=TRUE)
      }
      if(input$data_file_type == "Row1"){
        normalized_counts <- count
      }
      if(input$data_file_type == "Row1_count"){
        if(input$count_file_type == "Norm") {
          normalized_counts <- count
        }else{
          collist <- collist_bw_pair()
          group <- data.frame(con = factor(collist))
          dds <- dds()
          contrast <- c("con", unique(collist))
          normalized_counts <- counts(dds, normalized=TRUE)
        }
      }
      return(normalized_counts)
    }
  })
  
  observeEvent(bws(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="raw_count_panel")
  }))
  observeEvent(bam(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="raw_count_panel")
  }))
  observeEvent(peak_call_files(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="peak_call_files_panel")
  }))
  observeEvent(bws_count(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="raw_count_panel")
  }))
  
  data_degcount <- reactive({
    data <- deg_result()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- collist_bw_pair()
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      Row.names <- NULL
      log2FoldChange <- NULL
      value <- NULL
      data <- merge(data,count, by=0)
      data <- dplyr::filter(data, apply(data[,8:(7 + Cond_1 + Cond_2)],1,mean) > input$basemean)
      data$log2FoldChange <- -1 * data$log2FoldChange
      if(input$Genomic_region == "Promoter"){
        if(input$Species != "not selected"){
          my.symbols <- data$Row.names
          gene_IDs<-id_convert(my.symbols,input$Species,type="SYMBOL_double")
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data <- merge(data, gene_IDs, by="Row.names")
        }
      }
      genenames <- as.vector(data$Row.names)
      return(data)
    }
  })
  
  data_degcount2 <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- collist_bw_pair()
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- dplyr::filter(data, abs(data$log2FoldChange) > log2(input$fc))
      if(nrow(data2) != 0){
        data2$group <- paste0(unique(collist)[2],"_high")
        data2$group[data2$log2FoldChange < 0] <- paste0(unique(collist)[1],"_high")
        data3 <- dplyr::filter(data2, abs(data2$padj) < input$fdr)
        return(data3)
      }else{return(NULL)}
    }
  })
  
  deg_result_anno <- reactive({
    if(input$Species == "not selected") validate("Please select 'Species'.")
    if(!is.null(deg_result())){
      withProgress(message = "preparing annotation",{
        if(input$data_file_type == "Row1_count" && input$Genomic_region == "Genome-wide"){
          data <- range_changer(bws_count())
          data <- with(data, GRanges(seqnames = chr,ranges = IRanges(start=start,
                                                                     end=end)))
        }else data <- promoter_region()
        data2 <- ChIPseeker::annotatePeak(data, TxDb= txdb())
        return(data2)
      })
    }
  })
  deg_result_anno2 <- reactive({
    if(!is.null(deg_result_anno())){
      data <- as.data.frame(ChIPseeker::as.GRanges(deg_result_anno()))
      Row.name <- paste0(data$seqnames,":",data$start,"-",data$end)
      data$locus <- Row.name
      my.symbols <- data$geneId
      gene_IDs<-id_convert(my.symbols,input$Species,type="ENTREZID")
      colnames(gene_IDs) <- c("geneId","NearestGene")
      data <- merge(data, gene_IDs, by="geneId")
      data <- data %>% distinct(locus, .keep_all = T)
      rownames(data)<-data$locus
      return(data)
    }
  })
  
  deg_result2 <- reactive({
    if(!is.null(deg_result_anno2())){
      if(input$Genomic_region == "Genome-wide"){
        anno <- deg_result_anno2()
        anno <- dplyr::arrange(anno, locus)
        data <- data.frame(NearestGene = anno$NearestGene,
                           locus = anno$locus)
        deg_result <- deg_result()
        deg_result$locus <- rownames(deg_result())
        deg_result <- dplyr::arrange(deg_result, locus)
        data2 <- merge(data,deg_result,by="locus")
        rownames(data2) <- data2$locus
        data2 <- data2[, - which(colnames(data2) == "locus")]
        return(data2)
      }
    }
  })
  pre_DEG_result <- reactive({
    if(!is.null(deg_result())){
      if(dim(brush_info())[1] == 0){
        if(input$Genomic_region == "Promoter" || input$Species == "not selected"){
          deg_result()
        }else{
          deg_result2()
        }
      }else{
        if(input$Genomic_region == "Promoter" || input$Species == "not selected"){
          deg_result() %>% 
            dplyr::filter(rownames(.) %in% brush_info()$Row.names)
        }else{
          deg_result2()  %>%
            dplyr::filter(rownames(.) %in% brush_info()$Row.names)
        }
      }
    }
  })
  
  output$DEG_result <- DT::renderDT({
    if(!is.null(deg_result())){
      pre_DEG_result() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  
  output$download_pair_DEG_result = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("DAR_result-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("DAR_result-genomeWide",".txt")
      }},
    content = function(file){
      if(input$Genomic_region == "Promoter"){
        write.table(deg_result(), file, row.names = T, sep = "\t", quote = F)
      }else{
        if(input$Species != "not selected"){
          write.table(deg_result2(), file, row.names = T, sep = "\t", quote = F)
        }else write.table(deg_result(), file, row.names = T, sep = "\t", quote = F)
      }
    }
  )
  output$download_pair_norm_count = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("DESeq2_normCount-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("DESeq2_normCount-genomeWide", ".txt")
      }},
    content = function(file){write.table(deg_norm_count(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  
  # pair-wise PCA ------------------------------------------------------------------------------
  output$download_pair_PCA = downloadHandler(
    filename = "PCA-MDS-dendrogram.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(PCAplot(data = deg_norm_count(),plot=TRUE,legend=input$PCA_legend))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$PCA <- renderPlot({
    withProgress(message = "preparing result table",{
      if(is.null(deg_norm_count())){
        return(NULL)
      }else{
        print(PCAplot(data = deg_norm_count(),plot=TRUE,legend=input$PCA_legend))
      }
    })
  })
  
  output$pair_PCA_data <- DT::renderDataTable({
    PCAplot(data = deg_norm_count(),plot=FALSE)
  })
  
  output$download_pair_PCA_table = downloadHandler(
    filename ="PCA_table.txt",
    content = function(file){write.table( PCAplot(data = deg_norm_count(),plot=FALSE), 
                                          file, row.names = T, sep = "\t", quote = F)}
  )
  ###GOI----------------------
  output$volcano_x <- renderUI({
    if(!is.null(data_degcount())){
      withProgress(message = "Preparing volcano plot",{
        data <- as.data.frame(data_degcount())
        data <- na.omit(data)
        min <- floor(min(data$log2FoldChange))
        max <- ceiling(max(data$log2FoldChange))
        sliderInput("xrange","X_axis range:",min = min-1,
                    max=max+1, step = 1,
                    value = c(min, max))
      })
    }
  })
  output$volcano_y <- renderUI({
    if(!is.null(data_degcount())){
      data <- as.data.frame(data_degcount())
      data$padj[data$padj == 0] <- 10^(-300)
      data <- na.omit(data)
      max <- ceiling(max(-log10(data$padj)))
      sliderInput("yrange","Y_axis range:",min = 0, max= max+1, step = 1,
                  value = max)
    }
  })
  brush_info <- reactive({
    if(!is.null(input$xrange) && !is.null(input$yrange) && !is.null(data_degcount())){
      data <- as.data.frame(data_degcount())
      data$padj[data$padj == 0] <- 10^(-300)
      data$minusLog10padj<- -log10(data$padj)
      return(brushedPoints(data, input$plot1_brush,xvar = "log2FoldChange",yvar="minusLog10padj"))
    }else return(data.frame())
  })
  
  pair_volcano <- reactive({
    if(!is.null(input$xrange) && !is.null(input$yrange)){
      withProgress(message = "volcano plot",{
        data <- as.data.frame(data_degcount())
        count <- deg_norm_count()
        collist <- collist_bw_pair()
        data <- data %>% dplyr::filter(!is.na(padj)) 
        if(dim(brush_info())[1]!=0){
          label_data <- brush_info()$Row.names
        }else{
          if(!is.null(input$DEG_result_rows_selected)){
            label_data <- rownames(pre_DEG_result()[input$DEG_result_rows_selected,])
          }else label_data <- NULL
        }
        data$color <- "NS"
        data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr] <- paste(paste0(unique(collist)[1],"_high:"), length(data$Row.names[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr]))
        data$color[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr] <- paste(paste0(unique(collist)[2],"_high"),length(data$Row.names[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr]))
        data$padj[data$padj <= 10^(-300)] <- 10^(-300)
        if(!is.null(label_data)) {
          for(name in label_data){
            data$color[data$Row.names == name] <- "GOI"
          }
          Color <- c("blue","green","darkgray","red")
          data$color <- factor(data$color, levels = c(paste(paste0(unique(collist)[1],"_high:"), length(data$Row.names[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr])),
                                                      "GOI","NS", paste(paste0(unique(collist)[2],"_high"),length(data$Row.names[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr]))))
          if(length(data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr]) == 0) Color <- c("green","darkgray","red")
        }else{
          Color <- c("blue","darkgray","red")
          data$color <- factor(data$color, levels = c(paste(paste0(unique(collist)[1],"_high:"), length(data$Row.names[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr])),
                                                      "NS", paste(paste0(unique(collist)[2],"_high"),length(data$Row.names[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr]))))
          if(length(data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr]) == 0) Color <- c("darkgray","red")
        }
        data$minusLog10padj<--log10(data$padj)
        v <- ggplot(data, aes(x = log2FoldChange, y = minusLog10padj)) + ggrastr::geom_point_rast(aes(color = color),size = 0.4)
        v <- v  + geom_vline(xintercept = c(-log2(input$fc), log2(input$fc)), linetype = c(2, 2), color = c("black", "black")) +
          geom_hline(yintercept = c(-log10(input$fdr)), linetype = 2, color = c("black"))
        v <- v +theme_bw()+ scale_color_manual(values = Color)+
          theme(legend.position = "top" , legend.title = element_blank(),
                axis.text.x= ggplot2::element_text(size = 12),
                axis.text.y= ggplot2::element_text(size = 12),
                text = ggplot2::element_text(size = 12),
                title = ggplot2::element_text(size = 12)) +
          xlab("log2 fold change") + ylab("-log10(padj)") +
          xlim(input$xrange)+
          ylim(c(0, input$yrange))
        if(!is.null(label_data)) {
          v <- v + geom_point(data=dplyr::filter(data, Row.names == label_data),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, Row.names == label_data), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), fontface = "bold.italic")
        }
        return(v)
        incProgress(1)
      })
    }else return(NULL)
  })
  
  output$volcano1 <- renderPlot({
    if(!is.null(input$xrange)){
      if(is.null(bw_count())){
        return(NULL)
      }else{
        pair_volcano()
      }
    }
  })
  
  
  volcano_heatmap <- reactive({
    vol <- as.grob(pair_volcano())
    return(plot_grid(vol, rel_widths = c(2, 1)))
  })
  
  output$download_pair_volcano = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("Volcano-promoter(", -input$upstream,"-",input$downstream,").pdf")
      }else{
        paste0("Volcano-genomeWide",".pdf")
      }},
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 4
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_volcano())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##trackplot----------
  ref <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species))
    return(ref)
  })
  
  
  
  observeEvent(input$DEG_result_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_pair_DEG", open="Trackplot_panel")
  }))
  observeEvent(input$DEG_result_rows_selected,({
    if(input$Species != "not selected"){
      if(!is.null(goi_gene_position()) && !is.null(goi_promoter_position())){
        y <- goi_promoter_position()
        gene_position <- goi_gene_position()
        start_position <- min(c(y$start,gene_position$start))
        end_position <- max(c(y$end,gene_position$end))
        updateSliderInput(session,"igv_uprange","Range (x-axis):",
                          value = c(start_position,end_position),
                          step = 100, min = start_position - 10000, max = end_position + 10000)
      }
    }
  }))
  
  output$igv_uprange <- renderUI({
    if(input$data_file_type != "Row2"){
      if(is.null(bws())) validate("Bigwig files are required.")
    }
    if(!is.null(goi_gene_position()) && !is.null(goi_promoter_position())){
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("igv_uprange","Range (x-axis):",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$igv_ylim <- renderUI({
    numericInput("igv_ylim","Max peak intensity (y-axis):", value = 1, min = 0)
  })
  
  gtrack <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=0.8)
      return(gtrack)
    })
  })
  
  goi_promoter_position<- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      library(Gviz)
      if(input$Genomic_region == "Promoter"){
        label_data <- rownames(pre_DEG_result()[input$DEG_result_rows_selected,])
        gene_IDs<-id_convert(label_data,input$Species,type="SYMBOL_single")
        y <- as.data.frame(subset(promoter_region(), gene_id %in% gene_IDs))
      }else{
        label_data <- rownames(pre_DEG_result()[input$DEG_result_rows_selected,])
        y <- dplyr::filter(deg_result_anno2(),locus == label_data)
      }
      return(y)
    }
  })
  
  goi_gene_position <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      if(input$Genomic_region == "Promoter"){
        label_data <- rownames(pre_DEG_result()[input$DEG_result_rows_selected,])
      }else{
        label_data <- pre_DEG_result()[input$DEG_result_rows_selected,]$NearestGene
      }
      gene_IDs<- id_convert(label_data,input$Species,type="SYMBOL_single")
      gene_position <- as.data.frame(subset(gene_position(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  output$trackplot_additional <- renderUI({
    if(!is.null(bws())){
      if(length(list.files("./Volume/")) > 0){
        list <- list.files("Volume",full.names = T,recursive=T)
        bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
        selectInput("trackplot_additional1",label="Select additional bigwig files",choices = bw,multiple=T)
      }else{
        fileInput("trackplot_additional1",
                  "Select additional bigwig files",
                  accept = c("bw","BigWig"),
                  multiple = TRUE,
                  width = "80%")
      }
    }
  })
  pre_track_additional_files <-reactive({
    if(!is.null(input$trackplot_additional1)){
      if(length(list.files("./Volume/")) > 0){
        files <- input$trackplot_additional1
        name <- gsub(".+\\/","",input$trackplot_additional1)
        names(files) <- gsub("\\..+$", "", name)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$trackplot_additional1[, 1])){
          file <- input$trackplot_additional1[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$trackplot_additional1[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
      }
      return(files)
    }
  })
  track_additional_files <- reactive({
    return(bws_ordering(bws=pre_track_additional_files(),sample_order=input$sample_order_pair_track))
  })
  
  data_track <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      return(data_trac(y=goi_promoter_position(),gene_position=goi_gene_position(),
                       gen=ref(),txdb=txdb(),org=org1(),filetype=input$data_file_type,
                       bw_files=bws_order(),bam_files=bam(),
                       track_additional_files=track_additional_files()))
    }
  })
  
  highlight_trackplot <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      library(Gviz)
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      chr <- gene_position$seqnames
      df <- data_track()
      if(y$start < input$igv_uprange[2] && y$end > input$igv_uprange[1]){
        withProgress(message = "Highlighting the selected region",{
          ht <- HighlightTrack(trackList = df,
                               start = y$start, width = y$width,
                               chromosome = chr, alpha=0.5)
        })
      }else ht <- NULL
      return(ht)
    }
  })
  goi_trackplot <- reactive({
    if(!is.null(input$DEG_result_rows_selected) &&
       !is.null(goi_gene_position()) && 
       !is.null(input$igv_uprange)){
      library(Gviz)
      if(!is.null(highlight_trackplot())){
        plot<- plotTracks(list(gtrack(), highlight_trackplot()),
                          from = input$igv_uprange[1], 
                          to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                          type="hist",add=TRUE)
      }else{
        df <- data_track()
        df[["gtrack"]] <- gtrack()
        plot<- plotTracks(df,
                          from = input$igv_uprange[1], 
                          to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                          type="hist",add=TRUE)
      }
      return(plot)
    }
  })
  output$trackplot_goi <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$DEG_result_rows_selected)  &&
         !input$Species == "not selected"){
        goi_trackplot()
      }
    })
  })
  
  output$download_pair_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        print(class(goi_trackplot()))
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(highlight_trackplot())){
          plotTracks(list(gtrack(), highlight_trackplot()),
                     from = input$igv_uprange[1], 
                     to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                     type="hist",add=TRUE)
        }else{
          df <- data_track()
          df[["gtrack"]] <- gtrack()
          plotTracks(df,
                     from = input$igv_uprange[1], 
                     to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                     type="hist",add=TRUE)
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  ##Functional enrichment analysis-------------
  output$Gene_set <- renderUI({
    selectInput('Gene_set', 'Gene Set', gene_set_list)
  })
  data_degcount_up <- reactive({
    if(!is.null(deg_result())){
      if(input$Genomic_region == "Genome-wide"){
        data2 <- deg_result2()
        up_all <- dplyr::filter(data2, log2FoldChange < -log2(input$fc) & padj < input$fdr)
      }else{
        data <- deg_result() %>% dplyr::filter(log2FoldChange < -log2(input$fc) & padj < input$fdr)
        if(dim(data)[1] != 0){
        up <- symbol2gene_id(data,org1())
        up2 <- subset(promoter_region(), gene_id %in% up$gene_id) %>% as.data.frame()
        up_all <- data.frame(chr = up2$seqnames, start = up2$start, end = up2$end)
        }else up_all <- NULL
      }
      return(up_all)
    }
  })
  data_degcount_up_bed <-reactive({
    if(input$Genomic_region == "Genome-wide"){
      if(input$Species != "not selected"){
        data <- range_changer(data_degcount_up())
        data2 <- data.frame(chr = data$chr, start = data$start, end = data$end) 
      }else{
        data <- deg_result() %>% dplyr::filter(log2FoldChange < -log2(input$fc) & padj < input$fdr)
        data <- range_changer(data)
        data2 <- data.frame(chr = data$chr, start = data$start, end = data$end) 
      }
    }else{
      data <- deg_result() %>% dplyr::filter(log2FoldChange < -log2(input$fc) & padj < input$fdr)
      if(dim(data)[1] != 0){
      up <- symbol2gene_id(data,org1())
      up2 <- subset(promoter_region(), gene_id %in% up$gene_id) %>% as.data.frame()
      data2 <- data.frame(chr = up2$seqnames, start = up2$start, end = up2$end)
      }else data2 <-NULL
    }
    return(data2)
  })
  data_degcount_down_bed <-reactive({
    if(input$Genomic_region == "Genome-wide"){
      if(input$Species != "not selected"){
        data <- range_changer(data_degcount_down())
        data2 <- data.frame(chr = data$chr, start = data$start, end = data$end)
      }else{
        data <- deg_result()
        data$color <- "NS"
        data$color[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr] <- "down"
        data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr] <- "up"
        data <- data %>% dplyr::filter(color == "down")
        data <- range_changer(data)
        data2 <- data.frame(chr = data$chr, start = data$start, end = data$end) 
      }
    }else{
      data <- deg_result() %>% dplyr::filter(log2FoldChange > log2(input$fc) & padj < input$fdr)
      if(dim(data)[1] != 0){
      down <- symbol2gene_id(data,org1())
      down2 <- subset(promoter_region(), gene_id %in% down$gene_id) %>% as.data.frame()
      data2 <- data.frame(chr = down2$seqnames, start = down2$start, end = down2$end)
      }else data2 <- NULL
    }
    return(data2)
  })
  output$DEG_up_result <- DT::renderDT({
    data_degcount_up_bed()
  })
  output$DEG_down_result <- DT::renderDT({
    data_degcount_down_bed()
  })
  output$download_pair_DEG_up_result = downloadHandler(
    filename = function() {
      paste0(unique(collist_bw_pair())[2],"_high.bed")
    },
    content = function(file){write.table(data_degcount_up_bed(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  output$download_pair_DEG_down_result = downloadHandler(
    filename = function() {
      paste0(unique(collist_bw_pair())[1],"_high.bed")
    },
    content = function(file){write.table(data_degcount_down_bed(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  
  data_degcount_down <- reactive({
    if(!is.null(deg_result())){
      if(input$Genomic_region == "Genome-wide"){
        data2 <- deg_result2()
        down_all <- dplyr::filter(data2, log2FoldChange > log2(input$fc) & padj < input$fdr)
      }else{
        data <- deg_result() %>% dplyr::filter(log2FoldChange > log2(input$fc) & padj < input$fdr)
        if(dim(data)[1] != 0){
        down <- symbol2gene_id(data,org1())
        down2 <- subset(promoter_region(), gene_id %in% down$gene_id) %>% as.data.frame()
        down_all <- data.frame(chr = down2$seqnames, start = down2$start, end = down2$end)
        }else down_all <- NULL
      }
      return(down_all)
    }
  })
  Hallmark_set <- reactive({
    if(!is.null(input$Gene_set)){
      return(GeneList_for_enrichment(Species = input$Species, Gene_set = input$Gene_set, org = org1()))
    }
  })
  
  data_degcount2_anno <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      data <- data_degcount2()
      rownames(data) <- data$Row.names
      data <- data[, -1]
      anno <- deg_result_anno2()
      colnames(anno)[1] <- "ENTREZID"
      data2 <- merge(data,anno,by=0)
      return(data2)
    }
  })
  
  enrichment_enricher <- reactive({
    data3 <- data_degcount2()
    collist <- collist_bw_pair()
    print("enrichment_enricher start")
    if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(data3)){
      if(input$Genomic_region == "Promoter"){
        withProgress(message = "enrichment analysis",{
          H_t2g <- Hallmark_set()
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          em_up <- try(clusterProfiler::enricher(dplyr::filter(data3, group == paste0(unique(collist)[2],"_high"))$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
          em_down <- try(clusterProfiler::enricher(dplyr::filter(data3, group == paste0(unique(collist)[1],"_high"))$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
          df <- list()
          df[[paste0(unique(collist)[2],"_high")]] <- em_up
          df[[paste0(unique(collist)[1],"_high")]] <- em_down
          for(name in names(df)){
            if (length(as.data.frame(df[[name]])$ID) == 0) {
              df[[name]] <- NULL
            } else{
              df[[name]] <- clusterProfiler::setReadable(df[[name]], org1(), 'ENTREZID')
            }
          }
          incProgress(1)
          return(df)
        })
      }else{ 
        data3 <- deg_result_anno2()
        H_t2g <- Hallmark_set()
        source <- ref_for_GREAT(input$Species)
        data_degcount_up2 <- dplyr::filter(data3, locus %in% rownames(data_degcount_up()))
        data_degcount_up3 <- with(data_degcount_up2, GRanges(seqnames = seqnames, 
                                                             ranges = IRanges(start,end),
                                                             ENTREZID = geneId))
        res_up = rGREAT::great(gr = data_degcount_up3,gene_sets = gene_list_for_enrichment_genome(H_t2g,input$Species),source)
        data_degcount_down2 <- dplyr::filter(data3, locus %in% rownames(data_degcount_down()))
        data_degcount_down3 <- with(data_degcount_down2, GRanges(seqnames = seqnames, 
                                                                 ranges = IRanges(start,end),
                                                                 ENTREZID = geneId))
        res_down = rGREAT::great(data_degcount_down3,gene_list_for_enrichment_genome(H_t2g,input$Species), source)
        df <- list()
        df[[paste0(unique(collist)[2],"_high")]] <- res_up
        df[[paste0(unique(collist)[1],"_high")]] <- res_down
        incProgress(1)
        return(df)
      }
    }else{return(NULL)}
  })
  
  enrichment_1_1 <- reactive({
    df <- enrichment_enricher()
    if(!is.null(df)){
      if(input$Genomic_region == "Promoter"){
        if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(df)){
          data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for(name in names(df)){
            if(!is.null(df[[name]])) {
              group1 <- as.data.frame(df[[name]])
            }else group1 <- NULL
            group1$Group <- name
            data <- rbind(data, group1)
          }
          if(length(data$Description) != 0) {
            data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
            data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
            return(data)
          }else return(NULL)
        }else{return(NULL)}
      }else{
        if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(df)){
          data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
          colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                              "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
          for(name in names(df)){
            if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
              group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
              group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
              group1$Group <- name
            }else group1 <- NULL
            data <- rbind(data, group1)
          }
          if(length(data$id) != 0) {
            return(data)
          }else return(NULL)
        }else{return(NULL)}
      }}
  })
  
  pair_enrich1_H <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      count <- deg_norm_count()
      df <- enrichment_enricher()
      collist <- collist_bw_pair()
      if(input$Genomic_region == "Promoter"){
        data3 <- data_degcount2()
        H_t2g <- Hallmark_set()
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        if (is.null(df[[paste0(unique(collist)[2],"_high")]]) && is.null(df[[paste0(unique(collist)[1],"_high")]]))  {
          p1 <- NULL
        } else{
          data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          for(name in names(df)){
            if(!is.null(df[[name]])) {
              group1 <- as.data.frame(df[[name]])
              group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$ENTREZID),")",sep = "")
              if (length(group1$pvalue) > 5){
                group1 <- group1[sort(group1$pvalue, decreasing = F, index=T)$ix,]
                group1 <- group1[1:5,]
              }
            }else group1 <- NULL
            data <- rbind(data, group1)
          }
          if(length(data$Description) != 0){
            data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
            data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
            if ((length(data$Description) == 0) || length(which(!is.na(unique(data$qvalue))))==0) {
              p1 <- NULL
            } else{
              data$Description <- gsub("_", " ", data$Description)
              data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
              data$x <- gsub(":","", data$x)
              data <- dplyr::arrange(data, x)
              idx <- order(data[["x"]], decreasing = FALSE)
              data$Description <- factor(data$Description,
                                         levels=rev(unique(data$Description[idx])))
              p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="qvalue",size="GeneRatio"))+
                              geom_point() +
                              scale_color_continuous(low="red", high="blue",
                                                     guide=guide_colorbar(reverse=TRUE)) +
                              scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                              scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
            }}else p1 <- NULL
        }
      }else{
        data3 <- data_degcount2_anno()
        if (is.null(df[[paste0(unique(collist)[2],"_high")]]) && is.null(df[[paste0(unique(collist)[1],"_high")]]))  {
          p1 <- NULL
        } else{
          data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
          for(name in names(df)){
            if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
              group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
              group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$group),")",sep = "")
              if (length(group1$p_adjust_hyper) > 5){
                group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
                group1 <- group1[1:5,]
              }
            }else group1 <- NULL
            data <- rbind(data, group1)
          }
          colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                                 "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
          if(length(data$Description) != 0){
            data$Description <- gsub("_", " ", data$Description)
            data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
            data$x <- gsub(":","", data$x)
            data <- dplyr::arrange(data, x)
            idx <- order(data[["x"]], decreasing = FALSE)
            data$Description <- factor(data$Description,
                                       levels=rev(unique(data$Description[idx])))
            p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                            scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
          }else p1 <- NULL
        }
      }
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment1 <- renderPlot({
    dotplot_for_output(data = bw_count(),
                       plot_genelist = pair_enrich1_H(), Gene_set = input$Gene_set, 
                       Species = input$Species)
  })
  
  output$download_pair_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$Gene_set,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H()
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  pair_enrich2 <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      if(input$Genomic_region == "Promoter"){
        data3 <- data_degcount2()
      }else data3 <- data_degcount2_anno()
      count <- deg_norm_count()
      df <- enrichment_enricher()
      collist <- collist_bw_pair()
      upgene <- data3[data3$log2FoldChange > log(input$fc, 2),]
      downgene <- data3[data3$log2FoldChange < log(1/input$fc, 2),]
      p <- list()
      for(name in names(df)){
        if(length(as.data.frame(df[[name]])$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- clusterProfiler::setReadable(df[[name]], org1(), 'ENTREZID')
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          if(name == paste0(unique(collist)[2],"_high")) genes <- upgene
          if(name == paste0(unique(collist)[1],"_high")) genes <- downgene
          geneList <- genes$log2FoldChange
          names(geneList) = as.character(genes$ENTREZID)
          p2 <- try(as.grob(cnetplot(cnet1, foldChange=geneList,
                                     cex_label_gene = 0.7, cex_label_category = 0.75,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }
        }
        p[[name]] <- p2
      }
      p <- plot_grid(p[[paste0(unique(collist)[2],"_high")]], p[[paste0(unique(collist)[1],"_high")]], nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  
  pair_enrich_table <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1()), H_t2g = Hallmark_set(), Gene_set = input$Gene_set))
  })
  
  output$pair_enrichment_result <- DT::renderDataTable({
    as.data.frame(enrichment_1_1())
  })
  
  output$download_pair_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$Gene_set,".txt")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$whichGeneSet <- renderUI({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Up_or_Down) && !is.null(as.data.frame(enrichment_1_1())$id)){
      group <- input$Up_or_Down
      set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1()), Group == group)$id)
      selectInput('Pathway_list', 'Pathway list', set_list)
    }
  })
  
  region_gene_associate <- reactive({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Pathway_list)){
      set_list <- input$Pathway_list
      df <- enrichment_enricher()
      data3 <- data_degcount2()
      data <- data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0 &&
           length(as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))$seqnames) != 0) {
          group1 <- as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))
          group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$group),")",sep = "")
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      data <- apply(data,2,as.character)
      return(data)
    }
  })
  output$region_gene_associations <- DT::renderDT({
    if(input$Genomic_region == "Genome-wide"){
      region_gene_associate()
    }
  })
  output$whichGenes <- renderUI({
    if(input$Genomic_region == "Genome-wide"){
      if(!is.null(enrichment_enricher())){
        selectInput('Up_or_Down', 'Which group', names(enrichment_enricher()))
      }
    }
  })
  region_gene_associate_plot <- reactive({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Pathway_list) && 
       !is.null(input$Up_or_Down) && !is.na(input$Up_or_Down) && input$Up_or_Down != "" && 
       !is.null(enrichment_enricher())){
      set_list <- input$Pathway_list
      df <- enrichment_enricher()
      name <- input$Up_or_Down
      if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
        res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
      }else res <- NULL
      return(res)
    }
  })
  output$region_gene_associations_plot <- renderPlot({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Pathway_list) && 
       !is.null(input$Up_or_Down) && !is.na(input$Up_or_Down) && input$Up_or_Down != "" &&  
       !is.null(enrichment_enricher())){
      if(class(try(region_gene_associate_plot())) !="try-error") region_gene_associate_plot()
    }
  })
  output$download_region_gene_associations_plot = downloadHandler(
    filename = function(){
      paste0("region_gene_associations-",input$Pathway_list,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$pair_pdf_width
        set_list <- input$Pathway_list
        df <- enrichment_enricher()
        name <- input$Up_or_Down
        pdf(file, height = pdf_height, width = pdf_width)
        rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_region_gene_associations = downloadHandler(
    filename = function() {
      paste0("region_gene_associations",input$Pathway_list,".txt")
    },
    content = function(file){write.table(region_gene_associate(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  ##peak distribution-------
  input_peak_list <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      files<-list()
      name <- c()
      if(!is.null(peak_call_files())){
        a <- peak_call_files()
        Glist <- GRangesList()
        for(name in names(a)){
          Glist[[name]] <- a[[name]]
        }
        Glist[["Consensus_region"]] <- promoter_region()
        return(Glist)
      }else{
        if(input$goButton > 0 ){
          files[["data/peakcall/A_1_peaks.narrowPeak"]] <- "data/peakcall/A_1_peaks.narrowPeak"
          files[["data/peakcall/A_2_peaks.narrowPeak"]] <- "data/peakcall/A_2_peaks.narrowPeak"
          files[["data/peakcall/B_1_peaks.narrowPeak"]] <- "data/peakcall/B_1_peaks.narrowPeak"
          files[["data/peakcall/B_2_peaks.narrowPeak"]] <- "data/peakcall/B_2_peaks.narrowPeak"
          Glist <- files2GRangelist(files)
          Glist[["Consensus_region"]] <- promoter_region()
          return(Glist)
        }else return(NULL)
      }
    }
  })
  updistribution <- reactive({
    return(ChIPpeakAnno::genomicElementDistribution(deg_peak_list()[["Up"]], 
                                                    TxDb = txdb()))
  })
  output$input_peak_distribution <- renderPlot({
    if(input$Genomic_region == "Genome-wide"){
      if(!is.null(deg_peak_list()) && input$Species != "not selected"){
        withProgress(message = "Plot peak distribution of up peak",{
          updistribution()$plot
        })
      }
    }
  })
  deg_peak_list <- reactive({
    if(!is.null(data_degcount_up()) && !is.null(data_degcount_down())){
      up <- range_changer(data_degcount_up())
      down <- range_changer(data_degcount_down())
      print(head(up))
      Glist <- GRangesList()
      up2 <- with(up, GRanges(seqnames = chr,
                              ranges = IRanges(start=start,end=end),
                              NearestGene = NearestGene,
                              log2FoldChange = log2FoldChange,
                              padj = padj
      ))
      down2 <- with(down, GRanges(seqnames = chr,
                                  ranges = IRanges(start=start,end=end),
                                  NearestGene = NearestGene,
                                  log2FoldChange = log2FoldChange,
                                  padj = padj))
      Glist[["Up"]] <- up2
      Glist[["Down"]] <- down2
      return(Glist)
    }
  })  
  downdistribution <- reactive({
    return(ChIPpeakAnno::genomicElementDistribution(deg_peak_list()[["Down"]], 
                                                    TxDb = txdb()))
  })
  
  output$deg_peak_distribution <- renderPlot({
    withProgress(message = "Plot peak distribution of down peak",{
      if(input$Genomic_region == "Genome-wide"){
        if(!is.null(deg_peak_list()) && input$Species != "not selected"){
          downdistribution()$plot
        }
      }
    })
  })
  output$up_distribution_table <- renderDataTable({
    if(input$Genomic_region == "Genome-wide"){
      if(!is.null(deg_peak_list()) && input$Species != "not selected"){
        df <- as.data.frame(updistribution()$peaks)
        df$Group <- paste0(unique(collist_bw_pair())[2],"_high")
        df
      }
    }
  })
  output$down_distribution_table <- renderDataTable({
    if(input$Genomic_region == "Genome-wide"){
      if(!is.null(deg_peak_list()) && input$Species != "not selected"){
        df <- as.data.frame(downdistribution()$peaks)
        df$Group <- paste0(unique(collist_bw_pair())[1],"_high")
        df
      }
    }
  })
  output$download_input_peak_distribution_table = downloadHandler(
    filename = function() {
      paste0(unique(collist_bw_pair())[2],"_high_peak_distribution",".txt")
    },
    content = function(file) {
      df <- as.data.frame(updistribution()$peaks)
      df$Group <- paste0(unique(collist_bw_pair())[2],"_high")
      write.table(df, file, row.names = F, sep = "\t", quote = F)
    }
  )
  output$download_down_peak_distribution_table = downloadHandler(
    filename = function() {
      paste0(unique(collist_bw_pair())[1],"_high_peak_distribution",".txt")
    },
    content = function(file) {
      df <- as.data.frame(downdistribution()$peaks)
      df$Group <- paste0(unique(collist_bw_pair())[1],"_high")
      write.table(df, file, row.names = F, sep = "\t", quote = F)
    }
  )
  output$download_input_peak_distribution = downloadHandler(
    filename = function() {
      paste0(unique(collist_bw_pair())[2],"_high_peak_distribution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(updistribution())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_deg_peak_distribution = downloadHandler(
    filename = function() {
      paste0(unique(collist_bw_pair())[1],"_high_peak_distribution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(downdistribution())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ###Peak pattern comparison up--------
  output$peak_pattern_up_additional <- renderUI({
    if(!is.null(bws())){
      if(length(list.files("./Volume/")) > 0){
        list <- list.files("Volume",full.names = T,recursive=T)
        bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
        selectInput("peak_pattern_up_add",label="Select additional bigwig files",choices = bw,multiple=T)
      }else{
        fileInput("peak_pattern_up_add",
                  "Select additional bigwig files",
                  accept = c("bw","BigWig"),
                  multiple = TRUE,
                  width = "80%")
      }
    }
  })
  pre_up_additional <-reactive({
    if(!is.null(input$peak_pattern_up_add)){
      if(length(list.files("./Volume/")) > 0){
        files <- input$peak_pattern_up_add
        name <- gsub(".+\\/","",input$peak_pattern_up_add)
        names(files) <- gsub("\\..+$", "", name)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$peak_pattern_up_add[, 1])){
          file <- input$peak_pattern_up_add[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$peak_pattern_up_add[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
      }
      return(bigwig_breakline(files))
    }
  })
  up_additional <-reactive({
    return(bws_ordering(bws=pre_up_additional(),sample_order=input$sample_order_pair_pattern))
  })
  output$peak_pattern_up_heat_range <- renderUI({
    if(!is.null(deg_result())){
      withProgress(message = "Preparing peak pattern for up peak",{
        rg <- pair_pattern_range_up()
        sliderInput("peak_up_range","Intensity range",value=rg,min = 0,max=ceiling(rg*2),step=ceiling(rg*2)/100)
      })
    }
  })
  pair_pattern_range_up <- reactive({
    rg <- c()
    if(input$data_file_type != "Row2"){
      if(is.null(bws())) validate("Bigwig files are required.")
    }
    sig <- peak_pattern_function(grange=peak_up_grange(), files=bws_order(),
                                 additional=up_additional(),plot = FALSE)
    for(name in names(sig)){
      rg <- c(rg, mean(sig[[name]][,50]))
    }
    rg <- max(rg) + max(rg)*0.1
    return(rg)
  })
  pair_pattern_range_down <- reactive({
    rg <- c()
    if(input$data_file_type != "Row2"){
      if(is.null(bws())) validate("Bigwig files are required.")
    }
    sig <- peak_pattern_function(grange=peak_down_grange(), files=bws_order(),
                                 additional=up_additional(),plot = FALSE)
    for(name in names(sig)){
      rg <- c(rg, mean(sig[[name]][,50]))
    }
    rg <- max(rg) + max(rg)*0.1
    return(rg)
  })
  
  peak_up_grange <- reactive({
    if(!is.null(data_degcount_down())){
      if(input$Genomic_region == "Genome-wide"){
        up <- range_changer(data_degcount_up())
        up <- with(up, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
      }else{
        data <- deg_result() %>% dplyr::filter(log2FoldChange < -log2(input$fc) & padj < input$fdr)
        up <- symbol2gene_id(data,org1()) %>% distinct(gene_id, .keep_all = T)
        up <- subset(promoter_region(), gene_id %in% up$gene_id) 
      }
      return(up)
    }else return(NULL)
  })
  
  peak_up_alinedHeatmap <- reactive({
    if(!is.null(input$peak_up_range) && input$peak_up_range > 0){
      bigwig <- bigwig_breakline(bws_order())
      heatmap <- peak_pattern_function(grange=peak_up_grange(), files=bigwig,
                                       additional=up_additional(),rg = input$peak_up_range)
      return(heatmap)
    }
  })
  output$peak_pattern_up_heatmap <- renderPlot({
    withProgress(message = "feature aligned heatmap",{
      if(!is.null(deg_result())){
        plot_grid(peak_up_alinedHeatmap()[["heat"]])
      }
    })
  })
  output$peak_pattern_up_line <- renderPlot({
    withProgress(message = "feature aligned distribution",{
      if(!is.null(deg_result()) && !is.null(input$peak_up_range)){
        matplot(peak_up_alinedHeatmap()[["line"]],
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_up_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
      }
    })
  })
  
  output$download_peak_pattern_up_heatmap = downloadHandler(
    filename = function(){
      paste0(unique(collist_bw_pair())[2],"_high_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_up_alinedHeatmap()[["heat"]]))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_up_line = downloadHandler(
    filename = function(){
      paste0(unique(collist_bw_pair())[2],"_high_lineplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        matplot(peak_up_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_up_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
        dev.off()
        incProgress(1)
      })
    }
  )
  ###Peak pattern comparison down--------
  output$peak_pattern_down_heat_range <- renderUI({
    if(!is.null(deg_result())){
      withProgress(message = "Preparing peak pattern for down peak",{
        rg <- pair_pattern_range_down()
        sliderInput("peak_down_range","Intensity range",value=rg,min = 0,max=ceiling(rg*2),step=ceiling(rg*2)/100)
      })
    }
  })
  peak_down_grange <- reactive({
    if(!is.null(data_degcount_down())){
      if(input$Genomic_region == "Genome-wide"){
        down <- range_changer(data_degcount_down())
        down <- with(down, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
      }else{
        data <- deg_result() %>% dplyr::filter(log2FoldChange > log2(input$fc) & padj < input$fdr)
        down <- symbol2gene_id(data,org1()) %>% distinct(gene_id, .keep_all = T)
        down <- subset(promoter_region(), gene_id %in% down$gene_id)
      }
      return(down)
    }else return(NULL)
  })
  
  peak_down_alinedHeatmap <- reactive({
    if(!is.null(input$peak_down_range) && input$peak_down_range > 0){
      bigwig <- bigwig_breakline(bws_order())
      heatmap <- peak_pattern_function(grange=peak_down_grange(), files=bigwig,
                                       additional=up_additional(),rg = input$peak_down_range)
      return(heatmap)
    }
  })
  output$peak_pattern_down_heatmap <- renderPlot({
    withProgress(message = "feature aligned heatmap",{
      if(!is.null(deg_result())){
        plot_grid(peak_down_alinedHeatmap()[["heat"]])
      }
    })
  })
  output$peak_pattern_down_line <- renderPlot({
    withProgress(message = "feature aligned distribution",{
      if(!is.null(deg_result()) && !is.null(input$peak_down_range)){
        matplot(peak_down_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_down_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
      }
    })
  })
  output$download_peak_pattern_down_heatmap = downloadHandler(
    filename = function(){
      paste0(unique(collist_bw_pair())[1],"_high_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_down_alinedHeatmap()[["heat"]]))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_down_line = downloadHandler(
    filename = function(){
      paste0(unique(collist_bw_pair())[1],"_high_lineplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        matplot(peak_down_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_down_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##Regulatory potential------
  output$pairRNAseqresult <- renderUI({
    if(input$RNAseq_data_type == "Result") {
      label <- "Select RNA-seq DEG result file"
    }else{
      label <- "Select RNA-seq raw count file"
    }
    fileInput("pair_DEG_result",
              label,
              accept = c("txt","csv","xlsx"),
              multiple = TRUE,
              width = "80%")
  })
  pre_RNAseq_file <- reactive({
    withProgress(message = "Importing RNA-seq data, please wait",{
      tmp <- input$pair_DEG_result$datapath 
      if(is.null(input$pair_DEG_result) && input$goButton > 0 )  {
        if(input$RNAseq_data_type == "Result") tmp = "data/RNAseq.txt" else tmp = "data/RNAseq_count.txt"
      }
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") {
          df2 <- readxl::read_xlsx(tmp) 
          df2 <- as.data.frame(df2)
          df <- try(data.frame(row.names = df2[,1]),silent = T)
          if(class(df) != "try-error") {
            if(dim(df2)[2] == 2){
              df <- data.frame(row.names = df2[,1],a = df2[,2])
              colnames(df)[1] <- colnames(df2)[2]
            }else{
              rownames(df2) <- df2[,1]
              df <- df2[,-1]
              colnames(df) <- gsub("-",".",colnames(df))
            }
          }
        }
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
        rownames(df) = gsub("\"", "", rownames(df))
        if(length(colnames(df)) != 0){
          if(str_detect(colnames(df)[1], "^X\\.")){
            colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
          }
        }
        return(df)
      }
    })
  })
  output$pair_RNAseq_order <- renderUI({
    if(input$RNAseq_data_type != "Result"){
      selectInput("pair_RNAseq_order","Select samples:",
                  choices = colnames(pre_RNAseq_file()),selected = colnames(pre_RNAseq_file()),multiple = T)
    }
  })
  RNAseq_file <- reactive({
    if(input$RNAseq_data_type != "Result"){
      count <- pre_RNAseq_file()
      order <- input$pair_RNAseq_order
      data <- try(count[,order])
      if(length(data) == 1){
        if(class(data) == "try-error") validate("")
      }
      return(data)
    }
  })
  updateCounter_DEGanalysis <- reactiveValues(i = 0)
  
  observe({
    input$DEGanalysis_Button
    isolate({
      updateCounter_DEGanalysis$i <- updateCounter_DEGanalysis$i + 1
    })
  })
  
  #Restart
  observeEvent(RNAseq_file(), {
    isolate(updateCounter_DEGanalysis$i == 0)
    updateCounter_DEGanalysis <<- reactiveValues(i = 0)
  }) 
  pre_RNAseqDEG <- reactive({
    if(input$DEGanalysis_Button > 0 && updateCounter_DEGanalysis$i > 0){
      count <- RNAseq_file()
      if(!is.null(count)){
        collist <- gsub("\\_.+$", "", colnames(count))
        if(length(unique(collist)) == 2){
          group <- data.frame(con = factor(collist))
          dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
          dds$con <- factor(dds$con, levels = unique(collist))
          dds <- DESeq(dds)
          return(dds)
        }else validate(print(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 2 conditions.")))
      }
    }
  })
  RNAseqDEG <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if(input$RNAseq_data_type == "Result"){
        return(pre_RNAseq_file())
      }else{
        if(input$DEGanalysis_Button > 0 && updateCounter_DEGanalysis$i > 0){
          count <- RNAseq_file()
          collist <- gsub("\\_.+$", "", colnames(count))
          dds <- pre_RNAseqDEG()
          contrast <- c("con", unique(collist))
          res <- results(dds,  contrast = contrast)
          res <- as.data.frame(res)
          return(res)
        }
      }
    })
  })
  RNAseqDEG_norm <- reactive({
    if(input$RNAseq_data_type != "Result" && !is.null(pre_RNAseqDEG())){
      count <- RNAseq_file()
      collist <- gsub("\\_.+$", "", colnames(count))
      dds <- pre_RNAseqDEG()
      contrast <- c("con", unique(collist))
      normalized_counts <- counts(dds, normalized=TRUE)
      print(head(normalized_counts))
      return(normalized_counts)
    }
  })
  gene_type_pair_DEG_result <- reactive({
    return(gene_type(my.symbols=rownames(RNAseqDEG()),org=org1(),Species=input$Species))
  })
  output$pair_RNAseq_raw <- renderDataTable({
    if(input$RNAseq_data_type != "Result"){
      if(!is.null(RNAseq_file())){
        RNAseq_file() 
      }
    }
  })
  output$RNAseq_condition <- renderText({
    if(!is.null(RNAseqDEG())){
      paste0(input$RNAseq_cond1_pair, " vs ", input$RNAseq_cond2_pair, "\n",
             "Please confirm if the Log2FoldChange in the uploaded result file was calulated as Log2(",input$RNAseq_cond1_pair,"/",input$RNAseq_cond2_pair,").")
    }
  })
  output$pair_DEG_result <- renderDataTable({
    if(!is.null(RNAseqDEG())){
      RNAseqDEG() 
    }
  })
  RNAseq_name <- reactive({
    if(input$RNAseq_data_type != "Result"){
      if(!is.null(RNAseq_file())){
        collist <- gsub("\\_.+$", "", colnames(RNAseq_file()))
        cond1 <- paste0(unique(collist)[1],"_high")
        cond2 <- paste0(unique(collist)[2],"_high")
      }
    }else{
      cond1 <- paste0(input$RNAseq_cond1_pair,"_high")
      cond2 <- paste0(input$RNAseq_cond2_pair,"_high")
    }
    return(c(cond1,cond2))
  })
  
  RNAseqDEG_anno <- reactive({
    return(RNAseqDEG_ann(RNAdata=RNAseqDEG(),Species=input$Species,gene_type = gene_type_pair_DEG_result(),input_type=input$RNAseq_data_type))
  })
  
  output$RNAseq_RPmean <- renderText({
    if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
      paste0("Epigenome analysis: ",unique(collist_bw_pair())[1], " vs ", unique(collist_bw_pair())[2], "\n",
             "RP > 0 genes have a higher epigenetic potential in ", unique(collist_bw_pair())[2]," condition.\n",
             "RP < 0 genes have a higher epigenetic potential in ", unique(collist_bw_pair())[1]," condition.\n",
             "RP = 0 genes have no significant difference in epigenetic potential between ", unique(collist_bw_pair())[1], " and ", unique(collist_bw_pair())[2]," conditions.")
    }else{
      paste0("RP > 0 gene is associated with the indicated differential peak within ",input$peak_distance," kb from its TSS.\n",
             "RP = 0 gene is not associated with the indicated differential peak within ",input$peak_distance," kb from its TSS.") 
    }
  })
  
  mmAnno_pair <- reactive({
    mmAnno_list <- list()
    mmAnno_list[[paste0(unique(collist_bw_pair())[2],"_high")]] <- mmAnno(peak=peak_up_grange(),
                                                                          genomic_region=input$Genomic_region,
                                                                          txdb=txdb(),
                                                                          peak_distance=input$peak_distance,
                                                                          mode=input$RNAseq_mode,
                                                                          group_name=paste0(unique(collist_bw_pair())[2],"_high"),
                                                                          distribution=updistribution()$peaks,
                                                                          DAR=DAR_for_withRNAseq())
    mmAnno_list[[paste0(unique(collist_bw_pair())[1],"_high")]] <- mmAnno(peak=peak_down_grange(),
                                                                          genomic_region=input$Genomic_region,
                                                                          txdb=txdb(),
                                                                          peak_distance=input$peak_distance,
                                                                          mode=input$RNAseq_mode,
                                                                          group_name=paste0(unique(collist_bw_pair())[1],"_high"),
                                                                          distribution=downdistribution()$peaks,
                                                                          DAR=DAR_for_withRNAseq())
    return(mmAnno_list)
  })
  
  DAR_for_withRNAseq <- reactive({
    return(DAR_withRNAseq(DAR=pre_DEG_result(),genomic_region=input$Genomic_region,Species=input$Species))
  })
  
  RP <- reactive({
    withProgress(message = paste0("Calculating regulatory potential (RP_definition: ",input$RNAseq_mode_option,")"),{
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        mmAnno_up <- mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]]
        if(!is.null(mmAnno_up)) {
          result_geneRP_up <- modified_calcRP_TFHit(mmAnno = mmAnno_up,Txdb = txdb(),mode=input$RNAseq_mode_option)
          print(head(result_geneRP_up))
          colnames(result_geneRP_up)[2:4] <- paste0(colnames(result_geneRP_up)[2:4], "_up")
        }else result_geneRP_up <- NULL
        mmAnno_down <- mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]]
        if(!is.null(mmAnno_down)) {
          result_geneRP_down <- modified_calcRP_TFHit(mmAnno = mmAnno_down,Txdb = txdb(),mode=input$RNAseq_mode_option)
          colnames(result_geneRP_down)[2:4] <- paste0(colnames(result_geneRP_down)[2:4], "_down")
          print(head(result_geneRP_down))
        }else result_geneRP_down <- NULL
        if(!is.null(result_geneRP_up) && !is.null(result_geneRP_down)){
          result_geneRP <- merge(result_geneRP_up,result_geneRP_down,by="gene_id",all =T) 
          result_geneRP$sumRP <- apply(data.frame(result_geneRP$sumRP_up,
                                                  -result_geneRP$sumRP_down),1,sum,na.rm=TRUE)
        }
        if(!is.null(result_geneRP_up) && is.null(result_geneRP_down)) {
          result_geneRP <-result_geneRP_up
          result_geneRP$sumRP <- result_geneRP$sumRP_up
        }
        if(is.null(result_geneRP_up) && !is.null(result_geneRP_down)) {
          result_geneRP <-result_geneRP_down
          result_geneRP$sumRP <- -result_geneRP$sumRP_down
        }
        if(is.null(result_geneRP_up) && is.null(result_geneRP_down)){
          return(NULL)
        }else{
          result_geneRP <- result_geneRP %>% dplyr::arrange(-sumRP)
          result_geneRP$RP_rank <- rownames(result_geneRP) %>% as.numeric()
          return(result_geneRP)
        }
      }else{
        if(!is.null(mmAnno_pair())){
          RP_list <- list()
          for(name in names(mmAnno_pair())){
            data <- RP_f(mmAnno=mmAnno_pair()[[name]],txdb=txdb(),mode=input$RNAseq_mode_option)
            RP_list[[name]] <- data
          }
          return(RP_list)
        }
      }
      incProgress(1)
    })
  })
  regulatory_potential <- reactive({
    if(input$Species != "not selected" && !is.null(RP()) && !is.null(RNAseqDEG_anno())){
      data <- RNAseqDEG_anno()
      result_geneRP <- RP()
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        merge_data <- integrate_ChIP_RNA(
          result_geneRP = result_geneRP,
          result_geneDiff = data,lfc_threshold = log(input$DEG_fc,2),
          padj_threshold = input$DEG_fdr, name=RNAseq_name())
        return(merge_data)
      }else{
        RP_list <- list()
        for(name in names(result_geneRP)){
          data <- regulatory_potential_f(species=input$Species,data=RNAseqDEG_anno(),
                                         result_geneRP= result_geneRP[[name]],DEG_fc=input$DEG_fc,
                                         DEG_fdr=input$DEG_fdr,name=RNAseq_name())
          RP_list[[name]] <- data
        }
        return(RP_list)
      }
    }
  })
  
  output$ks_plot <- renderPlot({    
    if(!is.null(RNAseqDEG()) && 
       !is.null(regulatory_potential()) && !is.na(input$DEG_fdr) && 
       !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        regulatory_potential()
      }else{
        down <- regulatory_potential()[[RNAseq_name()[1]]] + ggtitle(RNAseq_name()[1])
        up <- regulatory_potential()[[RNAseq_name()[2]]] + ggtitle(RNAseq_name()[2])
        plot_grid(down,up,nrow = 1)
      }
    }
  })
  ks_tables <- reactive({
    if(!is.null(RNAseqDEG_anno()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && 
       input$Species != "not selected"  && !is.null(mmAnno_pair())){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        df <- regulatory_potential()$statistics
        print(df)
      }else{
        df <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
        for(name in names(regulatory_potential())){
          df2 <- regulatory_potential()[[name]]$statistics
          df2$Peak <- name
          df <- rbind(df,df2)
        }
        df <- df %>% dplyr::select(Peak,everything())
      }
      return(df) 
    }
  })
  output$ks_plot_table <- DT::renderDataTable({
    if(!is.null(RNAseqDEG()) && 
       !is.null(regulatory_potential()) && !is.na(input$DEG_fdr) && 
       !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      ks_tables()
    }
  })
  output$DEG_fc <- renderUI({
    numericInput("DEG_fc","Fold change cutoff for RNA-seq data",
                 min=0,max=NA,value=1.5,step = 0.5)
  })
  output$DEG_fdr <- renderUI({
    numericInput("DEG_fdr","FDR cutoff for RNA-seq data",
                 min=0,max=1, value=0.05,step = 0.001)
  })
  output$RNAseqGroup <- renderUI({
    if(input$Species != "not selected" &&!is.null(mmAnno_pair())){
      if(!is.null(RP_all_table())){
        group <- unique(pre_RP_selected_table()$Group)
        selectInput("RNAseqGroup","Group (Epigenome:RNAseq)",
                    group,
                    multiple = FALSE)
      }
    }
  })
  pari_RNAseq_boxplot <- reactive({
    RNA <- RNAseqDEG_anno()
    if(is.null(RP()) || is.null(RNA)) validate("")
    withProgress(message = "Preparing boxplot",{
      cond1 <- gsub("\\_.+$", "", RNAseq_name()[1])
      cond2 <- gsub("\\_.+$", "", RNAseq_name()[2])
      RNA <- dplyr::filter(RNA, !is.na(gene_id))
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        data <- merge(RNA,RP(), by="gene_id",all=T)
        print(dim(data))
        data <- dplyr::filter(data, baseMean != 0)
        print(head(data))
        data$group <- "Others"
        data$group[data$sumRP > 1] <- "1 < RP"
        data$group[data$sumRP <= 1 & data$sumRP > 0.1] <- "0.1 < RP < 1"
        data$group[data$sumRP <= 0.1 & data$sumRP > 0.01] <- "0.01 < RP < 0.1"
        data$group[data$sumRP <= 0.01 & data$sumRP > 0] <- "0 < RP < 0.01"
        data$group[data$sumRP >= -0.01 & data$sumRP < 0] <- "-0.01 < RP < 0"
        data$group[data$sumRP >= -0.1 & data$sumRP < -0.01] <- "-0.1 < RP < -0.01"
        data$group[data$sumRP >= -1 & data$sumRP < -0.1] <- "-1 < RP < -0.1"
        data$group[data$sumRP < -1] <- "RP < -1"
        level <- NULL
        col <- NULL
        if(dim(dplyr::filter(data, group == "RP < -1"))[1] != 0) {
          level <- c(level,"RP < -1")
          col <- c(col,"blue")
        }
        if(dim(dplyr::filter(data, group == "-1 < RP < -0.1"))[1] != 0){
          level <- c(level,"-1 < RP < -0.1")
          col <- c(col,"#00BFC4")
        }
        if(dim(dplyr::filter(data, group == "-0.1 < RP < -0.01"))[1] != 0) {
          level <- c(level,"-0.1 < RP < -0.01")
          col <- c(col,"lightgreen")
        }
        if(dim(dplyr::filter(data, group == "-0.01 < RP < 0"))[1] != 0) {
          level <- c(level,"-0.01 < RP < 0")
          col <- c(col,"darkseagreen")
        }
        if(dim(dplyr::filter(data, group == "Others"))[1] != 0) {
          level <- c(level,"Others")
          col <- c(col,"gray")
        }
        if(dim(dplyr::filter(data, group == "0 < RP < 0.01"))[1] != 0) {
          level <- c(level,"0 < RP < 0.01")
          col <- c(col,"orange3")
        }
        if(dim(dplyr::filter(data, group == "0.01 < RP < 0.1"))[1] != 0) {
          level <- c(level,"0.01 < RP < 0.1")
          col <- c(col,"orange")
        }
        if(dim(dplyr::filter(data, group == "0.1 < RP < 1"))[1] != 0) {
          level <- c(level,"0.1 < RP < 1")
          col <- c(col,"#F8766D")
        }
        if(dim(dplyr::filter(data, group == "1 < RP"))[1] != 0) {
          level <- c(level,"1 < RP")
          col <- c(col,"red")
        }
        data$group <- factor(data$group,levels=level,ordered=TRUE)
        data$log10FoldChange <- log10(2^(data$log2FoldChange))
        collist <- unique(data$group)
        if (length(collist) >= 3){
          stat.test <- data %>% tukey_hsd(log10FoldChange ~ group)
          stat.test <- stat.test %>% add_significance("p.adj")
          stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
        }else{
          group1 <- dplyr::filter(data, group == collist[1])
          group2 <- dplyr::filter(data, group == collist[2])
          if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
            stat.test <- data %>% t_test(log10FoldChange ~ group)
            stat.test <- stat.test %>% add_significance()
            stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
          }else stat.test <- NULL
        }
        p <- try(ggpubr::ggboxplot(data, x = "group", y = "log10FoldChange",
                                   fill = "group", scales = "free", 
                                   xlab = FALSE, ylab = "log10FoldChange")+theme_bw(base_size = 15)+guides(fill=guide_legend(title="RP\nstatus"))+
                   xlab(NULL)+ylab(paste0("RNAseq log10(",cond2,"/",cond1,")"))+scale_fill_manual(values = col) + scale_x_discrete(labels = label_wrap_gen(8)))
        if(input$Statistics !="not selected") p <- p + stat_pvalue_manual(stat.test,hide.ns = T, size = 5)  
        stat.test <- stat.test[,2:9]
        genelist <- data.frame(Symbol = data$Symbol, group = data$group,RNAlog2FC = -1*data$log2FoldChange,
                               sumRP = data$sumRP, gene_id = data$gene_id)
        bed_list <- list()
        for(name in unique(genelist$group)){
          if(name != "Others"){
            table <- genelist %>% dplyr::filter(group == name)
            gene <- table$gene_id
            up_peak3 <- NULL
            if(sum(is.element(c("0 < RP < 0.01","0.01 < RP < 0.1","0.1 < RP < 1","1 < RP"),name)) == 1) peak_type <- paste0(unique(collist_bw_pair())[2],"_high") else peak_type <- paste0(unique(collist_bw_pair())[1],"_high")
            up_peak <- subset(mmAnno_pair()[[peak_type]], gene_id %in% gene)
            up_peak2 <- as.data.frame(up_peak)
            up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
            up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
            up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
            mcols(up_peak3) <- DataFrame(Group = name)
            bed_list[[name]] <- up_peak3
          }
        }
        
      }else{
        
        df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        genelist <- list()
        bed_list <- list()
        for(name in names(RP())){
          data <- merge(RNA,RP()[[name]], by="gene_id",all=T)
          data <- dplyr::filter(data, baseMean != 0)
          data$group <- "Others"
          data$group[data$sumRP > 2] <- "2 < RP"
          data$group[data$sumRP <= 2 & data$sumRP > 1] <- "1 < RP < 2"
          data$group[data$sumRP <= 1 & data$sumRP > 0.1] <- "0.1 < RP < 1"
          data$group[data$sumRP <= 0.1 & data$sumRP > 0.01] <- "0.01 < RP < 0.1"
          data$group[data$sumRP <= 0.01 & data$sumRP > 0] <- "0 < RP < 0.01"
          
          level <- NULL
          col <- NULL
          if(dim(dplyr::filter(data, group == "Others"))[1] != 0) {
            level <- c(level,"Others")
            col <- c(col,"gray")
          }
          if(dim(dplyr::filter(data, group == "0 < RP < 0.01"))[1] != 0) {
            level <- c(level,"0 < RP < 0.01")
            col <- c(col,"blue")
          }
          if(dim(dplyr::filter(data, group == "0.01 < RP < 0.1"))[1] != 0) {
            level <- c(level,"0.01 < RP < 0.1")
            col <- c(col,"#00BFC4")
          }
          if(dim(dplyr::filter(data, group == "0.1 < RP < 1"))[1] != 0) {
            level <- c(level,"0.1 < RP < 1")
            col <- c(col,"lightgreen")
          }
          if(dim(dplyr::filter(data, group == "1 < RP < 2"))[1] != 0) {
            level <- c(level,"1 < RP < 2")
            col <- c(col,"#F8766D")
          }
          if(dim(dplyr::filter(data, group == "2 < RP"))[1] != 0) {
            level <- c(level,"2 < RP")
            col <- c(col,"red")
          }
          data$group <- factor(data$group,levels=level,ordered=TRUE)
          data$intersection <- as.factor(name)
          df <- rbind(df,data)
          genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,RNAlog2FC = -1*data$log2FoldChange,
                                         sumRP = data$sumRP, gene_id = data$gene_id)
          bed_list1 <- list()
          for(name2 in unique(data$group)){
            if(name2 != "Others"){
              table <- genelist[[name]] %>% dplyr::filter(group == name2)
              gene <- table$gene_id
              up_peak3 <- NULL
              if(length(gene) != 0){
                up_peak <- subset(mmAnno_pair()[[name]], gene_id %in% gene)
                up_peak2 <- as.data.frame(up_peak)
                up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
                up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
                up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
                mcols(up_peak3) <- DataFrame(Group = paste0(name,":",name2))
                bed_list1[[name2]] <- up_peak3
              }
            }
          }
          bed_list[[name]] <- bed_list1
        }
        data <- df 
        data$log10FoldChange <- log10(2^data$log2FoldChange)
        data$log10FoldChange <- as.numeric(data$log10FoldChange)
        for(name in names(RP())){
          check <- data %>% dplyr::filter(intersection == name) %>% 
            dplyr::filter(group != "Others") %>% summarise(n())
          print(check)
          if(check <= 1) data <- data %>% dplyr::filter(intersection != name)
        }
        if(dim(data)[1] == 0) validate("boxplot: There are few genes with |RP| > 1")
        
        data$intersection <- gsub("-","-\n",data$intersection)
        collist <- unique(data$group)
        if (length(collist) >= 3){
          stat.test <- data %>% dplyr::group_by(intersection) %>% 
            tukey_hsd(log10FoldChange ~ group)
          stat.test <- stat.test %>% add_significance("p.adj")
          stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
        }else{
          group1 <- dplyr::filter(data, group == collist[1])
          group2 <- dplyr::filter(data, group == collist[2])
          if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
            stat.test <- data %>% dplyr::group_by(intersection) %>% 
              t_test(log10FoldChange ~ group)
            stat.test <- stat.test %>% add_significance()
            stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
          }else stat.test <- NULL
        }
        p <- try(ggpubr::ggboxplot(data, x = "group", y = "log10FoldChange",
                                   fill = "group", scales = "free", 
                                   xlab = FALSE, ylab = paste0("RNAseq log10(",cond2,"/",cond1,")"))+theme_bw(base_size = 15)+
                   xlab(NULL)+scale_fill_manual(values = col) + scale_x_discrete(labels = label_wrap_gen(8)))
        if(input$Statistics !="not selected") p <- p + stat_pvalue_manual(stat.test,hide.ns = T, size = 5)         
        p <- facet(p, facet.by = "intersection",
                   panel.labs.background = list(fill = "transparent", color = "transparent"),
                   scales = "free", short.panel.labs = T, panel.labs.font = list(size=15))
        stat.test <- stat.test[,-2]
        stat.test <- stat.test[,1:9]
        colnames(stat.test)[1] <- "Peak"
      }
      df <- list()
      df[["plot"]] <- p
      df[["statistical_test"]] <- stat.test
      df[["genelist"]] <- genelist
      df[["bedlist"]] <- bed_list
      print(bed_list)
      incProgress(1)
    })
    return(df)
  })
  output$RNAseq_boxplot_error <- renderText({
    if(!is.null(RNAseqDEG()) && !is.null(RP()) &&
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        RNA <- RNAseqDEG_anno()
        RNA <- dplyr::filter(RNA, !is.na(gene_id))
        data <- merge(RNA,RP(), by="gene_id",all=T)
        data$group <- "Others"
        data$group[data$sumRP > input$DEG_fdr] <- "1 < RP"
        data$group[data$sumRP < -input$DEG_fdr] <- "RP < -1"
        data$group <- factor(data$group,levels=c("Others","1 < RP","RP < -1"),ordered=TRUE)
        collist <- unique(data$group)
        if (length(collist) < 3){
          group1 <- dplyr::filter(data, group == collist[1])
          group2 <- dplyr::filter(data, group == collist[2])
          if(length(rownames(group1)) <= 1 || length(rownames(group2)) <= 1){
            print("boxplot: There are few genes with |RP| > 1")
          }
        }
      }
    }
  })
  RNAseq_popu <- reactive({
    if(!is.null(collist_bw_pair()) && !is.null(RP_all_table())){
      withProgress(message = "Preparing bar plot",{
        up_name <- RNAseq_name()[2]
        down_name <- RNAseq_name()[1]
        if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
          epi_up_name <- paste0(unique(collist_bw_pair())[2], "_high_only")
          epi_down_name <- paste0(unique(collist_bw_pair())[1], "_high_only")
          table <- RP_all_table() %>% dplyr::mutate(
            type =if_else(withUpPeakN > 0 & withDownPeakN > 0, "both",
                          if_else(withUpPeakN > 0 & withDownPeakN == 0, epi_up_name,
                                  if_else(withUpPeakN == 0 & withDownPeakN > 0, epi_down_name,"no")))
          )
          table$Group <- gsub(".+\\:","",table$Group)
          table$Group <- paste0(table$Group," genes")
          table$type <- factor(table$type,levels = c(epi_down_name,epi_up_name,"both"))
          p2 <- ggplot(table,aes(x = TotalPeakN, fill = type)) + geom_bar(position = "stack") +
            theme_bw(base_size = 15)+facet_wrap(~Group,scales = "free") +
            scale_fill_manual(values = c("#00BFC4","#F8766D","grey")) +
            xlab("Number of associated peaks")+guides(fill=guide_legend(title="associated_peak_type"))
        }else{
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          for(name in names(RP_all_table())){
            epi_name <- paste0(name, "_associated")
            table <- RP_all_table()[[name]] %>% dplyr::mutate(
              type =if_else(withPeakN > 0, epi_name, "not_associated")
            )
            table$Group <- gsub(".+\\:","",table$Group)
            table$Group <- paste0(table$Group," genes")
            table$type <- factor(table$type,levels = c(epi_name,"not_associated"))
            table$intersection <- name
            df <- rbind(df, table)
          }
          table <- df
          table$intersection <- gsub("-","-\n",table$intersection)
          table2 <- table %>% group_by(Group, intersection, withPeakN) %>%
            summarise(count = n(), .groups = "drop")
          p2 <- ggplot(table2,aes(x = withPeakN,y= count,fill = intersection)) +
            geom_col(position=position_dodge2(preserve = "single")) +
            theme_bw(base_size = 15)+facet_wrap(~Group,scales = "free",ncol = 3) +
            xlab("Number of associated peaks")+guides(fill=guide_legend(title="associated_\npeak_type"))
        }
        incProgress(1)
      })
      return(p2)
    }
  })
  output$int_boxplot <- renderPlot({
    if(!is.null(RNAseqDEG()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      pari_RNAseq_boxplot()[["plot"]]
    }
  })
  output$int_boxplot_table <- DT::renderDataTable({
    if(!is.null(RNAseqDEG()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      if(input$Statistics !="not selected"){
        pari_RNAseq_boxplot()[["statistical_test"]]
      }
    }
  })
  
  output$int_bar <- renderPlot({
    if(!is.null(RNAseqDEG()) && !is.null(RNAseq_popu()) && !is.null(ChIPseq_popu()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      gridExtra::grid.arrange(RNAseq_popu(), ChIPseq_popu(), ncol = 1)
    }
  })
  ChIPseq_popu <- reactive({
    RNA <- RNAseqDEG_anno()
    mmano_up <- as.data.frame(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]])
    mmano_down <- as.data.frame(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])
    if(!is.null(RNA) && !is.null(mmano_up) && !is.null(mmano_down) && !is.null(RNAseq_name())){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        mmano <- rbind(mmano_up,mmano_down)
        mmano$locus <- paste0(mmano$seqnames,":",mmano$start,"-",mmano$end)
        up_name <- paste0(RNAseq_name()[2]," gene only")
        down_name <- paste0(RNAseq_name()[1]," gene only")
        merge <- merge(mmano,RNA, by = "gene_id")
        merge <- merge  %>% dplyr::mutate(Up_gene = if_else(padj < input$DEG_fdr & log2FoldChange > log(input$DEG_fc,2), 1, 0),
                                          Down_gene = if_else(padj < input$DEG_fdr & log2FoldChange < -log(input$DEG_fc,2), 1, 0),
                                          NS_gene = if_else(padj >= input$DEG_fdr | baseMean == 0 | is.na(padj) |
                                                              (padj <= input$DEG_fdr & abs(log2FoldChange) <= log(input$DEG_fc,2)), 1, 0))
        merge[is.na(merge)] <- 0
        merge2 <- merge %>% group_by(locus,Group) %>% summarise(Total_associated_gene = sum(Up_gene)+sum(Down_gene)+sum(NS_gene),Group = Group,
                                                                Up_gene = sum(Up_gene), Down_gene = sum(Down_gene), NS_gene = sum(NS_gene))
        table <- merge2 %>% dplyr::mutate(type = if_else(Up_gene > 0 & Down_gene == 0 & NS_gene == 0, up_name,
                                                         if_else(Up_gene == 0 & Down_gene > 0 & NS_gene == 0, down_name,
                                                                 if_else(Up_gene == 0 & Down_gene == 0 & NS_gene > 0, "NS gene only", "Multiple type of genes"))))
        table$type <- factor(table$type,levels = c(down_name,up_name,"NS gene only","Multiple type of genes"))
        table$Group <- paste0(table$Group," peak")
        if(max(table$Total_associated_gene) < 5) table$Total_associated_gene <- as.character(table$Total_associated_gene) else table$Total_associated_gene <- as.numeric(table$Total_associated_gene)
        col <- c("#00BFC4","#F8766D","grey","black")
        p2 <- ggplot(table,aes(x = Total_associated_gene, fill = type)) + geom_bar(position = "stack") +
          theme_bw(base_size = 15)+facet_wrap(~Group,scales = "free") +
          scale_fill_manual(values = col) + 
          ylab("Number of peaks")+
          xlab("Number of associated genes")+guides(fill=guide_legend(title="associated_gene_type"))
      }else{
        df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        for(name in names(mmAnno_pair())){
          mmano <- as.data.frame(mmAnno_pair()[[name]])
          mmano$locus <- paste0(mmano$seqnames,":",mmano$start,"-",mmano$end)
          merge <- merge(mmano,RNA, by = "gene_id")
          up_name <- paste0(RNAseq_name()[2]," gene only")
          down_name <- paste0(RNAseq_name()[1]," gene only")
          merge <- merge  %>% dplyr::mutate(Up_gene = if_else(padj < input$DEG_fdr & log2FoldChange > log(input$DEG_fc,2), 1, 0),
                                            Down_gene = if_else(padj < input$DEG_fdr & log2FoldChange < -log(input$DEG_fc,2), 1, 0),
                                            NS_gene = if_else(padj >= input$DEG_fdr | baseMean == 0 | is.na(padj) |
                                                                (padj <= input$DEG_fdr & abs(log2FoldChange) <= log(input$DEG_fc,2)), 1, 0))
          merge[is.na(merge)] <- 0
          merge2 <- merge %>% group_by(locus,Group) %>% summarise(Total_associated_gene = sum(Up_gene)+sum(Down_gene)+sum(NS_gene),Group = Group,
                                                                  Up_gene = sum(Up_gene), Down_gene = sum(Down_gene), NS_gene = sum(NS_gene))
          table <- merge2 %>% dplyr::mutate(type = if_else(Up_gene > 0 & Down_gene == 0 & NS_gene == 0, up_name,
                                                           if_else(Up_gene == 0 & Down_gene > 0 & NS_gene == 0, down_name,
                                                                   if_else(Up_gene == 0 & Down_gene == 0 & NS_gene > 0, "NS gene only", "Multiple type of genes"))))
          table$type <- factor(table$type,levels = c(down_name,up_name,"NS gene only","Multiple type of genes"))
          table$Group <- paste0(table$Group," peak")
          if(max(table$Total_associated_gene) < 5) table$Total_associated_gene <- as.character(table$Total_associated_gene) else table$Total_associated_gene <- as.numeric(table$Total_associated_gene)
          table$RNA_group <- name
          df <- rbind(df, table)
        }
        table <- df
        table$RNA_group <- gsub("-","\n",table$RNA_group)
        
        p2 <- ggplot(table,aes(x = Total_associated_gene, fill = type)) + geom_bar(position = "stack") +
          theme_bw(base_size = 15)+facet_wrap(~RNA_group,scales = "free") + ylab("Number of peaks") +
          xlab("Number of associated genes")+guides(fill=guide_legend(title="associated_gene_type"))
        col <- c("#00BFC4","#F8766D","grey","black")
        p2 <- p2 +
          scale_fill_manual(values = col) 
      }
      return(p2)
    }
  })
  RP_all_table <- reactive({
    if(!is.null(RNAseqDEG()) && !is.null(regulatory_potential()) &&
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        target_result <- regulatory_potential()$data
        epi_up_name <- paste0(unique(collist_bw_pair())[2], "_high")
        epi_down_name <- paste0(unique(collist_bw_pair())[1], "_high")
        target_result$epigenome_category <- epi_up_name
        target_result$epigenome_category[target_result$sumRP < 0] <- epi_down_name
        table <- NULL
        if(str_detect(target_result$gene_id[1], "FBgn")){
          symbol <- id_convert(my.symbols = target_result$gene_id,Species = input$Species,type = "ENSEMBL")
          symbol <- symbol %>% distinct(ENSEMBL, .keep_all = TRUE)
          symbol <- symbol$SYMBOL
        }else symbol <- target_result$Symbol
        if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]]) && !is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])) {
          table <- data.frame(Symbol = symbol,
                              Group = paste0(target_result$epigenome_category,":",target_result$gene_category),
                              RNA_log2FC = -target_result$log2FoldChange,
                              RNA_padj = target_result$padj,
                              regulatory_potential = target_result$sumRP,
                              withUpPeakN = target_result$withPeakN_up,
                              withDownPeakN = target_result$withPeakN_down,
                              gene_id = target_result$gene_id)
        }else{
          if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]]) && is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])){
            table <- data.frame(Symbol = symbol,
                                Group = paste0(target_result$epigenome_category,":",target_result$gene_category),
                                RNA_log2FC = -target_result$log2FoldChange,
                                RNA_padj = target_result$padj,
                                regulatory_potential = target_result$sumRP,
                                withUpPeakN = target_result$withPeakN_up,
                                gene_id = target_result$gene_id)
          }
          if(is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]]) && !is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])){
            table <- data.frame(Symbol = symbol,
                                Group = paste0(target_result$epigenome_category,":",target_result$gene_category),
                                RNA_log2FC = -target_result$log2FoldChange,
                                RNA_padj = target_result$padj,
                                regulatory_potential = target_result$sumRP,
                                withDownPeakN = target_result$withPeakN_down,
                                gene_id = target_result$gene_id)
          }
        }
        if(!is.null(table)) {
          table[is.na(table)] <- 0
          if(!is.null(table$withUpPeakN) && !is.null(table$withDownPeakN)) table$TotalPeakN <- table$withUpPeakN + table$withDownPeakN
          if(is.null(table$withUpPeakN) && !is.null(table$withDownPeakN)) table$TotalPeakN <- table$withDownPeakN
          if(!is.null(table$withUpPeakN) && is.null(table$withDownPeakN)) table$TotalPeakN <- table$withUpPeakN
        }
      }else{
        table_list <- list()
        for(name in names(RP())){
          target_result <- regulatory_potential()[[name]]$data
          target_result$epigenome_category <- name
          table <- NULL
          if(str_detect(target_result$gene_id[1], "FBgn")){
            symbol <- id_convert(my.symbols = target_result$gene_id,Species = input$Species,type = "ENSEMBL")
            symbol <- symbol %>% distinct(ENSEMBL, .keep_all = TRUE)
            symbol <- symbol$SYMBOL
          }else symbol <- target_result$Symbol
          if(!is.null(mmAnno_pair()[[name]])) {
            table <- data.frame(Symbol = symbol,
                                Group = paste0(target_result$epigenome_category,":",target_result$gene_category),
                                RNA_log2FC = -target_result$log2FoldChange,
                                RNA_padj = target_result$padj,
                                regulatory_potential = target_result$sumRP,
                                withPeakN = target_result$withPeakN,
                                gene_id = target_result$gene_id)
            table_list[[name]] <- table
          }
        }
        table <- table_list
      }
      return(table)
    }
  })
  
  pre_RP_selected_table <- reactive({
    if(!is.null(RNAseqDEG()) && !is.null(regulatory_potential()) &&
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        table <- RP_all_table()
      }else{
        df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
        table <- RP_all_table()
        for(name in names(table)){
          df <- rbind(df, table[[name]])
        }
        table <- df
      }
      return(table)
    }
  })
  RP_selected_table <- reactive({
    table <- pre_RP_selected_table() %>% dplyr::filter(Group == input$RNAseqGroup)
    return(table)
  })
  output$RP_table <- renderDT({
    if(!is.null(input$RNAseqGroup)){
      RP_selected_table() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  
  
  output$download_pairintbox <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"boxplot",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined") pdf_width <- 10 else pdf_width <- 12
        }else pdf_width <- input$pair_pdf_width
        fs <- c()
        dir.create("boxplot/",showWarnings = FALSE)
        boxplot_table <- paste0("boxplot/","boxplot.txt")
        boxplot <- paste0("boxplot/","boxplot.pdf")
        
        fs <- c(fs,boxplot_table,boxplot)
        pdf(boxplot, height = pdf_height, width = pdf_width)
        print(pari_RNAseq_boxplot()[["plot"]])
        dev.off()
        write.table(pari_RNAseq_boxplot()[["statistical_test"]],boxplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
        dir.create("boxplot/bed/",showWarnings = FALSE)
        if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
          boxplot_genelist <- paste0("boxplot/","boxplot_genelist.txt")
          fs <- c(fs,boxplot_genelist)
          write.table(pari_RNAseq_boxplot()[["genelist"]],boxplot_genelist,col.names = T,row.names = F,sep = "\t",quote = F)
          for(name in names(pari_RNAseq_boxplot()[["bedlist"]])){
            bed <- pari_RNAseq_boxplot()[["bedlist"]][[name]]
            bed_name <- paste0("boxplot/bed/",name,".bed")
            bed_name <- gsub(" < ","....",bed_name)
            fs <- c(fs, bed_name)
            write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
          }
        }else{
          dir.create("boxplot/gene_list/",showWarnings = FALSE)
          dir.create("boxplot/bed/",showWarnings = FALSE)
          genelist_up <- paste0("boxplot/gene_list/",RNAseq_name()[2],".txt")
          genelist_down <- paste0("boxplot/gene_list/",RNAseq_name()[1],".txt")
          up <- pari_RNAseq_boxplot()[["genelist"]][[RNAseq_name()[2]]]
          down <- pari_RNAseq_boxplot()[["genelist"]][[RNAseq_name()[1]]]
          fs <- c(fs,genelist_up,genelist_down)
          write.table(up,genelist_up,col.names = T,row.names = F,sep = "\t",quote = F)
          write.table(down,genelist_down,col.names = T,row.names = F,sep = "\t",quote = F)
          
          for(name in names(pari_RNAseq_boxplot()[["bedlist"]])){
            dir.create(paste0("boxplot/bed/",name),showWarnings = FALSE)
            for(file in names(pari_RNAseq_boxplot()[["bedlist"]][[name]])){
              bed <- pari_RNAseq_boxplot()[["bedlist"]][[name]][[file]]
              bed_name <- paste0("boxplot/bed/",name,"/",file,".bed")
              bed_name <- gsub(" < ","....",bed_name)
              fs <- c(fs, bed_name)
              write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
            }
          }
        }
        
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    },
    contentType = "application/zip"
  )
  
  output$download_pairintbar = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        gridExtra::grid.arrange(RNAseq_popu(), ChIPseq_popu(), ncol = 1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_pairKSplot = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"KSplot",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        fs <- c()
        dir.create("KSplot/",showWarnings = FALSE)
        KSplot_table <- paste0("KSplot/","KSplot.txt")
        fs <- c(fs,KSplot_table)
        write.table(ks_tables(),KSplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
        if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
          KSplot <- paste0("KSplot/","KSplot.pdf")
          fs <- c(fs,KSplot)
          pdf(KSplot, height = pdf_height, width = pdf_width)
          print(regulatory_potential())
          dev.off()
        }else{
          KSplot_up <- paste0("KSplot/",RNAseq_name()[2],".pdf")
          KSplot_down <- paste0("KSplot/",RNAseq_name()[1],".pdf")
          fs <- c(fs,KSplot_up,KSplot_down)
          down <- regulatory_potential()[[RNAseq_name()[1]]] + ggtitle(RNAseq_name()[1])
          up <- regulatory_potential()[[RNAseq_name()[2]]] + ggtitle(RNAseq_name()[2])
          pdf(KSplot_up, height = pdf_height, width = pdf_width)
          print(up)
          dev.off()
          pdf(KSplot_down, height = pdf_height, width = pdf_width)
          print(down)
          dev.off()
        }
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    },
    contentType = "application/zip"
  )
  output$download_RP_table = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"RP_table",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download",{
        fs <- c()
        dir.create("RP_table/",showWarnings = FALSE)
        dir.create("RP_table/selected_table(epigenome--RNA)/",showWarnings = FALSE)
        dir.create("RP_table/selected_bed(epigenome--RNA)/",showWarnings = FALSE)
        RP_summary <- paste0("RP_table/summary.txt")
        fs <- c(fs,RP_summary)
        write.table(pre_RP_selected_table(), RP_summary, row.names = F, sep = "\t", quote = F)
        for(name in unique(pre_RP_selected_table()$Group)){
          RP_selected <- paste0("RP_table/selected_table(epigenome--RNA)/",name,".txt")
          RP_selected_bed <- paste0("RP_table/selected_bed(epigenome--RNA)/",name,".bed")
          RP_selected <- gsub(":","--",RP_selected)
          RP_selected_bed <- gsub(":","--",RP_selected_bed)
          fs <- c(fs,RP_selected,RP_selected_bed)
          table <- pre_RP_selected_table() %>% dplyr::filter(Group == name)
          write.table(table, RP_selected, row.names = F, sep = "\t", quote = F)
          
          gene <- table$gene_id
          peak <- gsub("\\:.+$", "", name)
          y <- NULL
          if(!is.null(mmAnno_pair())) {
            up_peak <- subset(mmAnno_pair()[[peak]], gene_id %in% gene)
            up_peak2 <- as.data.frame(up_peak)
            up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
            up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
            up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
            mcols(up_peak3) <- DataFrame(Group = name)
            y <- as.data.frame(up_peak3)
          }
          write.table(y, RP_selected_bed, row.names = F, col.names = F,sep = "\t", quote = F)
        }
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    },
    contentType = "application/zip"
  )
  
  output$download_selected_RP_table = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup,"_table.txt")
    },
    content = function(file){write.table(RP_selected_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  int_selected_bed <- reactive({
    if(!is.null(RP_selected_table())){
      gene <- RP_selected_table()$gene_id
      type <- gsub("\\:.+$","", input$RNAseqGroup)
      epi_up_name <- paste0(unique(collist_bw_pair())[2], "_high")
      epi_down_name <- paste0(unique(collist_bw_pair())[1], "_high")
      if(type == epi_up_name) {
        peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]], gene_id %in% gene)
        up_peak2 <- as.data.frame(peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        y <- as.data.frame(up_peak3)
      }
      if(type == epi_down_name) {
        peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]], gene_id %in% gene)
        down_peak2 <- as.data.frame(peak)
        down_peak2$Row.names <- paste0(down_peak2$seqnames,":",down_peak2$start,"-",down_peak2$end)
        down_peak2 <- down_peak2 %>% distinct(Row.names, .keep_all = T)
        down_peak3 <- with(down_peak2,GRanges(seqnames,IRanges(start,end)))
        y <- as.data.frame(down_peak3)
      }
      return(y)
    }
  })
  output$download_selected_int_peak = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup,"_table.bed")
    },
    content = function(file){write.table(int_selected_bed(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  observeEvent(input$pair_DEG_result, ({
    updateCollapse(session,id =  "input_collapse_pair_RP", open="KS_panel")
  }))
  
  
  
  ##Integrative trackplot-------
  observeEvent(input$RP_table_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_pair_RP", open="int_Trackplot_panel")
  }))
  output$int_igv_uprange <- renderUI({
    if(!is.null(int_goi_gene_position()) && !is.null(int_goi_promoter_position())){
      y <- int_goi_promoter_position()
      gene_position <- int_goi_gene_position()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("int_igv_uprange","Range (x-axis):",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim <- renderUI({
    numericInput("int_igv_ylim","Max peak intensity (y-axis):", value = 2, min = 0)
  })
  
  int_goi_promoter_position<- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table()[input$RP_table_rows_selected,]$gene_id
      y <- NULL
      if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]])) {
        up_peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]], gene_id %in% gene)
      }
      if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])) {
        down_peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]], gene_id %in% gene)
        y <- as.data.frame(down_peak)
      }
      if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]]) && !is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])) {
        peak <- c(up_peak,down_peak)
        y <- as.data.frame(peak)
      }else{
        if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]])) y <- as.data.frame(up_peak)
        if(!is.null(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]])) y <- as.data.frame(down_peak)
      }
      print(head(y))
      return(y)
    }
  })
  
  int_goi_gene_position <- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      gene <- RP_selected_table()[input$RP_table_rows_selected,]$gene_id
      gene_position <- as.data.frame(subset(gene_position(), gene_id %in% gene))
      return(gene_position)
    }
  })
  
  output$int_trackplot_additional <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("int_trackplot_additional1",label="Select additional bigwig files",choices = bw,multiple=T)
    }else{
      fileInput("int_trackplot_additional1",
                "Select additional bigwig files",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  pre_int_track_additional_files <-reactive({
    if(!is.null(input$int_trackplot_additional1)){
      if(length(list.files("./Volume/")) > 0){
        files <- input$int_trackplot_additional1
        name <- gsub(".+\\/","",input$int_trackplot_additional1)
        names(files) <- gsub("\\..+$", "", name)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$int_trackplot_additional1[, 1])){
          file <- input$int_trackplot_additional1[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$int_trackplot_additional1[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
      }
      return(files)
    }
  })
  int_track_additional_files <- reactive({
    return(bws_ordering(bws=pre_int_track_additional_files(),sample_order=input$sample_order_pair_track_int))
  })
  
  int_data_track <- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      return(data_trac(y=int_goi_promoter_position(),gene_position=int_goi_gene_position(),
                       gen=ref(),txdb=txdb(),org=org1(),filetype=input$data_file_type,
                       bw_files=bws_order(),bam_files=bam(),
                       track_additional_files=int_track_additional_files()))
    }
  })
  
  int_highlight_trackplot <- reactive({
    if(!is.null(input$RP_table_rows_selected) && !is.null(input$int_igv_uprange)){
      library(Gviz)
      y <- int_goi_promoter_position()
      gene_position <- int_goi_gene_position()
      chr <- gene_position$seqnames
      df <- int_data_track()
      start <-c()
      width <- c()
      col <- c()
      fill <- c()
      for(i in 1:length(rownames(y))){
        if(y[i,]$start < input$int_igv_uprange[2] && y[i,]$end > input$int_igv_uprange[1]){
          start <- c(start, y[i,]$start)
          width <- c(width, y[i,]$width)
          print(y[i,]$Group)
          if(y[i,]$Group == paste0(unique(collist_bw_pair())[2], "_high")) {
            col <- c(col,"#E41A1C")
            fill <- c(fill,"#FBB4AE")
          }
          if(y[i,]$Group == paste0(unique(collist_bw_pair())[1], "_high")) {
            col <- c(col,"#377EB8")
            fill <- c(fill,"#B3CDE3")
          }
        }
      }
      if(length(start) != 0){
        ht <- HighlightTrack(trackList = df,
                             start = start, width = width,col = col,fill = fill,
                             chromosome = chr, alpha=0.5) 
      }else ht <- NULL
      return(ht)
    }
  })
  int_goi_trackplot <- reactive({
    if(!is.null(input$RP_table_rows_selected) &&
       !is.null(int_goi_promoter_position()) && 
       !is.null(int_goi_gene_position()) && 
       !is.null(gtrack()) &&
       !is.null(input$int_igv_uprange)){
      library(Gviz)
      if(!is.null(int_highlight_trackplot())){
        plot<- plotTracks(list(gtrack(), int_highlight_trackplot()),
                          from = input$int_igv_uprange[1], 
                          to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                          type="hist")
      }else{
        df <- int_data_track()
        df[["gtrack"]] <- gtrack()
        plot<- plotTracks(df,
                          from = input$int_igv_uprange[1], 
                          to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                          type="hist")
      }
      return(plot)
    }
  })
  output$int_trackplot_goi <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$RP_table_rows_selected)){
        int_goi_trackplot()
      }
    })
  })
  output$download_pair_int_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(int_highlight_trackplot())){
          plotTracks(list(gtrack(), int_highlight_trackplot()),
                     from = input$int_igv_uprange[1], 
                     to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                     type="hist")
        }else{
          df <- int_data_track()
          df[["gtrack"]] <- gtrack()
          plotTracks(df,
                     from = input$int_igv_uprange[1], 
                     to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                     type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  #Integrative functional analysis--------
  int_Hallmark_set <- reactive({
    if(!is.null(input$intGeneset)){
      return(GeneList_for_enrichment(Species = input$Species, Gene_set = input$intGeneset, org = org1()))
    }
  })
  
  order_for_intGroup <- reactive({
    if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
      order <- c("RP < -1","-1 < RP < -0.1","-0.1 < RP < -0.01","-0.01 < RP < 0","0 < RP < 0.01","0.01 < RP < 0.1","0.1 < RP < 1","1 < RP")
    }else{
      order <- c("0 < RP < 0.01","0.01 < RP < 0.1","0.1 < RP < 1","1 < RP < 2","2 < RP")
    }
    order <- factor(order,levels = order,ordered = TRUE)
    return(order)
  })
  output$intGroup_for_RPstatus <- renderUI({
    if(input$withRNAseq_pair_enrich_type == "boxplot"){
      if(input$RNAseq_mode_option == "fcw_separate" || input$RNAseq_mode_option == "Classical"){
        selectInput("intGroup_for_RPstatus","Peak",names(pari_RNAseq_boxplot()[["genelist"]]),multiple = F)
      }
    }
  })
  output$intGroup <- renderUI({
    if(!is.null(RP_all_table())){
      if(input$withRNAseq_pair_enrich_type == "boxplot"){
        selectInput("intGroup","Group",order_for_intGroup(),selected = order_for_intGroup(),multiple = T)
      }else{
        group <- unique(pre_RP_selected_table()$Group)
        selectInput("intGroup","Group (Epigenome:RNAseq)",group,selected = group,multiple = T)
      }
    }
  })
  output$intGeneset <- renderUI({
    selectInput('intGeneset', 'Gene Set', gene_set_list)
  })
  withRNAseq_enrichment_analysit_genelist <- reactive({
    if(input$withRNAseq_pair_enrich_type == "boxplot"){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        df <- pari_RNAseq_boxplot()[["genelist"]]
      }else{
        df <- pari_RNAseq_boxplot()[["genelist"]][[input$intGroup_for_RPstatus]]
      }
      colnames(df)[2] <- "Group"
      df <- df %>% dplyr::filter(Group != "Others")
      return(df)
    }
  })
  selected_int_group <- reactive({
    group <- input$intGroup
    df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
    colnames(df) <- c("ENTREZID","Group")
    for(name in group){
      if(input$withRNAseq_pair_enrich_type == "boxplot"){
        table <- withRNAseq_enrichment_analysit_genelist()
      }else{
        table <- pre_RP_selected_table()
      }
      table <- table %>% dplyr::filter(Group == name)
      if(dim(table)[1] != 0){
        entrezid <- table$gene_id
        if(str_detect(table$gene_id[1], "FBgn")){
          my.symbols <- gsub("\\..*","", table$gene_id)
          gene_IDs<-AnnotationDbi::select(org(input$Species),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","ENTREZID"))
          colnames(gene_IDs) <- c("gene_id","ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(ENTREZID, .keep_all = T)
          gene_IDs <- na.omit(gene_IDs)
          table <- merge(table,gene_IDs,by="gene_id")
          entrezid <- table$ENTREZID
        }
        df2 <- data.frame(ENTREZID = entrezid, Group = table$Group)
        df <- rbind(df,df2)
      }
    }
    return(df)
  })
  int_enrich <- reactive({
    return(enrich_viewer_forMulti2(data3 = selected_int_group(), Species = input$Species, org = org1(),
                                   H_t2g = int_Hallmark_set(),Gene_set = input$intGeneset))
  })
  int_enrich_list <- reactive({
    return(enrich_gene_list(data = selected_int_group(),
                            Gene_set = input$intGeneset, org = org1(), H_t2g = int_Hallmark_set()))
  })
  int_enrich_plot <- reactive({
    return(enrich_genelist(data = selected_int_group(),type = "withRNAseq",
                           enrich_gene_list = int_enrich_list(),group_order = input$intGroup))
  })
  output$int_enrichment1 <- renderPlot({
    dotplot_for_output(data = int_enrich(),
                       plot_genelist = int_enrich_plot(), Gene_set = input$intGeneset, 
                       Species = input$Species)
  })
  int_enrich_table <- reactive({
    return(enrich_for_table(data = as.data.frame(int_enrich()), H_t2g = int_Hallmark_set(), Gene_set = input$intGeneset))
  })
  output$int_enrichment_result <- DT::renderDataTable({
    int_enrich_table()
  })
  output$download_pair_int_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$intGeneset,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = int_enrich(),
                           plot_genelist = int_enrich_plot(), Gene_set = input$intGeneset, 
                           Species = input$Species)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_pair_int_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$intGeneset,".txt")
    },
    content = function(file){write.table(int_enrich_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  #withRNAseq combined heatmap-------------
  output$Group_integrated_heatmap <- renderUI({
    if(!is.null(pre_RP_selected_table())){
      selectInput("Group_integrated_heatmap","Group (Epigenome:RNAseq)",unique(pre_RP_selected_table()$Group),multiple = T)
    }
  })
  with_uniqueID_DAR<- reactive({
    group <- input$Group_integrated_heatmap
    list <- list()
    if(is.null(group)) validate("Select groups of interest.")
    for(name in group){
      table <- pre_RP_selected_table() %>% dplyr::filter(Group == name)
      gene <- table$gene_id
      if(input$Genomic_region == "Promoter"){
        tss <- promoters(genes(txdb()),upstream = 0,downstream = 1)
        peak <- subset(tss, gene_id %in% gene)
      }else{
        peak_type <- gsub("\\:.+$", "", name)
        peak <- subset(mmAnno_pair()[[peak_type]], gene_id %in% gene)
      }
      locus <- as.data.frame(peak)
      locus$locus <- paste0(locus$seqnames,":",locus$start,"-",locus$end)
      names(peak) <- paste0(locus$gene_id,"_",locus$locus)
      list[[gsub(":",":\n",name)]] <- peak
    }
    return(list)
  })
  with_integrate_h <- reactive({
    h <- batch_heatmap(files2 = with_uniqueID_DAR(),files_bw = bws_order(),type=input$Genomic_region,withRNAseq=TRUE)
    return(h)
  })
  with_integrated_legend <- reactive({
    lgd <- lgd(files2 = with_uniqueID_DAR(),withRNAseq=TRUE)
    return(lgd)
  })
  with_integrated_heatlist <- reactive({
    if(input$with_integrated_heatmapButton > 0 && updateCounter_with_int$i > 0){
      ht_list <- NULL
      if(!is.null(with_integrate_h())) ht_list <- ht_list + with_integrate_h()[["heatmap"]]
      if(!is.null(with_integrated_heatmap_add1())) ht_list <- ht_list + with_integrated_heatmap_add1()[["heatmap"]]
      if(!is.null(with_integrated_heatmap_add2())) ht_list <- ht_list + with_integrated_heatmap_add2()[["heatmap"]]
      if(!is.null(with_integrated_heatmap_add3())) ht_list <- ht_list + with_integrated_heatmap_add3()[["heatmap"]]
      if(!is.null(with_rnaseq_DEGs2())) ht_list <- ht_list + with_rnaseq_DEGs_heatmap()
      if(!is.null(with_rnaseq_count2())) ht_list <- ht_list + with_rnaseq_count_heatmap()
      return(ht_list)
    }else return(NULL)
  })
  updateCounter_with_int <- reactiveValues(i = 0)
  
  observe({
    input$with_integrated_heatmapButton
    isolate({
      updateCounter_with_int$i <- updateCounter_with_int$i + 1
    })
  })
  #Restart
  observeEvent(with_integrated_heatlist(), {
    isolate(updateCounter_with_int$i == 0)
    updateCounter_with_int <<- reactiveValues(i = 0)
  }) 
  with_integrated_heatlist_plot <- reactive({
    if(input$with_integrated_heatmapButton > 0 && !is.null(bws()) && 
       !is.null(with_integrated_heatlist())){
      return(draw(with_integrated_heatlist(),annotation_legend_list = list(with_integrated_legend()),
                  heatmap_legend_side = "bottom", ht_gap = unit(2, "mm")))
    }
  })
  output$with_integrated_heatmap <- renderPlot({
    if(input$data_file_type != "Row2"){
      if(is.null(bws())) validate("Bigwig files are required.")
    }
    if(input$with_integrated_heatmapButton > 0 && !is.null(bws()) && 
       !is.null(with_integrated_heatlist())){
      print(with_integrated_heatlist_plot())
    }
  })
  with_gene_type_rnaseq_DEGs <- reactive({
    files <- list()
    files[["deg"]] <- RNAseqDEG()
    return(gene_type_for_integrated_heatmap(files=files,Species=input$Species,org=org1()))
  })
  with_rnaseq_DEGs2 <- reactive({
    files <- list()
    files[["deg"]] <- RNAseqDEG()
    return(rnaseqDEGs_for_integrated_heatmap(files = files,Species=input$Species,gene_type=with_gene_type_rnaseq_DEGs()))
  })
  with_gene_type_rnaseq_counts <- reactive({
    if(input$RNAseq_data_type != "Result"){
      files <- list()
      files[["count"]] <- RNAseqDEG_norm()
      print(files)
      return(gene_type_for_integrated_heatmap(files=files,Species=input$Species,org=org1()))
    }
  })
  with_rnaseq_count2 <- reactive({
    if(input$RNAseq_data_type != "Result"){
      files <- list()
      files[["count"]] <- RNAseqDEG_norm()
      print(with_gene_type_rnaseq_counts())
      return(rnaseqCounts_for_integrated_heatmap(files = files,Species=input$Species,gene_type=with_gene_type_rnaseq_counts()))
    }
  })
  with_rnaseq_count_heatmap <- reactive({
    rna <-  with_rnaseq_count2()
    if(!is.null(rna)){
      group <- input$Group_integrated_heatmap
      list <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
      for(name in group){
        table <- pre_RP_selected_table() %>% dplyr::filter(Group == name)
        gene <- table$gene_id
        if(input$Genomic_region == "Promoter"){
          tss <- promoters(genes(txdb()),upstream = 0,downstream = 1)
          peak <- subset(tss, gene_id %in% gene)
        }else{
          type <- gsub("\\:.+$","", name)
          if(type == paste0(unique(collist_bw_pair())[2], "_high")) {
            peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]], gene_id %in% gene)
          }
          if(type == paste0(unique(collist_bw_pair())[1], "_high")) {
            peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]], gene_id %in% gene)
          }
        }
        locus <- as.data.frame(peak)
        locus$locus <- paste0(locus$seqnames,":",locus$start,"-",locus$end)
        list <-rbind(list, locus)
      }
      rna$gene_id <- rownames(rna)
      print(head(rna))
      m_z <- merge(rna,list,by="gene_id",all=T)
      m_z <- m_z %>% dplyr::filter(!is.na(locus))
      rownames(m_z) <- paste0(m_z$gene_id,"_",m_z$locus)
      m_z <- m_z[,2:length(colnames(rna))]
      m_z[is.na(m_z)] <- 0
      print(head(m_z))
      mat <- with_integrate_h()[["mat"]]
      
      cond <- gsub(".+\\.", "", colnames(m_z))
      cond <- gsub("\\_.+$", "", cond)
      cond <- factor(cond, levels = unique(cond), ordered = TRUE)
      cond_color <- rainbow_hcl(length(unique(cond)),c=100)
      names(cond_color) <- unique(cond)
      if(length(names(rnaseq_count())) == 1){
        file_name <- NULL
        file_name_color <- NULL
      }else{
        file_name <- gsub("\\..+$", "", colnames(m_z))
        file_name <- factor(file_name, levels = unique(file_name), ordered = TRUE)
        file_name_color <- rainbow_hcl(length(file_name))
        names(file_name_color) <- file_name
      }
      withProgress(message = "Heatmap of RNA-seq count data",{
        ht <- Heatmap(as.matrix(m_z)[rownames(mat),],name = "RNA-seq\nz-score", column_order = colnames(as.matrix(m_z)[rownames(mat),]),
                      top_annotation = HeatmapAnnotation(files = file_name, condition = cond,
                                                         col = list(files = file_name_color,
                                                                    condition = cond_color)),
                      show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "cm"),use_raster = TRUE)
      })
      return(ht)
    }
  })
  with_rnaseq_DEGs_heatmap <- reactive({
    rna <-  with_rnaseq_DEGs2()
    group <- input$Group_integrated_heatmap
    list <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    for(name in group){
      table <- pre_RP_selected_table() %>% dplyr::filter(Group == name)
      gene <- table$gene_id
      if(input$Genomic_region == "Promoter"){
        tss <- promoters(genes(txdb()),upstream = 0,downstream = 1)
        peak <- subset(tss, gene_id %in% gene)
      }else{
        type <- gsub("\\:.+$","", name)
        if(type == paste0(unique(collist_bw_pair())[2], "_high")) {
          peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[2],"_high")]], gene_id %in% gene)
        }
        if(type == paste0(unique(collist_bw_pair())[1], "_high")) {
          peak <- subset(mmAnno_pair()[[paste0(unique(collist_bw_pair())[1],"_high")]], gene_id %in% gene)
        }
      }
      locus <- as.data.frame(peak)
      locus$locus <- paste0(locus$seqnames,":",locus$start,"-",locus$end)
      list <-rbind(list,locus)
    }
    rna$gene_id <- rownames(rna)
    m_z <- merge(rna,list,by="gene_id",all=T)
    m_z <- m_z %>% dplyr::filter(!is.na(locus))
    rownames(m_z) <- paste0(m_z$gene_id,"_",m_z$locus)
    m_z[is.na(m_z)] <- 0
    m_z <- as.data.frame(m_z[,1:length(colnames(rna))])
    for(i in 1:length(colnames(rna))){
      m_z[,i] <- -1 * as.numeric(m_z[,i])
    }
    m_z[m_z > 5]<- 5
    m_z[m_z < -5]<- -5
    print(head(m_z))
    mat <- with_integrate_h()[["mat"]]
    colnames(m_z) <- gsub("\\.log2F.+$", "", colnames(m_z))
    withProgress(message = "Heatmap of RNA-seq log2FoldChange",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),2:length(colnames(rna))],name = "RNA-seq\nlog2FC", 
                    show_row_names = FALSE, width = unit(2.5, "cm"), column_names_gp = grid::gpar(fontsize = 9),
                    use_raster = TRUE,column_names_side = "top",show_column_dend = FALSE,
                    col = c("blue","white","gold"))
      return(ht)
    })
  })
  output$with_integrated_bw1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("with_integrated_bw_1",
                  "Option: Select additional bigwig files (blue)",
                  choices = bw,multiple=T)
    }else{
      fileInput("with_integrated_bw_1",
                "Option: Select additional bigwig files (blue)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$with_integrated_bw2 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("with_integrated_bw_2",
                  "Option: Select additional bigwig files (green)",
                  choices = bw,multiple=T)
    }else{
      fileInput("with_integrated_bw_2",
                "Option: Select additional bigwig files (green)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$with_integrated_bw3 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("with_integrated_bw_3",
                  "Option: Select additional bigwig files (purple)",
                  choices = bw,multiple=T)
    }else{
      fileInput("with_integrated_bw_3",
                "Option: Select additional bigwig files (purple)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  with_pre_integrated_additional1 <-reactive({
    return(pre_integrated_additional(input$with_integrated_bw_1))
  })
  with_integrated_additional1 <- reactive({
    return(bws_ordering(bws=with_pre_integrated_additional1(),sample_order=input$with_sample_order_pair_comb1))
  })
  with_pre_integrated_additional2 <-reactive({
    return(pre_integrated_additional(input$with_integrated_bw_2))
  })
  with_integrated_additional2 <- reactive({
    return(bws_ordering(bws=with_pre_integrated_additional2(),sample_order=input$with_sample_order_pair_comb2))
  })
  with_pre_integrated_additional3 <-reactive({
    return(pre_integrated_additional(input$with_integrated_bw_3))
  })
  with_integrated_additional3 <- reactive({
    return(bws_ordering(bws=with_pre_integrated_additional3(),sample_order=input$with_sample_order_pair_comb3))
  })
  with_integrated_heatmap_add1 <- reactive({
    if(!is.null(with_integrated_additional1())){
      h <- batch_heatmap(files2 = with_uniqueID_DAR(),files_bw = with_integrated_additional1(),
                         color = c("white","blue"),signal = "blue",type=input$Genomic_region)
      return(h)
    }
  })
  with_integrated_heatmap_add2 <- reactive({
    if(!is.null(with_integrated_additional2())){
      h <- batch_heatmap(files2 = with_uniqueID_DAR(),files_bw = with_integrated_additional2(),
                         color = c("white","green"),signal = "green",type=input$Genomic_region)
      return(h)
    }
  })
  with_integrated_heatmap_add3 <- reactive({
    if(!is.null(with_integrated_additional3())){
      h <- batch_heatmap(files2 = with_uniqueID_DAR(),files_bw = with_integrated_additional3(),
                         color = c("white", "purple"),signal = "purple",type=input$Genomic_region)
      return(h)
    }
  })
  output$download_with_integrated_heatmap = downloadHandler(
    filename = "withRNAseq_conbined_heatmap.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(with_integrated_heatlist_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #combined heatmap--------------
  uniqueID_DAR<- reactive({
    if(input$Genomic_region == "Promoter"){
      data <- data_degcount2()
      data <- data %>% dplyr::filter(!is.na(ENTREZID))
      tss <- as.data.frame(promoters(genes(txdb()),upstream = 0,downstream = 1))
      tss$locus <- paste0(tss$seqnames,":",tss$start,"-",tss$end)
      colnames(tss)[6]<-"ENTREZID"
      up <- dplyr::filter(tss, ENTREZID %in% dplyr::filter(data, group == paste0(unique(collist_bw_pair())[2],"_high"))$ENTREZID)
      down <- dplyr::filter(tss, ENTREZID %in% dplyr::filter(data, group == paste0(unique(collist_bw_pair())[1],"_high"))$ENTREZID)
      print(head(up))
    }else{
      data <- data_degcount2_anno()
      data <- data %>% dplyr::filter(!is.na(ENTREZID))
      up <- dplyr::filter(data, locus %in% rownames(data_degcount_up()))
      down <- dplyr::filter(data, locus %in% rownames(data_degcount_down()))
    }
    up_gr <- with(up, GRanges(seqnames = seqnames, 
                              ranges = IRanges(start,end),
                              strand = strand))
    down_gr <- with(down, GRanges(seqnames = seqnames, 
                                  ranges = IRanges(start,end),
                                  strand = strand))
    names(up_gr) <- paste0(up$ENTREZID,"_",up$locus)
    names(down_gr) <- paste0(down$ENTREZID,"_",down$locus)
    list <- list()
    list[[paste0(unique(collist_bw_pair())[2],"_high")]] <- up_gr
    list[[paste0(unique(collist_bw_pair())[1],"_high")]] <- down_gr
    print(list)
    return(list)
  })
  integrate_h <- reactive({
    h <- batch_heatmap(files2 = uniqueID_DAR(),files_bw = bws_order(),type=input$Genomic_region)
    return(h)
  })
  integrated_legend <- reactive({
    lgd <- lgd(files2 = uniqueID_DAR())
    return(lgd)
  })
  
  integrated_heatlist <- reactive({
    if(input$integrated_heatmapButton > 0 && updateCounter_int$i > 0){
      ht_list <- NULL
      if(!is.null(integrate_h())) ht_list <- ht_list + integrate_h()[["heatmap"]]
      if(!is.null(integrated_heatmap_add1())) ht_list <- ht_list + integrated_heatmap_add1()[["heatmap"]]
      if(!is.null(integrated_heatmap_add2())) ht_list <- ht_list + integrated_heatmap_add2()[["heatmap"]]
      if(!is.null(integrated_heatmap_add3())) ht_list <- ht_list + integrated_heatmap_add3()[["heatmap"]]
      if(!is.null(rnaseq_DEGs2())) ht_list <- ht_list + rnaseq_DEGs_heatmap()
      if(!is.null(rnaseq_count2())) ht_list <- ht_list + rnaseq_count_heatmap()
      return(ht_list)
    }else return(NULL)
  })
  updateCounter_int <- reactiveValues(i = 0)
  
  observe({
    input$integrated_heatmapButton
    isolate({
      updateCounter_int$i <- updateCounter_int$i + 1
    })
  })
  
  
  #Restart
  observeEvent(integrated_heatlist(), {
    isolate(updateCounter_int$i == 0)
    updateCounter_int <<- reactiveValues(i = 0)
  }) 
  integrated_heatlist_plot <- reactive({
    if(input$integrated_heatmapButton > 0 && !is.null(bws()) && !is.null(deg_result()) && 
       !is.null(integrated_heatlist())){
      return(draw(integrated_heatlist(),annotation_legend_list = list(integrated_legend()),
                  heatmap_legend_side = "bottom", ht_gap = unit(2, "mm")))
    }
  })
  output$integrated_heatmap <- renderPlot({
    if(input$data_file_type != "Row2"){
      if(is.null(bws())) validate("Bigwig files are required.")
    }
    if(input$integrated_heatmapButton > 0 && !is.null(bws()) && !is.null(deg_result()) && 
       !is.null(integrated_heatlist())){
      integrated_heatlist_plot()
    }
  })
  output$download_integrated_heatmap = downloadHandler(
    filename = "Integrated_heatmap.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(integrated_heatlist_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$rnaseq_count <- renderUI({
    if(input$Species != "not selected"){
      fileInput("pair_rnaseq_count",
                "Select RNA-seq normalized count files",
                accept = c("txt","csv","xlsx"),
                multiple = TRUE,
                width = "90%")
    }
  })
  output$rnaseq_DEGs <- renderUI({
    if(input$Species != "not selected"){
      fileInput("pair_rnaseq_DEGs",
                "Select RNA-seq DEG result files",
                accept = c("txt","csv","xlsx"),
                multiple = TRUE,
                width = "90%")
    }
  })
  rnaseq_count <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      tmp <- input$pair_rnaseq_count
      upload = list()
      if(is.null(input$pair_rnaseq_count) && input$goButton > 0 )  {
        tmp = "data/RNAseq_count.txt"
        upload[["rna"]]<- read_df(tmp)
        return(upload)
      }else if(is.null(tmp)) {
        return(NULL)
      }else{
        return(read_dfs(tmp))
      }
    })
  })
  rnaseq_DEGs <- reactive({
    withProgress(message = "Importing DEG result files, please wait",{
      tmp <- input$pair_rnaseq_DEGs
      upload = list()
      if(is.null(input$pair_rnaseq_DEGs) && input$goButton > 0 )  {
        tmp = "data/RNAseq.txt"
        upload[["rna"]]<- read_df(tmp)
        return(upload)
      }else if(is.null(tmp)) {
        return(NULL)
      }else{
        return(read_dfs(tmp))
      }
    })
  })
  output$pre_zscoring <- renderUI({
    if(length(names(rnaseq_count())) > 1){
      selectInput("pre_zscoring","pre_zscoring for multiple normalized count data",choices = c(TRUE,FALSE),selected = TRUE)
    }
  })
  gene_type_rnaseq_DEGs <- reactive({
    files <- rnaseq_DEGs()
    return(gene_type_for_integrated_heatmap(files=files,Species=input$Species,org=org1()))
  })
  rnaseq_DEGs2 <- reactive({
    files <- rnaseq_DEGs()
    return(rnaseqDEGs_for_integrated_heatmap(files = files,Species=input$Species,gene_type=gene_type_rnaseq_DEGs()))
  })
  gene_type_rnaseq_counts <- reactive({
    files <- rnaseq_count()
    return(gene_type_for_integrated_heatmap(files=files,Species=input$Species,org=org1()))
  })
  rnaseq_count2 <- reactive({
    files <- rnaseq_count()
    return(rnaseqCounts_for_integrated_heatmap(files = files,Species=input$Species,gene_type=gene_type_rnaseq_counts(),pre_zscoring = input$pre_zscoring))
  })
  observeEvent(input$pair_rnaseq_DEGs, ({
    updateCollapse(session,id =  "z-score_count", open="Uploaded_DEGs")
  }))
  observeEvent(input$pair_rnaseq_count, ({
    updateCollapse(session,id =  "z-score_count", open="z-score_multiple_count_panel")
  }))
  output$rnaseq_count_output <- renderDataTable({
    if(input$Species != "not selected" && !is.null(rnaseq_count())){
      rnaseq_count2()
    }
  })
  output$rnaseq_DEGs_output <- renderDataTable({
    if(input$Species != "not selected" && !is.null(rnaseq_DEGs())){
      rnaseq_DEGs2()
    }
  })
  
  output$Normalized_Count_matrix <- DT::renderDataTable({
    deg_norm_count()
  })
  rnaseq_count_heatmap <- reactive({
    rna <-  rnaseq_count2()
    if(input$Genomic_region == "Promoter"){
      data <- data_degcount2()
      data <- data %>% dplyr::filter(!is.na(ENTREZID))
      tss <- as.data.frame(promoters(genes(txdb()),upstream = 0,downstream = 1))
      tss$locus <- paste0(tss$seqnames,":",tss$start,"-",tss$end)
      colnames(tss)[6]<-"ENTREZID"
      up <- dplyr::filter(tss, ENTREZID %in% dplyr::filter(data, group == paste0(unique(collist_bw_pair())[2],"_high"))$ENTREZID)
      down <- dplyr::filter(tss, ENTREZID %in% dplyr::filter(data, group == paste0(unique(collist_bw_pair())[1],"_high"))$ENTREZID)
    }else{
      data <- data_degcount2_anno() %>% dplyr::filter(!is.na(ENTREZID))
      up <- dplyr::filter(data, locus %in% rownames(data_degcount_up()))
      down <- dplyr::filter(data, locus %in% rownames(data_degcount_down()))
    }
    rna$ENTREZID <- rownames(rna)
    up_m <- merge(rna,up,by="ENTREZID",all=T)
    up_m <- up_m %>% dplyr::filter(!is.na(locus))
    rownames(up_m) <- paste0(up_m$ENTREZID,"_",up_m$locus)
    up_m <- up_m[,2:length(colnames(rna))]
    up_m[is.na(up_m)] <- 0
    
    down_m <- merge(rna,down,by="ENTREZID",all=T)
    down_m <- down_m %>% dplyr::filter(!is.na(locus))
    rownames(down_m) <- paste0(down_m$ENTREZID,"_",down_m$locus)
    down_m <- down_m[,2:length(colnames(rna))]
    down_m[is.na(down_m)] <- 0
    
    m_z <- rbind(up_m,down_m)
    mat <- integrate_h()[["mat"]]
    
    cond <- gsub(".+\\.", "", colnames(m_z))
    cond <- gsub("\\_.+$", "", cond)
    cond <- factor(cond, levels = unique(cond), ordered = TRUE)
    cond_color <- rainbow_hcl(length(unique(cond)),c=100)
    names(cond_color) <- unique(cond)
    if(length(names(rnaseq_count())) == 1){
      file_name <- NULL
      file_name_color <- NULL
    }else{
      file_name <- gsub("\\..+$", "", colnames(m_z))
      file_name <- factor(file_name, levels = unique(file_name), ordered = TRUE)
      file_name_color <- rainbow_hcl(length(file_name))
      names(file_name_color) <- file_name
    }
    print(head(m_z))
    withProgress(message = "Heatmap of RNA-seq count data",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),],name = "RNA-seq\nz-score", 
                    top_annotation = HeatmapAnnotation(files = file_name, condition = cond,
                                                       col = list(files = file_name_color,
                                                                  condition = cond_color)),
                    show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "cm"),use_raster = TRUE)
    })
    return(ht)
  })
  rnaseq_DEGs_heatmap <- reactive({
    rna <-  rnaseq_DEGs2()
    if(input$Genomic_region == "Promoter"){
      data <- data_degcount2()
      data <- data %>% dplyr::filter(!is.na(ENTREZID))
      data <- data[,- which(colnames(data) == "log2FoldChange")]
      tss <- as.data.frame(promoters(genes(txdb()),upstream = 0,downstream = 1))
      tss$locus <- paste0(tss$seqnames,":",tss$start,"-",tss$end)
      colnames(tss)[6]<-"ENTREZID"
      up <- dplyr::filter(tss, ENTREZID %in% dplyr::filter(data, group == paste0(unique(collist_bw_pair())[2],"_high"))$ENTREZID)
      down <- dplyr::filter(tss, ENTREZID %in% dplyr::filter(data, group == paste0(unique(collist_bw_pair())[1],"_high"))$ENTREZID)
    }else{
      data <- data_degcount2_anno() %>% dplyr::filter(!is.na(ENTREZID))
      data <- data[,- which(colnames(data) == "log2FoldChange")]
      up <- dplyr::filter(data, locus %in% rownames(data_degcount_up()))
      down <- dplyr::filter(data, locus %in% rownames(data_degcount_down()))
    }
    
    rna$ENTREZID <- rownames(rna)
    up_m <- merge(rna,up,by="ENTREZID",all=T)
    up_m <- up_m %>% dplyr::filter(!is.na(locus))
    rownames(up_m) <- paste0(up_m$ENTREZID,"_",up_m$locus)
    
    down_m <- merge(rna,down,by="ENTREZID",all=T)
    down_m <- down_m %>% dplyr::filter(!is.na(locus))
    rownames(down_m) <- paste0(down_m$ENTREZID,"_",down_m$locus)
    m_z <- rbind(up_m,down_m)
    m_z[is.na(m_z)] <- 0
    m_z <- as.data.frame(m_z[,1:length(colnames(rna))])
    for(i in 1:length(colnames(rna))){
      m_z[,i] <- -1 * as.numeric(m_z[,i])
    }
    m_z[m_z > 5]<- 5
    m_z[m_z < -5]<- -5
    print(head(m_z))
    mat <- integrate_h()[["mat"]]
    colnames(m_z) <- gsub("\\.log2F.+$", "", colnames(m_z))
    withProgress(message = "Heatmap of RNA-seq log2FoldChange",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),2:length(colnames(rna))],name = "RNA-seq\nlog2FC", 
                    show_row_names = FALSE, width = unit(2.5, "cm"), column_names_gp = grid::gpar(fontsize = 9),
                    use_raster = TRUE,column_names_side = "top",show_column_dend = FALSE,
                    col = c("blue","white","gold"))
      return(ht)
    })
  })
  
  output$integrated_bw1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_1",
                  "Option: Select additional bigwig files (blue)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_1",
                "Option: Select additional bigwig files (blue)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw2 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_2",
                  "Option: Select additional bigwig files (green)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_2",
                "Option: Select additional bigwig files (green)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw3 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_3",
                  "Option: Select additional bigwig files (purple)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_3",
                "Option: Select additional bigwig files (purple)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  pre_integrated_additional1 <-reactive({
    return(pre_integrated_additional(input$integrated_bw_1))
  })
  integrated_additional1 <- reactive({
    return(bws_ordering(bws=pre_integrated_additional1(),sample_order=input$sample_order_pair_comb1))
  })
  pre_integrated_additional2 <-reactive({
    return(pre_integrated_additional(input$integrated_bw_2))
  })
  integrated_additional2 <- reactive({
    return(bws_ordering(bws=pre_integrated_additional2(),sample_order=input$sample_order_pair_comb2))
  })
  pre_integrated_additional3 <-reactive({
    return(pre_integrated_additional(input$integrated_bw_3))
  })
  integrated_additional3 <- reactive({
    return(bws_ordering(bws=pre_integrated_additional3(),sample_order=input$sample_order_pair_comb3))
  })
  integrated_heatmap_add1 <- reactive({
    if(!is.null(integrated_additional1())){
      h <- batch_heatmap(files2 = uniqueID_DAR(),files_bw = integrated_additional1(),
                         color = c("white","darkblue"),signal = "darkblue",type=input$Genomic_region)
      return(h)
    }
  })
  integrated_heatmap_add2 <- reactive({
    if(!is.null(integrated_additional2())){
      h <- batch_heatmap(files2 = uniqueID_DAR(),files_bw = integrated_additional2(),
                         color = c("white","darkgreen"),signal = "green",type=input$Genomic_region)
      return(h)
    }
  })
  integrated_heatmap_add3 <- reactive({
    if(!is.null(integrated_additional3())){
      h <- batch_heatmap(files2 = uniqueID_DAR(),files_bw = integrated_additional3(),
                         color = c("white", "purple"),signal = "purple",type=input$Genomic_region)
      return(h)
    }
  })
  
  ##Pair-wise Report---------
  input_list_data_pair <- reactive({
    if(input$data_file_type == "Row1" || input$data_file_type == "Row1_count"){
      if(!is.null(bws())) bw_files = as.data.frame(names(bws_order())) else bw_files = NA
    }else{
      bw_files = as.data.frame(names(bam()))
    }
    print("1")
    if(!is.null(peak_call_files())) peak_files = as.data.frame(names(peak_call_files())) else peak_files = NA
    print("2")
    if(input$data_file_type == "Row1" || input$data_file_type == "Row2"){
      if(input$Genomic_region == "Genome-wide"){
        list <- data.frame(bw = bw_files[,1], peak_call = peak_files[,1])
        if(input$data_file_type == "Row2") colnames(list)[1] <- "Bam"
      }else{
        list <- data.frame(bw = bw_files[,1])
      }
    }
    if(input$data_file_type == "Row1_count"){
      if(!is.null(bws())) list <- data.frame(bw = bw_files[,1])
    }else list <- NULL
    return(list)
  })
  
  output$pair_report = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"pairwiseDAR_report",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download, please wait for a few minutes.",{
        fs <- c()
        dir.create("DAR_result",showWarnings = FALSE)
        dir.create("Input",showWarnings = FALSE)
        dir.create("Clustering/",showWarnings = FALSE)
        DEG <- "DAR_result/DAR_result.txt"
        up <- paste0("DAR_result/",unique(collist_bw_pair())[2],"_high.bed")
        down <- paste0("DAR_result/",unique(collist_bw_pair())[1],"_high.bed")
        count <- "Input/count.txt"
        bed <- "Input/filtered_merged_peak_call.bed"
        PCA <- "Clustering/clustering.pdf"
        PCA_table <- "Clustering/pca.txt"
        
        fs <- c(DEG, up,down,count,bed,PCA,PCA_table)
        print(fs)
        if(input$Species != "not selected") process_num <- 9 else process_num <- 3
        pdf(PCA, height = 3.5, width = 9)
        print(PCAplot(data = deg_norm_count(),plot=TRUE,legend=input$PCA_legend))
        dev.off()
        write.table(PCAplot(data = deg_norm_count(),plot=FALSE), PCA_table, row.names = T, sep = "\t", quote = F)
        write.table(deg_norm_count(), count, row.names = T, sep = "\t", quote = F)
        if(!is.null(input_list_data_pair())){
          input_list <- "Input/input_list.txt"
          fs <- c(fs, input_list)
          write.table(input_list_data_pair(), input_list, row.names = F, sep = "\t", quote = F)
        }
        withProgress(message = "DAR_result",{
          if(input$Genomic_region == "Promoter" || input$Species == "not selected"){
            write.table(deg_result(), DEG, row.names = T, sep = "\t", quote = F)
          }else{
            write.table(deg_result2(), DEG, row.names = T, sep = "\t", quote = F) 
          }
          write.table(data_degcount_up_bed(), up, row.names = F, col.names = F,sep = "\t", quote = F)
          write.table(data_degcount_down_bed(), down, row.names = F, col.names = F,sep = "\t", quote = F)
          if(input$data_file_type == "Row1_count" && input$Genomic_region == "Genome-wide"){
            data <- range_changer(bws_count())
            data <- with(data, GRanges(seqnames = chr,ranges = IRanges(start=start,
                                                                       end=end)))
          }else data <- promoter_region()
          write.table(as.data.frame(data), bed, row.names = F,col.names=F, sep = "\t", quote = F)
        })
        incProgress(1/process_num)
        if(!is.null(input$peak_up_range) && !is.null(bws())){
          withProgress(message = "Peak pattern",{
            dir.create("peak_pattern",showWarnings = FALSE)
            up_pattern <- paste0("peak_pattern/",unique(collist_bw_pair())[2],"_high_lineplot.pdf")
            down_pattern <- paste0("peak_pattern/",unique(collist_bw_pair())[1],"_high_lineplot.pdf")
            up_pattern_heatmap <- paste0("peak_pattern/",unique(collist_bw_pair())[2],"_high_heatmap.pdf")
            down_pattern_heatmap <- paste0("peak_pattern/",unique(collist_bw_pair())[1],"_high_heatmap.pdf")
            fs <- c(fs,up_pattern,down_pattern, up_pattern_heatmap,down_pattern_heatmap)
            pdf(up_pattern, height = 5, width = 5)
            matplot(peak_up_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                    type="l",ylab="density",lty=1,xaxt="n",)
            axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
            legend("topright", legend=colnames(peak_up_alinedHeatmap()[["line"]]), col=1:6,
                   lty=1, lwd=1)
            dev.off()
            incProgress(1/4)
            pdf(down_pattern, height = 5, width = 5)
            matplot(peak_down_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                    type="l",ylab="density",lty=1,xaxt="n")
            axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
            legend("topright", legend=colnames(peak_down_alinedHeatmap()[["line"]]), col=1:6,
                   lty=1, lwd=1)
            dev.off()
            incProgress(1/4)
            pdf(up_pattern_heatmap, height = 6, width = 6)
            print(plot_grid(peak_up_alinedHeatmap()[["heat"]]))
            dev.off()
            incProgress(1/4)
            pdf(down_pattern_heatmap, height = 6, width = 6)
            print(plot_grid(peak_down_alinedHeatmap()[["heat"]]))
            dev.off()
            incProgress(1/4)
          })
        }
        incProgress(1/process_num)
        if(!is.null(input$xrange)){
          withProgress(message = "volcano plot",{
            volcano <- "DAR_result/volcano_plot.pdf"
            fs <- c(fs, volcano)
            pdf(volcano, height = 4, width = 4)
            print(pair_volcano())
            dev.off()
          })
        }
        incProgress(1/process_num)
        if(input$Species != "not selected"){
          if(input$Genomic_region == "Genome-wide"){
            dir.create("peak_distribution",showWarnings = FALSE)
            up_distribution <- paste0("peak_distribution/",unique(collist_bw_pair())[2],"_high_annotation.pdf")
            up_distribution_table <- paste0("peak_distribution/",unique(collist_bw_pair())[2],"_high_annotation.txt")
            down_distribution <- paste0("peak_distribution/",unique(collist_bw_pair())[1],"_high_annotation.pdf")
            down_distribution_table <- paste0("peak_distribution/",unique(collist_bw_pair())[1],"_high_annotation.txt")
            fs <- c(fs, up_distribution,up_distribution_table,down_distribution,down_distribution_table)
            withProgress(message = "Peak distribution",{
              
              df <- as.data.frame(updistribution()$peaks)
              df$Group <- paste0(unique(collist_bw_pair())[2],"_high")
              write.table(df, up_distribution_table, row.names = F,col.names=T, sep = "\t", quote = F)
              df <- as.data.frame(downdistribution()$peaks)
              df$Group <- paste0(unique(collist_bw_pair())[1],"_high")
              write.table(df, down_distribution_table, row.names = F,col.names=T, sep = "\t", quote = F)
              
              pdf(up_distribution, height = 6, width = 10)
              print(updistribution()$plot)
              dev.off()
              incProgress(1/2)
              pdf(down_distribution, height = 6, width = 10)
              print(downdistribution()$plot)
              dev.off()
              incProgress(1/2)
            })
          }
          incProgress(1/process_num)
          if(!is.null(input$Gene_set) && input$Genomic_region == "Genome-wide"){
            dir.create("GREAT",showWarnings = FALSE)
            dotplot <- paste0("GREAT/dotplot_",input$Gene_set,".pdf")
            enrich_table <- paste0("GREAT/enrichment_",input$Gene_set,".txt")
            fs <- c(fs, dotplot,enrich_table)
            withProgress(message = "GREAT",{
              write.table(as.data.frame(enrichment_1_1()), enrich_table, row.names = F, sep = "\t", quote = F)
              p1 <- pair_enrich1_H()
              pdf(dotplot, height = 5, width = 6)
              print(plot_grid(p1))
              dev.off()
              if(!is.null(input$Pathway_list) && 
                 !is.null(input$Up_or_Down) && !is.null(enrichment_enricher())){
                region_associate_plot <- paste0("GREAT/",input$Gene_set,"_",input$Pathway_list,"_",input$Up_or_Down,".pdf")
                region_associate_table <- paste0("GREAT/",input$Gene_set,"_",input$Pathway_list,"_",input$Up_or_Down,".txt")
                fs <- c(fs, region_associate_plot,region_associate_table)
                set_list <- input$Pathway_list
                df <- enrichment_enricher()
                name <- input$Up_or_Down
                write.table(region_gene_associate(), region_associate_table, row.names = F, sep = "\t", quote = F)
                pdf(region_associate_plot, height = 5, width = 12)
                rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
                dev.off()
              }
            })
          }
          incProgress(1/process_num)
          if(input$motifButton > 0 && !is.null(enrich_motif()) && 
             !is.null(input$homer_unknown)){
            withProgress(message = "HOMER",{
              path_list <- enrich_motif()
              base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
              for(name in names(path_list)){
                files <-list.files(path_list[[name]],pattern = "*.*")
                for(i in 1:length(files)){
                  data <- paste0(path_list[[name]],"/",files[[i]])
                  fs <- c(fs, data)
                }
              }
              motif <- paste0(base_dir,"/homer_dotplot",".pdf")
              p1 <- homer_Motifplot(df = enrich_motif(),showCategory = input$homer_showCategory)
              fs <- c(fs, motif)
              pdf(motif, height = 6, width = 7)
              print(p1)
              dev.off()
            })
          }else {
            base_dir <- NULL
          }
          incProgress(1/process_num)
          if(!is.null(RNAseqDEG()) && !is.na(input$DEG_fdr) && 
             !is.na(input$DEG_fc)){
            withProgress(message = "withRNAseq",{
              dirname_withRNA <- paste0("withRNAseq-",input$RNAseq_mode,"_range-",input$peak_distance,"kb_fc",input$DEG_fc,"_fdr",input$DEG_fdr,"_RNAseq-",input$pair_DEG_result$name,"/")
              dir.create(dirname_withRNA,showWarnings = FALSE)
              RNAseq_barplot <- paste0(dirname_withRNA,"barplot.pdf")
              fs <- c(fs, RNAseq_barplot)
              pdf(RNAseq_barplot, height = 5, width = 10)
              gridExtra::grid.arrange(RNAseq_popu(), ChIPseq_popu(), ncol = 1)
              dev.off()
              dir.create(paste0(dirname_withRNA,"KSplot"),showWarnings = FALSE)
              KSplot_table <- paste0(dirname_withRNA,"KSplot/KSplot.txt")
              fs <- c(fs,KSplot_table)
              write.table(ks_tables(),KSplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
              if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
                KSplot <- paste0(dirname_withRNA,"KSplot/KSplot.pdf")
                fs <- c(fs,KSplot)
                pdf(KSplot, height = 5, width = 7)
                print(regulatory_potential())
                dev.off()
              }else{
                KSplot_up <- paste0(dirname_withRNA,"KSplot/",RNAseq_name()[2],".pdf")
                KSplot_down <- paste0(dirname_withRNA,"KSplot/",RNAseq_name()[1],".pdf")
                fs <- c(fs,KSplot_up,KSplot_down)
                down <- regulatory_potential()[[RNAseq_name()[1]]] + ggtitle(RNAseq_name()[1])
                up <- regulatory_potential()[[RNAseq_name()[2]]] + ggtitle(RNAseq_name()[2])
                pdf(KSplot_up, height = 5, width = 7)
                print(up)
                dev.off()
                pdf(KSplot_down, height = 5, width = 7)
                print(down)
                dev.off()
              }
              print("ksplot finish")
              dir.create(paste0(dirname_withRNA,"boxplot"),showWarnings = FALSE)
              boxplot_table <- paste0(dirname_withRNA,"boxplot/boxplot.txt")
              boxplot <- paste0(dirname_withRNA,"boxplot/boxplot.pdf")
              fs <- c(fs,boxplot_table,boxplot)
              if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined") pdf_width <- 7 else pdf_width <- 12
              pdf(boxplot, height = 5, width = pdf_width)
              print(pari_RNAseq_boxplot()[["plot"]])
              dev.off()
              write.table(pari_RNAseq_boxplot()[["statistical_test"]],boxplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
              
              
              dir.create(paste0(dirname_withRNA,"boxplot/bed/"),showWarnings = FALSE)
              if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
                boxplot_genelist <- paste0(dirname_withRNA,"boxplot/","boxplot_genelist.txt")
                fs <- c(fs,boxplot_genelist)
                write.table(pari_RNAseq_boxplot()[["genelist"]],boxplot_genelist,col.names = T,row.names = F,sep = "\t",quote = F)
                for(name in names(pari_RNAseq_boxplot()[["bedlist"]])){
                  bed <- pari_RNAseq_boxplot()[["bedlist"]][[name]]
                  bed_name <- paste0(dirname_withRNA,"boxplot/bed/",name,".bed")
                  bed_name <- gsub(" < ","....",bed_name)
                  fs <- c(fs, bed_name)
                  write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
                }
              }else{
                dir.create(paste0(dirname_withRNA,"boxplot/gene_list/"),showWarnings = FALSE)
                genelist_up <- paste0(dirname_withRNA,"boxplot/gene_list/",RNAseq_name()[2],".txt")
                genelist_down <- paste0(dirname_withRNA,"boxplot/gene_list/",RNAseq_name()[1],".txt")
                up <- pari_RNAseq_boxplot()[["genelist"]][[RNAseq_name()[2]]]
                down <- pari_RNAseq_boxplot()[["genelist"]][[RNAseq_name()[1]]]
                fs <- c(fs,genelist_up,genelist_down)
                write.table(up,genelist_up,col.names = T,row.names = F,sep = "\t",quote = F)
                write.table(down,genelist_down,col.names = T,row.names = F,sep = "\t",quote = F)
                
                for(name in names(pari_RNAseq_boxplot()[["bedlist"]])){
                  dir.create(paste0(dirname_withRNA,"boxplot/bed/",name),showWarnings = FALSE)
                  for(file in names(pari_RNAseq_boxplot()[["bedlist"]][[name]])){
                    bed <- pari_RNAseq_boxplot()[["bedlist"]][[name]][[file]]
                    bed_name <- paste0(dirname_withRNA,"boxplot/bed/",name,"/",file,".bed")
                    bed_name <- gsub(" < ","....",bed_name)
                    fs <- c(fs, bed_name)
                    write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
                  }
                }
              }
              
              print("boxplot finish")
              
              
              dir.create(paste0(dirname_withRNA,"RP_table"),showWarnings = FALSE)
              RP_summary <- paste0(dirname_withRNA,"RP_table/summary.txt")
              fs <- c(fs,RP_summary)
              write.table(pre_RP_selected_table(), RP_summary, row.names = F, sep = "\t", quote = F)
              print("RP_summary finish")
              dir.create(paste0(dirname_withRNA,"RP_table/selected_table(epigenome--RNA)"),showWarnings = FALSE)
              dir.create(paste0(dirname_withRNA,"RP_table/selected_bed(epigenome--RNA)"),showWarnings = FALSE)
              for(name in unique(pre_RP_selected_table()$Group)){
                RP_selected <- paste0(dirname_withRNA,"RP_table/selected_table(epigenome--RNA)/",name,".txt")
                RP_selected_bed <- paste0(dirname_withRNA,"RP_table/selected_bed(epigenome--RNA)/",name,".bed")
                RP_selected <- gsub(":","--",RP_selected)
                RP_selected_bed <- gsub(":","--",RP_selected_bed)
                fs <- c(fs,RP_selected,RP_selected_bed)
                table <- pre_RP_selected_table() %>% dplyr::filter(Group == name)
                write.table(table, RP_selected, row.names = F, sep = "\t", quote = F)
                print("RP_selected finish")
                gene <- table$gene_id
                peak <- gsub("\\:.+$", "", name)
                y <- NULL
                if(!is.null(mmAnno_pair())) {
                  up_peak <- subset(mmAnno_pair()[[peak]], gene_id %in% gene)
                  up_peak2 <- as.data.frame(up_peak)
                  up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
                  up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
                  up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
                  mcols(up_peak3) <- DataFrame(Group = name)
                  y <- as.data.frame(up_peak3)
                }
                write.table(y, RP_selected_bed, row.names = F, col.names = F,sep = "\t", quote = F)
              }
              
              if(!is.null(input$RP_table_rows_selected) &&
                 !is.null(int_goi_promoter_position()) && 
                 !is.null(int_goi_gene_position()) && 
                 !is.null(input$int_igv_uprange) && !is.null(bws())){
                print(RP_selected_table()[input$RP_table_rows_selected,])
                gene <- RP_selected_table()[input$RP_table_rows_selected,]$Symbol
                inttrack <- paste0(dirname_withRNA,gene,".pdf")
                fs <- c(fs, inttrack)
                pdf(inttrack, height = 4, width = 7)
                if(!is.null(int_highlight_trackplot())){
                  plotTracks(list(gtrack(), int_highlight_trackplot()),
                             from = input$int_igv_uprange[1], 
                             to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                             type="hist")
                }else{
                  df <- int_data_track()
                  df[["gtrack"]] <- gtrack()
                  plotTracks(df,
                             from = input$int_igv_uprange[1], 
                             to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                             type="hist")
                }
                dev.off()
              }
              if(!is.null(input$intGeneset) && !is.null(input$intGroup)){
                dir.create(paste0(dirname_withRNA,"enrichment_analysis/"),showWarnings = FALSE)
                intdotplot <- paste0(dirname_withRNA,"enrichment_analysis/dotplot_",input$intGeneset,".pdf")
                intenrichtable <- paste0(dirname_withRNA,"enrichment_analysis/enrichment_",input$intGeneset,".txt")
                fs <- c(fs, intdotplot,intenrichtable)
                pdf(intdotplot, height = 5, width = 8)
                dotplot_for_output(data = int_enrich(),
                                   plot_genelist = int_enrich_plot(), Gene_set = input$intGeneset, 
                                   Species = input$Species)
                dev.off()
                write.table(int_enrich_table(), intenrichtable, row.names = F, sep = "\t", quote = F)
              }
              if(!is.null(input$Group_integrated_heatmap) && !is.null(with_integrated_heatlist_plot())){
                dir.create(paste0(dirname_withRNA,"combined_heatmap/"),showWarnings = FALSE)
                intcomheat <- paste0(dirname_withRNA,"combined_heatmap/combined_heatmap.pdf")
                fs <- c(fs, intcomheat)
                pdf(intcomheat, height = 6, width = 10)
                with_integrated_heatlist_plot()
                dev.off()
              }
              if(input$with_motifButton > 0 && !is.null(with_enrich_motif()) && 
                 !is.null(input$with_homer_unknown)){
                withProgress(message = "HOMER",{
                  with_path_list <- with_enrich_motif()
                  with_base_dir <- gsub("\\/.+$", "", path_list[[names(with_path_list)[1]]])
                  for(name in names(with_path_list)){
                    files <-list.files(with_path_list[[name]],pattern = "*.*")
                    for(i in 1:length(files)){
                      data <- paste0(with_path_list[[name]],"/",files[[i]])
                      fs <- c(fs, data)
                    }
                  }
                  with_motif <- paste0(dirname_withRNA,"/homer_dotplot",".pdf")
                  p1 <- homer_Motifplot(df = with_enrich_motif(),showCategory = input$with_homer_showCategory)
                  fs <- c(fs, with_motif)
                  pdf(with_motif, height = 6, width = 7)
                  print(p1)
                  dev.off()
                })
              }else {
                with_base_dir <- NULL
              }
            })
          }else {
            dirname_withRNA <- NULL
            gene <- NULL
          }
          incProgress(1/process_num)
          if(input$integrated_heatmapButton > 0 && !is.null(bws()) && !is.null(deg_result()) && 
             !is.null(integrated_heatlist())){
            dir.create("combined_heatmap",showWarnings = FALSE)
            intheatmap <- paste0("combined_heatmap/","combined_heatmap.pdf")
            fs <- c(fs, intheatmap)
            pdf(intheatmap, height = 6, width = 10)
            print(integrated_heatlist_plot())
            dev.off()
          }
        }else {
          base_dir <- NULL
          dirname_withRNA <- NULL
          gene <- NULL
        }
        report_name <- paste0(format(Sys.time(), "%Y%m%d_"),"pairwise_report",".docx")
        fs <- c(fs, report_name)
        tempReport <- file.path(tempdir(),"pair_report.Rmd")
        file.copy("pair_report.Rmd", tempReport, overwrite = TRUE)
        rmarkdown::render("pair_report.Rmd", output_format = "word_document", output_file = report_name,
                          params = list(input = input,
                                        collist_bw_pair = collist_bw_pair(),
                                        deg_norm_count = deg_norm_count(),
                                        input_list_data_pair = input_list_data_pair(),
                                        enrichment_1_1 = enrichment_1_1(),
                                        region_gene_associate = region_gene_associate(),
                                        enrich_motif = enrich_motif(),
                                        base_dir = base_dir,
                                        RNAseqDEG = RNAseqDEG(),
                                        dirname_withRNA = dirname_withRNA,
                                        RP_all_table = RP_all_table(),
                                        int_goi_promoter_position = int_goi_promoter_position(),
                                        int_goi_gene_position = int_goi_gene_position(),
                                        gene = gene,
                                        int_enrich_table = int_enrich_table(),
                                        integrated_heatlist = integrated_heatlist(),
                                        RNAseq_file = RNAseq_file(),
                                        with_integrated_heatlist_plot = with_integrated_heatlist_plot()), 
                          envir = new.env(parent = globalenv()),intermediates_dir = tempdir(),encoding="utf-8"
        )
      })
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  #Venn diagram -----------
  output$Spe_venn_distribution <- renderText({
    if(input$Species_venn == "not selected") print("Please select 'Species'")
  })
  observeEvent(input$goButton_venn,({
    updateSelectInput(session,inputId = "Species_venn","Species",species_list, selected = "Homo sapiens (hg19)")
  }))
  observeEvent(pre_bws_venn(),({
    updateSelectizeInput(session,inputId = "sample_order_venn","Sample order:",
                         choices =  names(pre_bws_venn()),
                         selected = names(pre_bws_venn()))
  }))
  output$sample_order_venn_comb2 <- renderUI({
    if(!is.null(pre_integrated_additional2_venn())){
      selectInput(inputId = "sample_order_venn_comb2","Sample order (blue):",
                  choices =  names(pre_integrated_additional2_venn()),
                  selected = names(pre_integrated_additional2_venn()),multiple = TRUE)
    }
  })
  output$sample_order_venn_comb3 <- renderUI({
    if(!is.null(pre_integrated_additional3_venn())){
      selectInput(inputId = "sample_order_venn_comb3","Sample order (green):",
                  choices =  names(pre_integrated_additional3_venn()),
                  selected = names(pre_integrated_additional3_venn()),multiple = TRUE)
    }
  })
  output$sample_order_venn_comb4 <- renderUI({
    if(!is.null(pre_integrated_additional4_venn())){
      selectInput(inputId = "sample_order_venn_comb4","Sample order (purple):",
                  choices =  names(pre_integrated_additional4_venn()),
                  selected = names(pre_integrated_additional4_venn()),multiple = TRUE)
    }
  })
  ##-----------
  output$peak_call_file_venn1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bed <- c(str_subset(list,".narrowPeak"),str_subset(list,".bed"))
      selectInput("peak_call_file_venn1",label=NULL,choices = bed,multiple=T)
    }else{
      fileInput("peak_call_file_venn1",
                NULL,
                accept = c("bed","narrowPeak"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$file_venn1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("file_venn1",label=NULL,choices = bw,multiple=T)
    }else{
      fileInput("file_venn1",
                NULL,
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$Spe_motif_venn <- renderText({
    if(input$Species_venn == "not selected") print("Please select 'Species'")
  })
  output$Spe_GREAT_venn <- renderText({
    if(input$Species_venn == "not selected") print("Please select 'Species'")
  })
  
  txdb_venn <- reactive({
    if(input$Species_venn != "not selected"){
      return(txdb_function(Species = input$Species_venn))
    }
  })
  org_venn <- reactive({
    return(org(Species = input$Species_venn))
  })
  gene_position_venn <- reactive({
    if(input$Species_venn != "not selected"){
      return(genes(txdb_venn()))
    }
  })
  
  updateCounter_vennStart <- reactiveValues(i = 0)
  
  observe({
    input$vennButton
    isolate({
      updateCounter_vennStart$i <- updateCounter_vennStart$i + 1
    })
  })
  
  
  #Restart
  observeEvent(Venn_peak_call_files(), {
    isolate(updateCounter_vennStart$i == 0)
    updateCounter_vennStart <<- reactiveValues(i = 0)
  }) 
  #venn-------
  Venn_peak_call_files <- reactive({
    if(is.null(input$peak_call_file_venn1)){
      if(input$goButton_venn > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      if(length(list.files("./Volume/")) > 0){
        files <- input$peak_call_file_venn1
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        names(files2) <- bed_name(input$peak_call_file_venn1)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$peak_call_file_venn1[, 1])){
          file <- input$peak_call_file_venn1[[nr, 'datapath']]
          name <- c(name, bed_name(input$peak_call_file_venn1[nr,]$name))
          files <- c(files,file)
        }
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        names(files2)<-name
      }
      return(files2)
    }
  })
  
  pre2_bws_venn <- reactive({
    if(is.null(input$file_venn1)){
      if(input$goButton_venn > 0 ){
        df<-list()
        df[["A_1"]] <- "data/bigwig/A_1.BigWig"
        df[["A_2"]] <- "data/bigwig/A_2.BigWig"
        df[["B_1"]] <- "data/bigwig/B_1.BigWig"
        df[["B_2"]] <- "data/bigwig/B_2.BigWig"
        return(df)
      }
      return(NULL)
    }else{
      if(length(list.files("./Volume/")) > 0){
        files <- input$file_venn1
        names(files) <- bigwig_name(input$file_venn1)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file_venn1[, 1])){
          file <- input$file_venn1[[nr, 'datapath']]
          name <- c(name, bigwig_name(input$file_venn1[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
      }
      return(files)
    }
  })
  pre_bws_venn <- debounce(pre2_bws_venn, 1000)
  bws_venn <- reactive({
    return(bws_ordering(bws=pre_bws_venn(),sample_order=input$sample_order_venn,additional=FALSE))
  })
  venn_overlap <- reactive({
    if(!is.null(Venn_peak_call_files()) && updateCounter_vennStart$i > 0){
      withProgress(message = "Preparing intersection, takes a few minutes",{
        ol <- ChIPpeakAnno::findOverlapsOfPeaks(Venn_peak_call_files(),connectedPeaks = "keepAll")
        names(ol$peaklist) <- gsub("///","-",names(ol$peaklist))
        return(ol)
      })
    }
  })
  make_venn <-reactive({
    df<-list()
    for(name in names(Venn_peak_call_files())){
      data <- as.data.frame(Venn_peak_call_files()[[name]])
      data$strand <- as.character(data$strand)
      df[[name]] <- data
    }
    return(df)
  })
  
  output$venn <- renderPlot({
    withProgress(message = "Venn diagram",{
      if(!is.null(Venn_peak_call_files()) && updateCounter_vennStart$i > 0){
        ggVennPeaks(make_venn(),label_size = 5, alpha = .2)
      }
    })
  })
  
  output$download_vennplot = downloadHandler(
    filename = function(){
      paste0("Venn_diagram",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(ggVennPeaks(make_venn(),label_size = 5, alpha = .2))
        dev.off()
        incProgress(1)
      })
    }
  )
  observeEvent(bws_venn(), ({
    updateCollapse(session,id =  "Lineplot of intersections", open="vennLineplot_intersection_panel")
  }))
  
  venn_batch_lineplot <- reactive({
    withProgress(message = "Line plot",{
      return(batch_lineplot(files2 = venn_overlap()$peaklist,files_bw = bws_venn()))
    })
  })
  
  output$venn_batch_bed_line <- renderPlot({
    if(!is.null(Venn_peak_call_files()) && !is.null(bws_venn())){
      if(!is.null(venn_overlap())){
        venn_batch_lineplot()[["bed"]]
      }
    }
  })
  output$venn_batch_bigwig_line <- renderPlot({
    if(!is.null(Venn_peak_call_files()) && !is.null(bws_venn())){
      if(!is.null(venn_overlap())){
        venn_batch_lineplot()[["bigwig"]]
      }
    }
  })
  
  output$download_venn_intersection_line = downloadHandler(
    filename = function(){
      paste0("Lineplot_intersection",".pdf")
    },
    content = function(file) {
      rowlist <- length(names(venn_overlap()$peaklist))
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)+3
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)+3
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(venn_batch_lineplot()[["bed"]])
        dev.off()
        incProgress(1)
      })
    }
  )
  output$venn_report = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"venn_report",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download, please wait for >5 minutes.",{
        fs <- c()
        print(fs)
        dir.create("intersection_bed",showWarnings = FALSE)
        dir.create("Input",showWarnings = FALSE)
        venn <- "venn.pdf"
        input_list <- "Input/input_bed_list.txt"
        fs <- c(venn,input_list)
        pdf(venn, height = 6, width = 6)
        print(ggVennPeaks(make_venn(),label_size = 5, alpha = .2))
        dev.off()
        write.table(input_list_data_venn()[["bed"]], input_list, row.names = F, sep = "\t", quote = F)
        for(name in names(venn_overlap()$peaklist)){
          intersect_bed <- paste0("intersection_bed/",name,".bed")
          fs <- c(fs, intersect_bed)
          venn_g <- venn_overlap()$peaklist[[name]]
          write.table(as.data.frame(venn_g)[1:3], 
                      intersect_bed, row.names = F, col.names = F,sep = "\t", quote = F)
        }
        process_num <- length(names(venn_overlap()$peaklist))
        if(input$Species_venn != "not selected"){
          dir.create("annotation",showWarnings = FALSE)
          distribution <- paste0("annotation","/distribution.pdf")
          fs <- c(fs, distribution)
          pdf(distribution, height = 4.5, width = 6)
          print(vendistribution())
          dev.off()
          for(name in names(vendistribution()$peaks)){
            annotation <- paste0("annotation/",name,".txt")
            annotation_pdf <- paste0("annotation/",name,".pdf")
            fs <- c(fs, annotation,annotation_pdf)
            pdf(annotation_pdf, height = 4.5, width = 6)
            print(donut_replot(dplyr::filter(vendistribution()$plot[[1]], source == name)))
            dev.off()
            write.table(apply(selected_annoData_table()[[name]][,1:10],2,as.character), 
                        annotation, row.names = F, col.names = T,sep = "\t", quote = F)
          }
        }
        line_plot_all_bed <- paste0("lineplot_bed.pdf")
        fs <- c(fs, line_plot_all_bed)
        rowlist <- length(names(venn_overlap()$peaklist))
        pdf(line_plot_all_bed, height = pdf_h(rowlist)+3, width = pdf_w(rowlist)+3)
        print(venn_batch_lineplot()[["bed"]])
        dev.off()
        
        if(!is.null(bws_venn())){
          line_plot_all_bw <- paste0("lineplot_bigwig.pdf")
          fs <- c(fs, line_plot_all_bw)
          rowlist <- length(names(bws_venn()))
          pdf(line_plot_all_bw, height = pdf_h(rowlist)+3, width = pdf_w(rowlist)+3)
          print(venn_batch_lineplot()[["bigwig"]])
          dev.off()
          if(!is.null(input$intersect_select2)){
            dir.create(input$intersect_select2,showWarnings = FALSE)
            input_bw_list <- "Input/input_bw_list.txt"
            heatmap <- paste0(input$intersect_select2,"/heatmap.pdf")
            lineplot <- paste0(input$intersect_select2,"/lineplot.pdf")
            fs <- c(fs, heatmap, lineplot,input_bw_list)
            write.table(input_list_data_venn()[["bw"]], input_bw_list, row.names = F, sep = "\t", quote = F)
            pdf(heatmap, height = 6, width = 6)
            print(plot_grid(peak_venn_alinedHeatmap()[["heat"]]))
            dev.off()
            pdf(lineplot, height = 5, width = 5)
            matplot(peak_venn_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                    type="l",ylab="density",lty=1,xaxt="n")
            axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
            legend("topright", legend=colnames(peak_venn_alinedHeatmap()[["line"]]), col=1:6,
                   lty=1, lwd=1)
            dev.off()
          }
        }
        if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && 
           !is.null(input$homer_unknown_venn) && input$Species_venn != "not selected"){
          path_list <- enrich_motif_venn()
          base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
          for(name in names(path_list)){
            files <-list.files(path_list[[name]],pattern = "*.*")
            for(i in 1:length(files)){
              data <- paste0(path_list[[name]],"/",files[[i]])
              fs <- c(fs, data)
            }
          }
          homer_plot <- paste0(base_dir,"/homer_dotplot.pdf")
          fs <- c(fs, homer_plot)
          pdf(homer_plot, height = 6, width = 8)
          print(venn_motif_plot())
          dev.off()
        }
        if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){
          dir.create("GREAT",showWarnings = FALSE)
          great_dotplot <- paste0("GREAT/dotplot_",input$Gene_set_venn,".pdf")
          great_table <- paste0("GREAT/enrichment_",input$Gene_set_venn,".txt")
          region_associate_plot <- paste0("GREAT/",input$Gene_set,"_",input$Pathway_list,"_",input$Up_or_Down,".pdf")
          great_regionplot <- paste0("GREAT/",input$Gene_set_venn,"_",input$Pathway_list_venn,"_",input$intersection_venn_fun,".pdf")
          great_regionplot_table <- paste0("GREAT/",input$Gene_set_venn,"_",input$Pathway_list_venn,"_",input$intersection_venn_fun,".txt")
          fs <- c(fs, great_dotplot,great_regionplot,great_regionplot_table)
          if(input$venn_pdf_height == 0){
            pdf_height <- 6
          }else pdf_height <- input$venn_pdf_height
          if(input$venn_pdf_width == 0){
            pdf_width <- 8
          }else pdf_width <- input$venn_pdf_width
          pdf(great_dotplot, height = 6, width = 8)
          print(plot_grid(pair_enrich1_H_venn()))
          dev.off()
          write.table(as.data.frame(enrichment_1_1_venn()), great_table, row.names = F, sep = "\t", quote = F)
          set_list <- input$Pathway_list_venn
          df <- enrichment_enricher_venn()
          name <- input$intersection_venn_fun
          pdf(great_regionplot, height = 5, width = 12)
          rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
          dev.off()
          write.table(region_gene_associate_venn(), great_regionplot_table, row.names = F, sep = "\t", quote = F)
        }
        if(!is.null(input$RNAseqGroup_venn) && 
           !is.null(input$peak_distance_venn && !is.null(pre_mmAnno_venn()))){
          dirname <- paste0("withRNAseq-",input$RNAseq_mode_venn,"_range-",input$peak_distance_venn,"kb_fc",input$DEG_fc_venn,"_fdr",input$DEG_fdr_venn,"_RNAseq-",input$venn_DEG_result$name,"/")
          dir.create(dirname,showWarnings = FALSE)
          
          dir.create(paste0(dirname,"boxplot/"),showWarnings = FALSE)
          if(input$RNAseq_data_type_venn != "List"){
            pdf_height <- pdf_h(input$venn_select_RNA)
            pdf_width <- pdf_w(input$venn_select_RNA)
            boxplot_table <- paste0(dirname,"boxplot/","boxplot.txt")
            boxplot <- paste0(dirname,"boxplot/","boxplot.pdf")
            fs <- c(fs,boxplot_table,boxplot)
            pdf(boxplot, height = pdf_height, width = pdf_width)
            print(RNAseq_boxplot_venn()[["plot"]])
            dev.off()
            write.table(RNAseq_boxplot_venn()[["statistical_test"]],boxplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
          }
          dir.create(paste0(dirname,"boxplot/gene_list/"),showWarnings = FALSE)
          dir.create(paste0(dirname,"boxplot/bed/"),showWarnings = FALSE)
          for(name in venn_select_RNA_debounce()){
            genelist <- paste0(dirname,"boxplot/gene_list/",name,".txt")
            up <- RNAseq_boxplot_venn()[["genelist"]][[name]]
            fs <- c(fs,genelist)
            write.table(up,genelist,col.names = T,row.names = F,sep = "\t",quote = F)
          }
          for(name in names(RNAseq_boxplot_venn()[["bedlist"]])){
            dir.create(paste0(dirname,"boxplot/bed/",name),showWarnings = FALSE)
            for(file in names(RNAseq_boxplot_venn()[["bedlist"]][[name]])){
              bed <- RNAseq_boxplot_venn()[["bedlist"]][[name]][[file]]
              bed_name <- paste0(dirname,"boxplot/bed/",name,"/",file,".bed")
              bed_name <- gsub(" < ","....",bed_name)
              fs <- c(fs, bed_name)
              write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
            }
          }
          RNAseq_barplot <- paste0(dirname,"barplot.pdf")
          RP_all <- paste0(dirname,"RP_summary.txt")
          fs <- c(fs,  RNAseq_barplot,RP_all)
          pdf(RNAseq_barplot, height = 7.5, width = 14)
          gridExtra::grid.arrange(RNAseq_popu_venn(), ChIPseq_popu_venn(), ncol = 1)
          dev.off()
          write.table(RP_all_table_venn(), RP_all, row.names = F, sep = "\t", quote = F)
          
          dir.create(paste0(dirname,"KSplot/"),showWarnings = FALSE)
          for(name in names(venn_overlap()$peaklist)){
            ksplot <- paste0(dirname,"KSplot/",name,".pdf")
            fs <- c(fs,ksplot)
            pdf(ksplot, height = 5, width = 7)
            print(regulatory_potential_venn()[[name]])
            dev.off()
          }
          dir.create(paste0(dirname,"selected_table(intersection--RNA)/"),showWarnings = FALSE)
          dir.create(paste0(dirname,"selected_bed(intersection--RNA)/"),showWarnings = FALSE)
          for(name in unique(RP_all_table_venn()$Group)){
            RP_selected <- paste0(dirname,"selected_table(intersection--RNA)/",name,".txt")
            RP_selected_bed <- paste0(dirname,"selected_bed(intersection--RNA)/",name,".bed")
            RP_selected <- gsub(":","--",RP_selected)
            RP_selected_bed <- gsub(":","--",RP_selected_bed)
            fs <- c(fs, RP_selected,RP_selected_bed)
            table <- RP_all_table_venn() %>% dplyr::filter(Group == name)
            write.table(table, RP_selected, row.names = F, sep = "\t", quote = F)
            gene <- table$gene_id
            y <- NULL
            peak <- gsub("\\:.+$", "", name)
            if(!is.null(pre_mmAnno_venn())) {
              up_peak <- subset(pre_mmAnno_venn()[[peak]], gene_id %in% gene)
              up_peak2 <- as.data.frame(up_peak)
              up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
              up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
              up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
              mcols(up_peak3) <- DataFrame(Group = name)
              y <- as.data.frame(up_peak3)
            }
            write.table(y, RP_selected_bed, row.names = F, col.names = F,sep = "\t", quote = F)
          }
          if(!is.null(input$intGeneset_venn) && !is.null(intGroup_venn_debounce())){
            dir.create(paste0(dirname,"enrichment_analysis/"),showWarnings = FALSE)
            intdotplot <- paste0(dirname,"enrichment_analysis/dotplot_",input$intGeneset_venn,".pdf")
            intenrichtable <- paste0(dirname,"enrichment_analysis/enrichment_",input$intGeneset_venn,".txt")
            fs <- c(fs, intdotplot,intenrichtable)
            pdf(intdotplot, height = 5, width = 8)
            dotplot_for_output(data = int_enrich_venn(),
                               plot_genelist = int_enrich_plot_venn(), Gene_set = input$intGeneset_venn, 
                               Species = input$Species_venn)
            dev.off()
            write.table(int_enrich_table_venn(), intenrichtable, row.names = F, sep = "\t", quote = F)
          }
        }
        if(input$integrated_heatmapButton_venn > 0 && !is.null(bws_venn()) && 
           !is.null(integrated_heatlist_venn()) && !is.null(input$venn_heatmap_group)){
          dir.create("combined_heatmap",showWarnings = FALSE)
          intheatmap <- paste0("combined_heatmap/","combined_heatmap.pdf")
          fs <- c(fs, intheatmap)
          pdf(intheatmap, height = 6, width = 10)
          print(integrated_heatlist_venn_plot())
          dev.off()
        }
        incProgress(1)
      })
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  output$download_venn_bigwig_line = downloadHandler(
    filename = function(){
      paste0("Lineplot_bigwig",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        rowlist <- length(names(bws_venn()))
        if(input$venn_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)+3
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)+3
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(venn_batch_lineplot()[["bigwig"]])
        dev.off()
        incProgress(1)
      })
    }
  )
  input_list_data_venn <- reactive({
    peak_files = as.data.frame(names(Venn_peak_call_files()))
    list <- list()
    list[["bed"]] <- data.frame(input_peak_call_files = peak_files[,1])
    if(!is.null(bws_venn())){
      list[["bw"]] <- data.frame(input_bigwig_files = as.data.frame(names(bws_venn()))[,1])
    }
    return(list)
  })
  output$select_file2_venn <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("select_file2_venn", "Select intersect",
                  c(names(venn_overlap()$peaklist)),multiple = F)
    }
  })
  selected_grange_venn <- reactive({
    if(!is.null(input$select_file2_venn)){
      if(input$select_file2_venn != "not selected"){
        return(venn_overlap()$peaklist[[input$select_file2_venn]])
      }
    }
  })
  output$selected_intersect <- DT::renderDT({
    if(!is.null(selected_grange_venn())){
      as.data.frame(selected_grange_venn())
    }
  })
  output$download_venn_all_bed = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"all_bed",".zip")
    },
    content = function(fname){
      fs <- c()
      dir.create("all_bed/",showWarnings = FALSE)
      for(name in names(venn_overlap()$peaklist)){
        file <- paste0("all_bed/",name,".bed")
        fs <- c(fs,file)
        df <- as.data.frame(venn_overlap()$peaklist[[name]])
        write.table(apply(df,2,as.character),file,col.names = F,row.names = F,sep = "\t",quote = F)
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  output$intersect_select <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("intersect_select", "Select intersect",
                  c(names(venn_overlap()$peaklist)),multiple = F)
    }
  })
  output$intersect_select <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("intersect_select", "Select intersect",
                  c(names(venn_overlap()$peaklist)),multiple = F)
    }
  })
  output$select_file2 <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("intersect_select2", "Select intersect",
                  c(names(venn_overlap()$peaklist)),multiple = F)
    }
  })
  
  vendistribution <- reactive({
    return(ChIPpeakAnno::genomicElementDistribution(GRangesList(venn_overlap()$peaklist), 
                                                    TxDb = txdb_venn()))
  })
  output$venn_peak_distribution <- renderPlot({
    withProgress(message = "Peak distribution",{
      if(!is.null(venn_overlap()) && input$Species_venn != "not selected"){
        vendistribution()$plot
      }
    })
  })
  output$selected_venn_peak_distribution <- renderPlot({
    if(!is.null(input$intersect_select) &&
       input$Species_venn != "not selected" ){
      withProgress(message = "Peak distribution",{
        if(!is.null(venn_overlap()) && input$Species_venn != "not selected"){
          donut_replot(dplyr::filter(vendistribution()$plot[[1]], source == input$intersect_select))
        }
      })
    }
  })
  
  
  selected_annoData_table <- reactive({
    withProgress(message = "Preparing annotation",{
      df <- list()
      for(name in names(vendistribution()$peaks)){
        overlaps.anno <- ChIPseeker::annotatePeak(vendistribution()$peaks[[name]],
                                                  TxDb = txdb_venn()) %>% as.data.frame()
        my.symbols <- overlaps.anno$geneId
        gene_IDs<-id_convert(my.symbols,input$Species_venn,type="ENTREZID")
        colnames(gene_IDs) <- c("geneId","NearestGene")
        data <- merge(gene_IDs,overlaps.anno,by="geneId")
        data <- data[,2:11] %>% distinct(peakNames, .keep_all = T) %>% as.data.frame() 
        df[[name]] <- data
      }
      return(df)
    })
  })
  
  output$selected_intersect_annotation <- DT::renderDT({
    if(!is.null(input$intersect_select) &&
       input$Species_venn != "not selected" ){
      selected_annoData_table()[[input$intersect_select]] %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  output$download_venn_peak_distribution = downloadHandler(
    filename = function(){
      paste0("Peak_distribution-",input$intersect_select,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 4.5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(vendistribution())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_selected_intersect_annotation_table <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"all_annotation",".zip")
    },
    content = function(fname){
      fs <- c()
      dir.create("all_annotation/",showWarnings = FALSE)
      for(name in names(selected_annoData_table())){
        file <- paste0("all_annotation/",name,".txt")
        file_pdf <- paste0("all_annotation/",name,".pdf")
        fs <- c(fs,file,file_pdf)
        df <- apply(selected_annoData_table()[[name]],2,as.character)
        pdf(file_pdf, height = 4.5, width = 6)
        print(donut_replot(dplyr::filter(vendistribution()$plot[[1]], source == name)))
        dev.off()
        write.table(df,file,col.names = T,row.names = F,sep = "\t",quote = F)
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  ###Venn Peak pattern comparison--------
  selected_grange2 <- reactive({
    if(!is.null(input$intersect_select2)){
      if(input$intersect_select2 != "not selected"){
        return(venn_overlap()$peaklist[[input$intersect_select2]])
      }
    }
  })
  output$peak_pattern_venn_heat_range <- renderUI({
    if(!is.null(input$intersect_select2) && 
       input$intersect_select2 != "not selected" && !is.null(bws_venn())){
      withProgress(message = "Preparing peak pattern",{
        rg <- venn_pattern_range()
        print(ceiling(rg*2))
        sliderInput("peak_venn_range","Intensity range",value=rg,min = 0,max=ceiling(rg*2),step=ceiling(rg*2)/100)
      })
    }
  })
  venn_pattern_range <- reactive({
    rg <- c()
    sig <- peak_pattern_function(grange=selected_grange2(), files=bws_venn(),plot = FALSE)
    for(name in names(sig)){
      rg <- c(rg, mean(sig[[name]][,50]))
    }
    rg <- max(rg) + max(rg)*0.1
    return(rg)
  })
  
  peak_venn_alinedHeatmap <- reactive({
    if(!is.null(input$peak_venn_range)){
      bigwig <- bigwig_breakline(bws_venn())
      heatmap <- peak_pattern_function(grange=selected_grange2(), files=bigwig,rg = input$peak_venn_range)
      return(heatmap)
    }
  })
  
  output$peak_pattern_venn_heatmap <- renderPlot({
    if(!is.null(input$intersect_select2) && 
       input$intersect_select2 != "not selected" && !is.null(bws_venn())){
      withProgress(message = "feature aligned heatmap",{
        plot_grid(peak_venn_alinedHeatmap()[["heat"]])
      })
    }
  })
  output$peak_pattern_venn_line <- renderPlot({
    if(!is.null(input$intersect_select2) && input$intersect_select2 != "not selected" && 
       !is.null(bws_venn()) && !is.null(input$peak_venn_range)){
      withProgress(message = "feature aligned distribution",{
        matplot(peak_venn_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_venn_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
      })
    }
  })
  output$download_peak_pattern_venn_heatmap = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_venn_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_venn_alinedHeatmap()[["heat"]]))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_venn_line = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_venn_alignedDistibution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        matplot(peak_venn_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_venn_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##venn_trackplot----------
  ref_venn <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species_venn))
    return(ref)
  })
  
  
  
  observeEvent(input$selected_intersect_annotation_rows_selected, ({
    updateCollapse(session,id =  "int_result_collapse_panel", open="Trackplot_venn_panel")
  }))
  observeEvent(input$selected_intersect_annotation_rows_selected,({
    if(input$Species_venn != "not selected"){
      if(!is.null(goi_promoter_position_venn()) && !is.null(goi_gene_position_venn())){
        y <- goi_promoter_position_venn()
        gene_position <- goi_gene_position_venn()
        start_position <- min(c(y$start,gene_position$start))
        end_position <- max(c(y$end,gene_position$end))
        updateSliderInput(session,"igv_venn_uprange","Range (x-axis):",value = c(start_position,end_position),
                          step = 100, min = start_position - 10000, max = end_position + 10000)
      }
    }
  }))
  
  output$igv_venn_uprange <- renderUI({
    if(!is.null(goi_promoter_position_venn()) && !is.null(goi_gene_position_venn())){
      y <- goi_promoter_position_venn()
      gene_position <- goi_gene_position_venn()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("igv_venn_uprange","Range (x-axis):",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$igv_venn_ylim <- renderUI({
    numericInput("igv_venn_ylim","Max peak intensity (y-axis):", value = 2, min = 0)
  })
  
  gtrack_venn <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=0.8)
      return(gtrack)
    })
  })
  
  goi_promoter_position_venn<- reactive({
    if(!is.null(input$selected_intersect_annotation_rows_selected)){
      library(Gviz)
      label_data <- selected_annoData_table()[[input$intersect_select]][input$selected_intersect_annotation_rows_selected,]
      return(label_data)
    }
  })
  
  goi_gene_position_venn <- reactive({
    if(!is.null(goi_promoter_position_venn())){
      label_data <- goi_promoter_position_venn()$NearestGene
      gene_IDs<- id_convert(label_data,input$Species_venn,type="SYMBOL_single")
      gene_position <- as.data.frame(subset(gene_position_venn(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  
  data_track_venn <-reactive({
    return(data_trac(y=goi_promoter_position_venn(),gene_position=goi_gene_position_venn(),
                     gen=ref_venn(),txdb=txdb_venn(),org=org_venn(),
                     bw_files=bws_venn(),track_additional_files=NULL))
  })
  
  highlight_trackplot_venn <- reactive({
    library(Gviz)
    y <- goi_promoter_position_venn()
    gene_position <- goi_gene_position_venn()
    chr <- gene_position$seqnames
    df <- data_track_venn()
    if(y$start < input$igv_venn_uprange[2] && y$end > input$igv_venn_uprange[1]){
      withProgress(message = "Highlighting the selected region",{
        ht <- HighlightTrack(trackList = df,
                             start = y$start, width = y$width,
                             chromosome = chr, alpha=0.5)
      })
    }else ht <- NULL
    return(ht)
  })
  goi_trackplot_venn <- reactive({
    withProgress(message = "Trackplot",{
      if(!is.null(input$selected_intersect_annotation_rows_selected) &&
         !is.null(input$igv_venn_uprange)){
        library(Gviz)
        if(!is.null(highlight_trackplot_venn())){
          plot<- plotTracks(list(gtrack_venn(), highlight_trackplot_venn()),
                            from = input$igv_venn_uprange[1], 
                            to = input$igv_venn_uprange[2],ylim=c(0,input$igv_venn_ylim),
                            type="hist")
        }else{
          df <- data_track_venn()
          df[["gtrack"]] <- gtrack_venn()
          plot<- plotTracks(df,
                            from = input$igv_venn_uprange[1], 
                            to = input$igv_venn_uprange[2],ylim=c(0,input$igv_venn_ylim),
                            type="hist") 
        }
        return(plot)
      }
    })
  })
  output$trackplot_venn_goi <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      goi_trackplot_venn()
    })
  })
  
  output$download_venn_trackplot = downloadHandler(
    filename = "Trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        plotTracks(list(gtrack_venn(), highlight_trackplot_venn()),
                   from = goi_gene_position_venn()$start-input$igv_venn_uprange, 
                   to = goi_gene_position_venn()$end+input$igv_venn_downrange,ylim=c(0,input$igv_venn_ylim),
                   type="hist")
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  #Venn functional analysis------
  output$Gene_set_venn <- renderUI({
    selectInput('Gene_set_venn', 'Gene Set', gene_set_list)
  })
  
  output$venn_whichGroup2 <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      if(!is.null(venn_overlap())){
        selectInput("venn_whichGroup2", "Select intersects", choices = c(names(venn_overlap()$peaklist)),multiple = TRUE)
      }
    }
  })
  
  Hallmark_set_venn <- reactive({
    if(!is.null(input$Gene_set_venn)){
      return(GeneList_for_enrichment(Species = input$Species_venn, Gene_set = input$Gene_set_venn, org = org_venn()))
    }
  })
  
  pre_selected_grange_venn_list <- reactive({
    if(!is.null(input$venn_whichGroup2)){
      Glist <- GRangesList()
      for(name in input$venn_whichGroup2){
        data <- venn_overlap()$peaklist[[name]]
        Glist[[name]] <- data
      }
      return(Glist)
    }
  })
  selected_grange_venn_list <- debounce(pre_selected_grange_venn_list,1000)
  
  enrichment_enricher_venn <- reactive({
    data <- selected_grange_venn_list()
    H_t2g <- Hallmark_set_venn()
    df <- list()
    source <- ref_for_GREAT(input$Species_venn)
    for(name in names(data)){
      res = rGREAT::great(gr = data[[name]],gene_sets = gene_list_for_enrichment_genome(H_t2g,input$Species_venn),source)
      df[[name]] <- res
    }
    return(df)
  })
  
  enrichment_1_1_venn <- reactive({
    df <- enrichment_enricher_venn()
    if(!is.null(input$Gene_set_venn) && input$Species_venn != "not selected" && !is.null(df)){
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                          "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
        }else group1 <- NULL
        if(length(group1$id)!=0){
          group1$Group <- name
          data <- rbind(data, group1)
        }
      }
      if(length(data$id) != 0) {
        return(data)
      }else return(NULL)
    }else{return(NULL)}
  })
  
  pair_enrich1_H_venn <- reactive({
    if(!is.null(input$Gene_set_venn) && input$Species_venn != "not selected"){
      df <- enrichment_enricher_venn()
      data3 <- selected_grange_venn_list()
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          if(length(group1$id)==0){
            group1 <- NULL
          }else{
            group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
            if (length(group1$p_adjust_hyper) > 5){
              group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
              group1 <- group1[1:5,]
            }
          }
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                             "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      if(length(data$Description) != 0){
        data$Group <- gsub("-", "- ", data$Group)
        for(i in 1:length(data$Group)){
          data$Group[i] <- paste(strwrap(data$Group[i], width = 15),collapse = "\n")
        }
        data$Group <- gsub(" \\(", "\n\\(", data$Group)
        data$Description <- gsub("_", " ", data$Description)
        data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
        data$x <- gsub("-","", data$x)
        data <- dplyr::arrange(data, x)
        idx <- order(data[["x"]], decreasing = FALSE)
        data$Description <- factor(data$Description,
                                   levels=rev(unique(data$Description[idx])))
        p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                        geom_point() +
                        scale_color_continuous(low="red", high="blue",
                                               guide=guide_colorbar(reverse=TRUE)) +
                        scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                        scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      }else p1 <- NULL
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment_venn <- renderPlot({
    if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){
      dotplot_for_output(data = selected_grange_venn_list(),
                         plot_genelist = pair_enrich1_H_venn(), Gene_set = input$Gene_set_venn, 
                         Species = input$Species_venn)
    }
  })
  
  
  output$download_venn_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$Gene_set_venn,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H_venn()
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_venn_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$Gene_set_venn,".txt")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1_venn()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  pair_enrich_table_venn <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1_venn()), 
                            H_t2g = Hallmark_set_venn(), Gene_set = input$Gene_set_venn))
  })
  
  output$venn_enrichment_result <- DT::renderDataTable({
    if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){
      as.data.frame(enrichment_1_1_venn())
    }
  })
  
  
  
  output$whichGeneSet_venn <- renderUI({
    if(!is.null(input$intersection_venn_fun) && dim(as.data.frame(enrichment_1_1_venn()))[1] != 0){
      group <- input$intersection_venn_fun
      if(dim(as.data.frame(enrichment_1_1_venn()))[1]!=0){
        set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1_venn()), Group == group)$id)
        selectInput('Pathway_list_venn', 'Pathway list', set_list)
      }
    }
  })
  
  region_gene_associate_venn <- reactive({
    if(!is.null(input$Pathway_list_venn)){
      set_list <- input$Pathway_list_venn
      df <- enrichment_enricher_venn()
      data3 <- selected_grange_venn_list()
      data <- data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))
          if(dim(group1)[1] != 0){
            group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
          }else group1 <- NULL
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      data <- apply(data,2,as.character)
      return(data)
    }
  })
  output$region_gene_venn_associations <- DT::renderDT({
    if(!is.null(input$Pathway_list_venn) && 
       !is.null(input$intersection_venn_fun) && input$Species_venn != "not selected"){
      region_gene_associate_venn()
    }
  })
  output$whichGenes_venn <- renderUI({
    if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){ 
      selectInput('intersection_venn_fun', 'Select intersection', names(enrichment_enricher_venn()))
    }
  })
  
  region_gene_associate_plot_venn <- reactive({
    set_list <- input$Pathway_list_venn
    df <- enrichment_enricher_venn()
    name <- input$intersection_venn_fun
    if(!is.null(df[[name]])) {
      res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
    }else res <- NULL
    return(res)
  })
  output$region_gene_associations_venn_plot <- renderPlot({
    if(!is.null(input$Pathway_list_venn) && 
       !is.null(input$intersection_venn_fun) && input$Species_venn != "not selected"){
      region_gene_associate_plot_venn()
    }
  })
  output$download_venn_region_gene_associations_plot = downloadHandler(
    filename = function(){
      paste0("region_gene_associations-",input$intersection_venn_fun,"-",input$Pathway_list_venn,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$venn_pdf_width
        set_list <- input$Pathway_list_venn
        df <- enrichment_enricher_venn()
        name <- input$intersection_venn_fun
        pdf(file, height = pdf_height, width = pdf_width)
        rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_venn_region_gene_associations = downloadHandler(
    filename = function() {
      paste0("region_gene_associations-",input$intersection_venn_fun,"-",input$Pathway_list_venn,".txt")
    },
    content = function(file){write.table(region_gene_associate_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  #Venn regulatory potential------
  output$intersect_RNA <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("intersect_RNA", "Select intersect",
                  c(names(venn_overlap()$peaklist)),
                  selected = c(names(venn_overlap()$peaklist)),multiple = F)
    }
  })
  output$venn_select_RNA <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      if(length(names(venn_overlap()$peaklist)) > 3) intersects <- c("") else intersects <- c(names(venn_overlap()$peaklist))
      selectInput("venn_select_RNA", "Select intersect",
                  c(names(venn_overlap()$peaklist)),
                  selected = intersects,multiple = T)
    }
  })
  react_venn_select_RNA_debounce <- reactive({
    if(!is.null(input$venn_select_RNA)) return(input$venn_select_RNA) else return(NULL)
  })
  venn_select_RNA_debounce <- debounce(react_venn_select_RNA_debounce, 1000)
  
  peak_venn_grange_RNA  <- reactive({
    if(!is.null(input$intersect_RNA)){
      return(venn_overlap()$peaklist[[input$intersect_RNA]])
    }
  })
  
  
  
  output$vennRNAseqresult <- renderUI({
    if(input$RNAseq_data_type_venn == "Result") {
      label <- "Select RNA-seq DEG result file"
    }else{
      if(input$RNAseq_data_type_venn == "List") label <- "Select Gene list files" else 
        label <- "Select RNA-seq raw count file"
    }
    fileInput("venn_DEG_result",
              label,
              accept = c("txt","csv","xlsx"),
              multiple = TRUE,
              width = "80%")
  })
  output$pre_genelist_input_choice_venn <- renderUI({
    if(!is.null(pre_RNAseq_file_venn())){
      list <- unique(pre_RNAseq_file_venn()[,2])
      if(length(list) > 20 && length(input$venn_DEG_result[, 1]) == 1){ 
        selectInput("pre_genelist_input_choice_venn","Group name",c("File name","Second column"),selected = "File name",multiple = F)
      }else return(NULL)
    }
  })
  output$genelist_input_choice_venn <- renderUI({
    if(!is.null(pre_RNAseq_file_venn())){
      list <- unique(pre_RNAseq_file_venn()[,2])
      if(!is.null(input$pre_genelist_input_choice_venn)){
        if(input$pre_genelist_input_choice_venn == "File name"){
          file_name <- gsub(paste0("\\.",tools::file_ext(input$venn_DEG_result[[1, 'datapath']]),"$"), "", input$venn_DEG_result[1,]$name)
          list <- file_name
        }
      }
      list <- sort(list)
      selectInput("genelist_input_choice_venn","Groups (RNA)",list, selected = list, multiple = TRUE)
    }
  })
  
  pre_RNAseq_file_venn <- reactive({
    withProgress(message = "Importing RNA-seq data, please wait",{
      tmp <- input$venn_DEG_result$datapath 
      if(is.null(input$venn_DEG_result) && input$goButton_venn > 0 )  {
        if(input$RNAseq_data_type_venn == "Result") tmp = "data/RNAseq.txt" else{
          if(input$RNAseq_data_type_venn == "List") tmp = "data/genelist.txt" else tmp = "data/RNAseq_count.txt"
        } 
      }
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(input$RNAseq_data_type_venn == "List"){
          upload = list()
          name = c()
          if(is.null(input$venn_DEG_result) && input$goButton_venn > 0 ){
            upload["genelist"] <- list(read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = ""))
            name <- "genelist"
          }else{
            for(nr in 1:length(input$venn_DEG_result[, 1])){
              if(tools::file_ext(input$venn_DEG_result[[nr, 'datapath']]) == "xlsx") {
                df2 <- read_xlsx(input$venn_DEG_result[[nr, 'datapath']])
                df2 <- as.data.frame(df2)
                df <- try(data.frame(row.names = df2[,1]),silent = T)
                if(class(df) != "try-error") {
                  if(dim(df2)[2] == 2){
                    df <- data.frame(row.names = df2[,1],a = df2[,2])
                    colnames(df)[1] <- colnames(df2)[2]
                  }else{
                    rownames(df2) <- df2[,1]
                    df <- df2[,-1]
                    colnames(df) <- gsub("-",".",colnames(df))
                  }
                }
              }
              if(tools::file_ext(input$venn_DEG_result[[nr, 'datapath']]) == "csv") df <- read.csv(input$venn_DEG_result[[nr, 'datapath']], header=TRUE, sep = ",",quote = "")
              if(tools::file_ext(input$venn_DEG_result[[nr, 'datapath']]) == "txt" || 
                 tools::file_ext(input$venn_DEG_result[[nr, 'datapath']]) == "tsv") df <- read.table(input$venn_DEG_result[[nr, 'datapath']], header=TRUE, sep = "\t",quote = "")
              file_name <- gsub(paste0("\\.",tools::file_ext(input$venn_DEG_result[[nr, 'datapath']]),"$"), "", input$venn_DEG_result[nr,]$name)
              name <- c(name, file_name)
              upload[nr] <- list(df)
            }
          }
          names(upload) <- name
          if(length(names(upload)) == 1){
            tmp <- upload[[name]]
            rownames(tmp) = gsub("\"", "", rownames(tmp))
            if(str_detect(colnames(tmp)[1], "^X\\.")){
              colnames(tmp) = str_sub(colnames(tmp), start = 3, end = -2) 
            }
            print(rownames(tmp)[1])
            print(dim(tmp)[2])
            if(rownames(tmp)[1] == 1){
              if(dim(tmp)[2] >= 2){
                tmp <- data.frame(Gene = tmp[,1], Group = tmp[,2])
              }else{
                tmp <- data.frame(Gene = tmp[,1], 
                                  Group = gsub(paste0("\\.",tools::file_ext(input$venn_DEG_result[[1, 'datapath']]),"$"), "", input$venn_DEG_result[1,]$name))
              }
            }else{
              tmp <- data.frame(Gene = rownames(tmp), Group = tmp[,1])
            }
          }else{
            df2 <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
            for(file in names(upload)){
              df <- upload[[file]]
              if(rownames(df)[1] == 1){
                df[,1] = gsub("\"", "", df[,1])
                df <- data.frame(Gene = df[,1], Group = file)
              }else{
                rownames(df) = gsub("\"", "", rownames(df))
                df <- data.frame(Gene = rownames(df), Group = file)
              }
              df2 <- rbind(df2,df)
            }
            tmp <- df2
          }
          tmp$Group <- gsub(":","-",tmp$Group)
          return(tmp)
        }else{
          if(tools::file_ext(tmp) == "xlsx") {
            df2 <- readxl::read_xlsx(tmp) 
            df2 <- as.data.frame(df2)
            df <- try(data.frame(row.names = df2[,1]),silent = T)
            if(class(df) != "try-error") {
              if(dim(df2)[2] == 2){
                df <- data.frame(row.names = df2[,1],a = df2[,2])
                colnames(df)[1] <- colnames(df2)[2]
              }else{
                rownames(df2) <- df2[,1]
                df <- df2[,-1]
                colnames(df) <- gsub("-",".",colnames(df))
              }
            }
          }
          if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
          if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
          rownames(df) = gsub("\"", "", rownames(df))
          if(length(colnames(df)) != 0){
            if(str_detect(colnames(df)[1], "^X\\.")){
              colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
            }
          }
          return(df)
        }
      }
    })
  })
  output$venn_RNAseq_order <- renderUI({
    if(!is.null(pre_RNAseq_file_venn())){
      if(input$RNAseq_data_type_venn != "Result" && input$RNAseq_data_type_venn != "List"){
        if(length(grep("log2FoldChange", colnames(pre_RNAseq_file_venn()))) != 0) validate("Uploaded data is a DEG_result file. Please Select 'DEG_Result' from 'Input'.")
        selectInput("venn_RNAseq_order","Select samples:",
                    choices = colnames(pre_RNAseq_file_venn()),selected = colnames(pre_RNAseq_file_venn()),multiple = T)
      }
    }
  })
  RNAseq_file_venn <- reactive({
    if(!is.null(pre_RNAseq_file_venn())){
      if(input$RNAseq_data_type_venn != "Result"){
        count <- pre_RNAseq_file_venn()
        order <- input$venn_RNAseq_order
        data <- try(count[,order])
        if(length(data) == 1){
          if(class(data) == "try-error") validate("")
        }
        return(data)
      }
    }
  })
  updateCounter_DEGanalysis_venn <- reactiveValues(i = 0)
  
  observe({
    input$DEGanalysis_Button_venn
    isolate({
      updateCounter_DEGanalysis_venn$i <- updateCounter_DEGanalysis_venn$i + 1
    })
  })
  
  #Restart
  observeEvent(RNAseq_file_venn(), {
    isolate(updateCounter_DEGanalysis_venn$i == 0)
    updateCounter_DEGanalysis_venn <<- reactiveValues(i = 0)
  }) 
  pre_RNAseqDEG_venn <- reactive({
    if(input$DEGanalysis_Button_venn > 0 && updateCounter_DEGanalysis_venn$i > 0){
      count <- RNAseq_file_venn()
      if(!is.null(count)){
        collist <- gsub("\\_.+$", "", colnames(count))
        if(length(unique(collist)) == 2){
          group <- data.frame(con = factor(collist))
          dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
          dds$con <- factor(dds$con, levels = unique(collist))
          dds <- DESeq(dds)
          return(dds)
        }else validate(print(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 2 conditions.")))
      }
    }
  })
  RNAseqDEG_venn <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if(input$RNAseq_data_type_venn == "Result"){
        if(length(grep("log2FoldChange", colnames(pre_RNAseq_file_venn()))) == 0) validate("Error: the uploaded file is in an unexpected format.\nThe file does not contain 'log2FoldChange' or 'padj'.")
        return(pre_RNAseq_file_venn())
      }else{
        if(input$RNAseq_data_type_venn != "List"){
          if(input$DEGanalysis_Button_venn > 0 && updateCounter_DEGanalysis_venn$i > 0){
            count <- RNAseq_file_venn()
            collist <- gsub("\\_.+$", "", colnames(count))
            dds <- pre_RNAseqDEG_venn()
            contrast <- c("con", unique(collist))
            res <- results(dds,  contrast = contrast)
            res <- as.data.frame(res)
            return(res)
          }
        }else {
          if(!is.null(input$genelist_input_choice_venn)){
            RNAdata <- pre_RNAseq_file_venn()
            if(!is.null(input$pre_genelist_input_choice_venn)){
              if(input$pre_genelist_input_choice_venn == "File name"){
                file_name <- gsub(paste0("\\.",tools::file_ext(input$venn_DEG_result[[1, 'datapath']]),"$"), "", input$venn_DEG_result[1,]$name)
                RNAdata$Group <- file_name
              }
            }
            group <- factor(input$genelist_input_choice_venn, levels = input$genelist_input_choice_venn)
            RNAdata <- RNAdata %>% dplyr::filter(Group %in% group)
            if(is.element(TRUE, duplicated(RNAdata$Gene)) == TRUE) validate("Duplication of gene names are not allowed.")
            rownames(RNAdata) <- RNAdata$Gene
            return(RNAdata)
          }
        }
      }
    })
  })
  RNAseqDEG_norm_venn <- reactive({
    if(input$RNAseq_data_type_venn != "Result" && input$RNAseq_data_type_venn != "List" && !is.null(pre_RNAseqDEG_venn())){
      count <- RNAseq_file_venn()
      collist <- gsub("\\_.+$", "", colnames(count))
      dds <- pre_RNAseqDEG_venn()
      contrast <- c("con", unique(collist))
      normalized_counts <- counts(dds, normalized=TRUE)
      return(normalized_counts)
    }
  })
  RNAseqDEG_anno_venn <- reactive({
    if(is.null(pre_RNAseq_file_venn())) validate ("Please upload 'RNA-seq data'")
    if(input$Species_venn == "not selected") validate ("Please select 'Species'")
    return(RNAseqDEG_ann(RNAdata=RNAseqDEG_venn(),Species=input$Species_venn,gene_type = gene_type_venn_DEG_result(),input_type=input$RNAseq_data_type_venn))
  })
  gene_type_venn_DEG_result <- reactive({
    if(is.null(pre_RNAseq_file_venn())) validate ("Please upload 'RNA-seq data'")
    if(input$Species_venn == "not selected") validate ("Please select 'Species'")
    return(gene_type(my.symbols=rownames(RNAseqDEG_venn()),org=org_venn(),Species=input$Species_venn))
  })
  output$venn_RNAseq_raw <- renderDataTable({
    if(input$RNAseq_data_type_venn != "Result"){
      if(!is.null(RNAseq_file_venn())){
        RNAseq_file_venn() 
      }
    }
  })
  output$RNAseq_condition_venn <- renderText({
    if(!is.null(pre_RNAseq_file_venn())){
      paste0(input$RNAseq_cond1_venn, " vs ", input$RNAseq_cond2_venn, "\n",
             "Please confirm if the Log2FoldChange in the uploaded result file was calulated as Log2(",input$RNAseq_cond1_venn,"/",input$RNAseq_cond2_venn,").")
    }
  })
  output$RNAseq_RPmean_venn <- renderText({
    paste0("RP > 0 gene is associated with the indicated intersection within ",input$peak_distance_venn," kb from its TSS.\n",
           "RP = 0 gene is not associated with the indicated intersection within ",input$peak_distance_venn," kb from its TSS.")
  })
  output$venn_DEG_result <- renderDataTable({
    if(!is.null(pre_RNAseq_file_venn())){
      RNAseqDEG_venn() 
    }
  })
  
  RNAseq_name_venn <- reactive({
    if(!is.null(pre_RNAseq_file_venn())){
      if(input$RNAseq_data_type_venn != "Result"){
        if(!is.null(RNAseq_file_venn())){
          if(input$RNAseq_data_type_venn == "List"){
            return(c(unique(RNAseqDEG_venn()$Group)))
          }else{
            collist <- gsub("\\_.+$", "", colnames(RNAseq_file_venn()))
            cond1 <- paste0(unique(collist)[1],"_high")
            cond2 <- paste0(unique(collist)[2],"_high")
            return(c(cond1, cond2))
          }
        }
      }else{
        cond1 <- paste0(input$RNAseq_cond1_venn,"_high")
        cond2 <- paste0(input$RNAseq_cond2_venn,"_high")
        return(c(cond1, cond2))
      }
    }
  })
  pre_mmAnno_venn <- reactive({
    if(input$Species_venn == "not selected") validate ("Please select 'Species'")
    mmAnno_list <- list()
    for(name in names(venn_overlap()$peaklist)){
      peak <- venn_overlap()$peaklist[[name]]
      mAnn <- mmAnno(peak=peak,
                     genomic_region="Genome-wide",
                     txdb=txdb_venn(),
                     peak_distance=input$peak_distance_venn,
                     mode=input$RNAseq_mode_venn,
                     group_name=name,
                     distribution=vendistribution()$peaks[[name]],
                     DAR=NULL)
      mmAnno_list[[name]] <- mAnn
    }
    return(mmAnno_list)
  })
  
  
  
  
  
  
  
  
  mmAnno_venn <- reactive({
    return(pre_mmAnno_venn()[[input$intersect_RNA]])
  })
  
  RP_venn <- reactive({
    withProgress(message = "Calculating regulatory potential",{
      if(!is.null(pre_mmAnno_venn())){
        RP_list <- list()
        for(name in names(pre_mmAnno_venn())){
          data <- RP_f(mmAnno=pre_mmAnno_venn()[[name]],txdb=txdb_venn())
          RP_list[[name]] <- data
        }
        return(RP_list)
      }
      incProgress(1)
    })
  })
  regulatory_potential_venn <- reactive({
    if(!is.null(RNAseqDEG_anno_venn()) && !is.null(RP_venn()) && 
       !is.null(RNAseq_name_venn())){
      RP_list <- list()
      for(name in names(RP_venn())){
        data <- regulatory_potential_f(species=input$Species_venn,data=RNAseqDEG_anno_venn(),
                                       result_geneRP= RP_venn()[[name]],DEG_fc=input$DEG_fc_venn,
                                       DEG_fdr=input$DEG_fdr_venn,name=RNAseq_name_venn())
        RP_list[[name]] <- data
      }
      return(RP_list)
    }
  })
  
  output$ks_plot_venn <- renderPlot({    
    if(!is.null(RNAseqDEG_anno_venn()) && !is.null(venn_select_RNA_debounce()) &&
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && 
       input$Species_venn != "not selected"  && !is.null(pre_mmAnno_venn())){
      regulatory_potential_venn()[[input$intersect_RNA]]
    }
  })
  ks_tables_venn <- reactive({
    if(!is.null(RNAseqDEG_anno_venn()) && !is.null(venn_select_RNA_debounce()) &&
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && 
       input$Species_venn != "not selected"  && !is.null(pre_mmAnno_venn())){
      df <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
      for(name in names(regulatory_potential_venn())){
        df2 <- regulatory_potential_venn()[[name]]$statistics
        df2$Peak <- name
        df <- rbind(df,df2)
      }
      df <- df %>% dplyr::select(Peak,everything())
      return(df)
    }
  })
  output$ks_table_venn <- renderDataTable({    
    if(!is.null(RNAseqDEG_anno_venn()) && !is.null(venn_select_RNA_debounce()) &&
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && 
       input$Species_venn != "not selected"  && !is.null(pre_mmAnno_venn())){
      ks_tables_venn()
    }
  })
  
  output$DEG_fc_venn <- renderUI({
    numericInput("DEG_fc_venn","Fold change cutoff for RNA-seq data",
                 min=0,max=NA,value=1.5,step = 0.5)
  })
  output$DEG_fdr_venn <- renderUI({
    numericInput("DEG_fdr_venn","FDR cutoff for RNA-seq data",
                 min=0,max=1, value=0.05,step = 0.001)
  })
  
  output$RNAseqGroup_venn <- renderUI({
    if(input$Species_venn != "not selected" && 
       !is.null(regulatory_potential_venn())){
      if(!is.null(RP_all_table_venn())){
        selectInput("RNAseqGroup_venn","Group (intersection:RNAseq)",
                    unique(RP_all_table_venn()$Group),
                    multiple = FALSE)
      }
    }
  })
  
  RNAseq_boxplot_venn <- reactive({
    RNA <- RNAseqDEG_anno_venn()
    cond1 <- gsub("\\_.+$", "", RNAseq_name_venn()[1])
    cond2 <- gsub("\\_.+$", "", RNAseq_name_venn()[2])
    
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    genelist <- list()
    bed_list <- list()
    for(name in venn_select_RNA_debounce()){
      data <- merge(RNA,RP_venn()[[name]], by="gene_id",all=T)
      if(input$RNAseq_data_type_venn != "List") data <- dplyr::filter(data, baseMean != 0)
      data$group <- "Others"
      data$group[data$sumRP > 2] <- "2 < RP"
      data$group[data$sumRP <= 2 & data$sumRP > 1] <- "1 < RP < 2"
      data$group[data$sumRP <= 1 & data$sumRP > 0.1] <- "0.1 < RP < 1"
      data$group[data$sumRP <= 0.1 & data$sumRP > 0.01] <- "0.01 < RP < 0.1"
      data$group[data$sumRP <= 0.01 & data$sumRP > 0] <- "0 < RP < 0.01"
      level <- NULL
      col <- NULL
      if(dim(dplyr::filter(data, group == "Others"))[1] != 0) {
        level <- c(level,"Others")
        col <- c(col,"gray")
      }
      if(dim(dplyr::filter(data, group == "0 < RP < 0.01"))[1] != 0) {
        level <- c(level,"0 < RP < 0.01")
        col <- c(col,"blue")
      }
      if(dim(dplyr::filter(data, group == "0.01 < RP < 0.1"))[1] != 0) {
        level <- c(level,"0.01 < RP < 0.1")
        col <- c(col,"#00BFC4")
      }
      if(dim(dplyr::filter(data, group == "0.1 < RP < 1"))[1] != 0) {
        level <- c(level,"0.1 < RP < 1")
        col <- c(col,"lightgreen")
      }
      if(dim(dplyr::filter(data, group == "1 < RP < 2"))[1] != 0) {
        level <- c(level,"1 < RP < 2")
        col <- c(col,"#F8766D")
      }
      if(dim(dplyr::filter(data, group == "2 < RP"))[1] != 0) {
        level <- c(level,"2 < RP")
        col <- c(col,"red")
      }
      data$group <- factor(data$group,levels=level,ordered=TRUE)
      data$intersection <- as.factor(name)
      df <- rbind(df,data)
      if(input$RNAseq_data_type_venn != "List"){
        genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,RNAlog2FC = -1*data$log2FoldChange,
                                       sumRP = data$sumRP, gene_id = data$gene_id)
      }else{
        genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,
                                       sumRP = data$sumRP, gene_id = data$gene_id)
      }
      bed_list1 <- list()
      for(name2 in unique(data$group)){
        if(name2 != "Others"){
          table <- genelist[[name]] %>% dplyr::filter(group == name2)
          gene <- table$gene_id
          if(length(gene) != 0){
            up_peak3 <- NULL
            up_peak <- subset(pre_mmAnno_venn()[[name]], gene_id %in% gene)
            up_peak2 <- as.data.frame(up_peak)
            up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
            up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
            up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
            mcols(up_peak3) <- DataFrame(Group = paste0(name,":",name2))
            bed_list1[[name2]] <- up_peak3
          }
        }
      }
      bed_list[[name]] <- bed_list1
    }
    data <- df 
    if(input$RNAseq_data_type_venn != "List"){
      data$log10FoldChange <- log10(2^data$log2FoldChange)
      data$log10FoldChange <- as.numeric(data$log10FoldChange)
      for(name in venn_select_RNA_debounce()){
        check <- data %>% dplyr::filter(intersection == name) %>% 
          dplyr::filter(group != "Others") %>% summarise(n())
        if(check <= 1) data <- data %>% dplyr::filter(intersection != name)
      }
      if(dim(data)[1] == 0) validate("boxplot: There are few genes with |RP| > 1")
      
      data$intersection <- gsub("-","-\n",data$intersection)
      collist <- unique(data$group)
      col <-c("gray","blue","#00BFC4","lightgreen","#F8766D","red")
      if (length(collist) >= 3){
        stat.test <- data %>% dplyr::group_by(intersection) %>% 
          tukey_hsd(log10FoldChange ~ group)
        stat.test <- stat.test %>% add_significance("p.adj")
        stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      }else{
        group1 <- dplyr::filter(data, group == collist[1])
        group2 <- dplyr::filter(data, group == collist[2])
        if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
          stat.test <- data %>% dplyr::group_by(intersection) %>% 
            t_test(log10FoldChange ~ group)
          stat.test <- stat.test %>% add_significance()
          stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
        }else stat.test <- NULL
      }
      p <- try(ggpubr::ggboxplot(data, x = "group", y = "log10FoldChange",
                                 fill = "group", scales = "free", 
                                 xlab = FALSE, ylab = paste0("RNAseq log10(",cond2,"/",cond1,")"))+theme_bw(base_size = 15)+
                 xlab(NULL)+scale_fill_manual(values = col) + scale_x_discrete(labels = label_wrap_gen(8)))
      if(input$Statistics_venn !="not selected") p <- p + stat_pvalue_manual(stat.test,hide.ns = T, size = 5)         
      p <- facet(p, facet.by = "intersection",
                 panel.labs.background = list(fill = "transparent", color = "transparent"),
                 scales = "free", short.panel.labs = T, panel.labs.font = list(size=15))
      stat.test <- stat.test[,-2]
      stat.test <- stat.test[,1:9]
      colnames(stat.test)[1] <- "Peak"
      df <- list()
      df[["plot"]] <- p
      df[["statistical_test"]] <- stat.test
    }else{
      df <- list()
    }
    df[["genelist"]] <- genelist
    df[["bedlist"]] <- bed_list
    print(bed_list)
    
    return(df)
  })
  
  
  RNAseq_popu_venn <- reactive({
    if(!is.null(pre_RP_all_table_venn()) && !is.null(RNAseq_name_venn())){
      withProgress(message = "Preparing barplot",{
        up_name <- RNAseq_name_venn()[2]
        down_name <- RNAseq_name_venn()[1]
        df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        for(name in venn_select_RNA_debounce()){
          epi_name <- paste0(name, "_associated")
          table <- pre_RP_all_table_venn()[[name]] %>% dplyr::mutate(
            type =if_else(withPeakN > 0, epi_name, "not_associated")
          )
          table$Group <- gsub(".+\\:","",table$Group)
          table$Group <- paste0(table$Group," genes")
          table$type <- factor(table$type,levels = c(epi_name,"not_associated"))
          table$intersection <- name
          df <- rbind(df, table)
        }
        table <- df
        table$intersection <- gsub("-","-\n",table$intersection)
        table2 <- table %>% group_by(Group, intersection, withPeakN) %>%
          summarise(count = n(), .groups = "drop")
        p <- ggplot(table2,aes(x = withPeakN,y= count,fill = intersection)) +
          geom_col(position=position_dodge2(preserve = "single")) +
          theme_bw(base_size = 15)+facet_wrap(~Group,scales = "free",ncol = 3) +
          xlab("Number of associated peaks")+guides(fill=guide_legend(title="associated_\npeak_type"))
        incProgress(1)
      })
      return(p)
    }
  })
  output$int_box_venn <- renderPlot({
    withProgress(message = "Boxplot",{
      if(!is.null(venn_select_RNA_debounce()) && !is.null(pre_RP_all_table_venn()) && 
         input$Species_venn != "not selected" && !is.null(pre_mmAnno_venn()) &&
         !is.null(RNAseqDEG_anno_venn())){
        if(input$RNAseq_data_type_venn == "List") validate("Boxplot is not available for this input type. Please use 'DEG_result' or 'Raw_count' data.")
        RNAseq_boxplot_venn()[["plot"]]
      }
    })
  })
  output$int_boxplot_table_venn <- DT::renderDataTable({
    if(!is.null(venn_select_RNA_debounce()) && !is.null(pre_RP_all_table_venn()) && 
       input$Species_venn != "not selected" && !is.null(pre_mmAnno_venn()) &&
       !is.null(RNAseqDEG_anno_venn())){
      if(input$Statistics_venn !="not selected"){
        if(input$RNAseq_data_type_venn == "List") validate("Boxplot is not available for this input type. Please use 'DEG_result' or 'Raw_count' data.")
        RNAseq_boxplot_venn()[["statistical_test"]]
      }
    }
  })
  output$bar_rna_venn <- renderPlot({
    if(!is.null(venn_select_RNA_debounce()) && !is.null(pre_RP_all_table_venn()) && 
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && !is.null(RNAseq_popu_venn()) &&
       input$Species_venn != "not selected" && !is.null(pre_mmAnno_venn()) &&
       !is.null(RNAseqDEG_anno_venn())){
      RNAseq_popu_venn()
    }
  })
  output$bar_chip_venn <- renderPlot({
    if(!is.null(venn_select_RNA_debounce()) && !is.null(pre_RP_all_table_venn()) && 
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && !is.null(ChIPseq_popu_venn()) &&
       input$Species_venn != "not selected" && !is.null(pre_mmAnno_venn()) &&
       !is.null(RNAseqDEG_anno_venn())){
      ChIPseq_popu_venn()
    }
  })
  ChIPseq_popu_venn <- reactive({
    RNA <- RNAseqDEG_anno_venn()
    if(!is.null(RNA) && !is.null(pre_mmAnno_venn())){
      df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
      for(name in venn_select_RNA_debounce()){
        mmano <- as.data.frame(pre_mmAnno_venn()[[name]])
        mmano$locus <- paste0(mmano$seqnames,":",mmano$start,"-",mmano$end)
        up_name <- paste0(RNAseq_name_venn()[2]," gene only")
        down_name <- paste0(RNAseq_name_venn()[1]," gene only")
        merge <- merge(mmano,RNA, by = "gene_id")
        if(input$RNAseq_data_type_venn != "List"){
          merge <- merge  %>% dplyr::mutate(Up_gene = if_else(padj < input$DEG_fdr_venn & log2FoldChange > log(input$DEG_fc_venn,2), 1, 0),
                                            Down_gene = if_else(padj < input$DEG_fdr_venn & log2FoldChange < -log(input$DEG_fc_venn,2), 1, 0),
                                            NS_gene = if_else(padj >= input$DEG_fdr_venn | baseMean == 0 | is.na(padj) |
                                                                (padj <= input$DEG_fdr_venn & abs(log2FoldChange) <= log(input$DEG_fc_venn,2)), 1, 0))
          merge[is.na(merge)] <- 0
          merge2 <- merge %>% group_by(locus,Group) %>% summarise(Total_associated_gene = sum(Up_gene)+sum(Down_gene)+sum(NS_gene),Group = Group,
                                                                  Up_gene = sum(Up_gene), Down_gene = sum(Down_gene), NS_gene = sum(NS_gene))
          table <- merge2 %>% dplyr::mutate(type = if_else(Up_gene > 0 & Down_gene == 0 & NS_gene == 0, up_name,
                                                           if_else(Up_gene == 0 & Down_gene > 0 & NS_gene == 0, down_name,
                                                                   if_else(Up_gene == 0 & Down_gene == 0 & NS_gene > 0, "NS gene only", "Multiple type of genes"))))
          table$type <- factor(table$type,levels = c(down_name,up_name,"NS gene only","Multiple type of genes"))
        }else{
          merge2 <- merge %>% group_by(locus, Group.y) %>% summarise(Count = n()) 
          merge2 <- spread(merge2,Group.y,Count)
          merge2[is.na(merge2)] <- 0
          namelist <- unique(merge$Group.y)
          merge2$Total_associated_gene <- rowSums(merge2[,-1],na.rm = TRUE)
          vec <- c()
          merge3 <- merge2[,-1]
          merge3 <- merge3[,-length(colnames(merge3))]
          for(i in 1:dim(merge3)[1]){
            type <- merge3[i,]
            type <- paste(names(type[which(type != 0)]),collapse = "&")
            vec <- c(vec,type) 
          }
          merge2$type <- vec
          table <- merge2
        }
        table$Group <- paste0(table$Group," peak")
        if(max(table$Total_associated_gene) < 5) table$Total_associated_gene <- as.character(table$Total_associated_gene) else table$Total_associated_gene <- as.numeric(table$Total_associated_gene)
        table$intersection <- name
        df <- rbind(df, table)
      }
      table <- df
      table$intersection <- gsub("-","\n",table$intersection)
      
      p2 <- ggplot(table,aes(x = Total_associated_gene, fill = type)) + geom_bar(position = "stack") +
        theme_bw(base_size = 15)+facet_wrap(~intersection,scales = "free") +
        ylab("Number of peaks") +
        xlab("Number of associated genes")+guides(fill=guide_legend(title="associated_gene_type"))
      if(input$RNAseq_data_type_venn != "List"){
        col <- c("#00BFC4","#F8766D","grey","black")
        p2 <- p2 + scale_fill_manual(values = col)
      } 
      return(p2)
    }
  })
  pre_RP_all_table_venn <- reactive({
    if(!is.null(RP_venn()) && !is.null(regulatory_potential_venn())){
      table_list <- list()
      for(name in names(RP_venn())){
        target_result <- regulatory_potential_venn()[[name]]$data
        target_result$epigenome_category <- name
        table <- NULL
        if(str_detect(target_result$gene_id[1], "FBgn")){
          symbol <- id_convert(my.symbols = target_result$gene_id,Species = input$Species_venn,type = "ENSEMBL")
          symbol <- symbol %>% distinct(ENSEMBL, .keep_all = TRUE)
          symbol <- symbol$SYMBOL
        }else symbol <- target_result$Symbol
        if(!is.null(pre_mmAnno_venn()[[name]])) {
          if(input$RNAseq_data_type_venn != "List"){
            table <- data.frame(Symbol = symbol,
                                Group = paste0(target_result$epigenome_category,":",target_result$gene_category),
                                RNA_log2FC = -target_result$log2FoldChange,
                                RNA_padj = target_result$padj,
                                regulatory_potential = target_result$sumRP,
                                withPeakN = target_result$withPeakN,
                                gene_id = target_result$gene_id)
          }else{
            my.symbols <- target_result$gene_id
            gene_IDs <- id_convert(my.symbols,Species = input$Species_venn,type = "ENTREZID")
            colnames(gene_IDs) <- c("gene_id","SYMBOL")
            gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T)
            target_result <- merge(target_result, gene_IDs, by="gene_id")
            table <- data.frame(Symbol = target_result$SYMBOL,
                                Group = paste0(target_result$epigenome_category,":",target_result$Group),
                                regulatory_potential = target_result$sumRP,
                                withPeakN = target_result$withPeakN,
                                gene_id = target_result$gene_id)
          }
          table_list[[name]] <- table
        }
      }
      return(table_list)
    }
  })
  RP_all_table_venn <- reactive({
    if(!is.null(pre_RP_all_table_venn())){
      df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
      table <- pre_RP_all_table_venn()
      for(name in names(table)){
        df <- rbind(df, table[[name]])
      }
      return(df)
    }
  })
  
  RP_selected_table_venn <- reactive({
    table <- RP_all_table_venn() %>% dplyr::filter(Group == input$RNAseqGroup_venn)
    return(table)
  })
  
  output$RP_table_venn <- renderDT({
    if(!is.null(input$RNAseqGroup_venn) && !is.null(venn_select_RNA_debounce()) &&
       !is.null(input$peak_distance_venn && !is.null(pre_mmAnno_venn())) &&
       !is.null(regulatory_potential_venn())){
      RP_selected_table_venn() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  output$download_vennintbar = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        gridExtra::grid.arrange(RNAseq_popu_venn(), ChIPseq_popu_venn(), ncol = 1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_vennintbox = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), "%Y%m%d_"),"boxplot",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- pdf_h(input$venn_select_RNA)
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- pdf_w(input$venn_select_RNA)
        }else pdf_width <- input$venn_pdf_width
        fs <- c()
        dir.create("boxplot/",showWarnings = FALSE)
        if(input$RNAseq_data_type_venn != "List"){
          boxplot_table <- paste0("boxplot/","boxplot.txt")
          boxplot <- paste0("boxplot/","boxplot.pdf")
          
          fs <- c(fs,boxplot_table,boxplot)
          pdf(boxplot, height = pdf_height, width = pdf_width)
          print(RNAseq_boxplot_venn()[["plot"]])
          dev.off()
          write.table(RNAseq_boxplot_venn()[["statistical_test"]],boxplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
        }
        dir.create("boxplot/gene_list/",showWarnings = FALSE)
        dir.create("boxplot/bed/",showWarnings = FALSE)
        for(name in venn_select_RNA_debounce()){
          genelist <- paste0("boxplot/gene_list/",name,".txt")
          up <- RNAseq_boxplot_venn()[["genelist"]][[name]]
          fs <- c(fs,genelist)
          write.table(up,genelist,col.names = T,row.names = F,sep = "\t",quote = F)
        }
        
        for(name in names(RNAseq_boxplot_venn()[["bedlist"]])){
          dir.create(paste0("boxplot/bed/",name),showWarnings = FALSE)
          for(file in names(RNAseq_boxplot_venn()[["bedlist"]][[name]])){
            bed <- RNAseq_boxplot_venn()[["bedlist"]][[name]][[file]]
            bed_name <- paste0("boxplot/bed/",name,"/",file,".bed")
            bed_name <- gsub(" < ","....",bed_name)
            print(bed_name)
            fs <- c(fs, bed_name)
            write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
          }
        }
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    },
    contentType = "application/zip"
  )
  
  
  output$download_vennKSplot = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), "%Y%m%d_"),"KSplot",".zip")
    },
    content = function(fname) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        fs <- c()
        dir.create("KSplot/",showWarnings = FALSE)
        for(name in names(regulatory_potential_venn())){
          ksplot <- paste0("KSplot/",name,".pdf")
          fs <- c(fs,ksplot)
          pdf(ksplot, height = pdf_height, width = pdf_width)
          print(regulatory_potential_venn()[[name]])
          dev.off()
        }
        kstable <- paste0("KSplot/ks_test.txt")
        fs <- c(fs,kstable)
        write.table(ks_tables_venn(), kstable, row.names = F, sep = "\t", quote = F)
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    }
  )
  output$download_RP_venn_table = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), "%Y%m%d_"),"RP_table",".zip")
    },
    content = function(fname) {
      withProgress(message = "Preparing download",{
        fs <- c()
        dir.create("RP_table/",showWarnings = FALSE)
        dir.create("RP_table/selected_table(epigenome--RNA)/",showWarnings = FALSE)
        dir.create("RP_table/selected_bed(epigenome--RNA)/",showWarnings = FALSE)
        RP_summary <- paste0("RP_table/summary.txt")
        fs <- c(fs,RP_summary)
        write.table(RP_all_table_venn(), RP_summary, row.names = F, sep = "\t", quote = F)
        for(name in unique(RP_all_table_venn()$Group)){
          RP_selected <- paste0("RP_table/selected_table(epigenome--RNA)/",name,".txt")
          RP_selected_bed <- paste0("RP_table/selected_bed(epigenome--RNA)/",name,".bed")
          RP_selected <- gsub(":","--",RP_selected)
          RP_selected_bed <- gsub(":","--",RP_selected_bed)
          fs <- c(fs,RP_selected,RP_selected_bed)
          table <- RP_all_table_venn() %>% dplyr::filter(Group == name)
          write.table(table, RP_selected, row.names = F, sep = "\t", quote = F)
          
          gene <- table$gene_id
          peak <- gsub("\\:.+$", "", name)
          if(!is.null(pre_mmAnno_venn())) {
            up_peak <- subset(pre_mmAnno_venn()[[peak]], gene_id %in% gene)
            up_peak2 <- as.data.frame(up_peak)
            up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
            up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
            up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
            mcols(up_peak3) <- DataFrame(Group = "up")
            y <- as.data.frame(up_peak3)
          }
          write.table(y, RP_selected_bed, row.names = F, col.names = F,sep = "\t", quote = F)
        }
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    }
  )
  
  observeEvent(input$venn_DEG_result, ({
    updateCollapse(session,id =  "input_collapse_venn_RP", open="KS_panel")
  }))
  #Venn Integrative trackplot-------
  observeEvent(input$RP_table_venn_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_venn_RP", open="int_Trackplot_panel")
  }))
  output$int_igv_uprange_venn <- renderUI({
    if(!is.null(int_goi_gene_position_venn()) && !is.null(int_goi_promoter_position_venn())){
      y <- int_goi_promoter_position_venn()
      gene_position <- int_goi_gene_position_venn()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("int_igv_uprange_venn","Range (x-axis):",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim_venn <- renderUI({
    numericInput("int_igv_ylim_venn","Max peak intensity (y-axis):", value = 2, min = 0)
  })
  
  int_goi_promoter_position_venn<- reactive({
    if(!is.null(input$RP_table_venn_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table_venn()[input$RP_table_venn_rows_selected,]$gene_id
      y <- NULL
      if(!is.null(pre_mmAnno_venn())) {
        up_peak <- subset(mmAnno_venn(), gene_id %in% gene)
        up_peak2 <- as.data.frame(up_peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      return(y)
    }
  })
  
  int_goi_gene_position_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected)){
      gene <- RP_selected_table_venn()[input$RP_table_venn_rows_selected,]$gene_id
      gene_position <- as.data.frame(subset(gene_position_venn(), gene_id %in% gene))
      return(gene_position)
    }
  })
  
  
  int_data_track_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected)){
      return(data_trac(y=int_goi_promoter_position_venn(),gene_position=int_goi_gene_position_venn(),
                       gen=ref_venn(),txdb=txdb_venn(),org=org_venn(),
                       bw_files=bws_venn(),
                       track_additional_files=NULL))
    }
  })
  
  int_highlight_trackplot_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected) && !is.null(input$int_igv_uprange_venn)){
      library(Gviz)
      y <- int_goi_promoter_position_venn()
      gene_position <- int_goi_gene_position_venn()
      chr <- gene_position$seqnames
      df <- int_data_track_venn()
      start <-c()
      width <- c()
      col <- c()
      fill <- c()
      for(i in 1:length(rownames(y))){
        if(y[i,]$start < input$int_igv_uprange_venn[2] && y[i,]$end > input$int_igv_uprange_venn[1]){
          start <- c(start, y[i,]$start)
          width <- c(width, y[i,]$width)
          if(y[i,]$Group == "up") {
            col <- c(col,"#E41A1C")
            fill <- c(fill,"#FBB4AE")
          }
          if(y[i,]$Group == "down") {
            col <- c(col,"#377EB8")
            fill <- c(fill,"#B3CDE3")
          }
        }
      }
      if(length(start) != 0){
        ht <- HighlightTrack(trackList = df,
                             start = start, width = width,col = col,fill = fill,
                             chromosome = chr, alpha=0.5) 
      }else ht <- NULL
      return(ht)
    }
  })
  int_goi_trackplot_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected) &&
       !is.null(int_goi_promoter_position_venn()) && 
       !is.null(int_goi_gene_position_venn()) && 
       !is.null(gtrack_venn()) &&
       !is.null(input$int_igv_uprange_venn)){
      library(Gviz)
      if(!is.null(int_highlight_trackplot_venn())){
        plot<- plotTracks(list(gtrack_venn(), int_highlight_trackplot_venn()),
                          from = input$int_igv_uprange_venn[1], 
                          to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                          type="hist")
      }else{
        df <- int_data_track_venn()
        df[["gtrack"]] <- gtrack_venn()
        plot<- plotTracks(df,
                          from = input$int_igv_uprange_venn[1], 
                          to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                          type="hist")
      }
      return(plot)
    }
  })
  output$int_trackplot_goi_venn <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$RP_table_venn_rows_selected)){
        int_goi_trackplot_venn()
      }
    })
  })
  output$download_venn_int_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(int_highlight_trackplot_venn())){
          plot<- plotTracks(list(gtrack_venn(), int_highlight_trackplot_venn()),
                            from = input$int_igv_uprange_venn[1], 
                            to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                            type="hist")
        }else{
          df <- int_data_track_venn()
          df[["gtrack"]] <- gtrack_venn()
          plot<- plotTracks(df,
                            from = input$int_igv_uprange_venn[1], 
                            to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                            type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  #Venn integrative functional enrichment-------
  int_Hallmark_set_venn <- reactive({
    if(!is.null(input$intGeneset_venn)){
      return(GeneList_for_enrichment(Species = input$Species_venn, Gene_set = input$intGeneset_venn, org = org_venn()))
    }
  })
  
  
  order_for_intGroup_venn <- reactive({
    order <- c("0 < RP < 0.01","0.01 < RP < 0.1","0.1 < RP < 1","1 < RP < 2","2 < RP")
    order <- factor(order,levels = order,ordered = TRUE)
    return(order)
  })
  output$intGroup_for_RPstatus_venn <- renderUI({
    if(input$withRNAseq_venn_enrich_type == "boxplot"){
      selectInput("intGroup_for_RPstatus_venn","Peak",names(RNAseq_boxplot_venn()[["genelist"]]),multiple = F)
    }
  })
  
  output$intGroup_venn <- renderUI({
    if(!is.null(RP_all_table_venn())){
      if(input$withRNAseq_venn_enrich_type == "boxplot"){
        if(input$RNAseq_data_type_venn == "List") validate("When using gene lists as input for RNA-seq data, RP-status-based enrichment analysis cannot be applied. \nPlease utilize the 'Relationship (Epigenome:RNAseq)' mode instead.")
        selectInput("intGroup_venn","Group (intersection:RNAseq)",
                    order_for_intGroup_venn(),selected=order_for_intGroup_venn(),multiple = T)
      }else{
        selectInput("intGroup_venn","Group (intersection:RNAseq)",
                    unique(RP_all_table_venn()$Group),multiple = T)
      }
    }
  })
  pre_intGroup_venn_debounce <- reactive({
    if(!is.null(input$intGroup_venn)) return(input$intGroup_venn) else return(NULL)
  })
  intGroup_venn_debounce <- debounce(pre_intGroup_venn_debounce, 1000)
  
  output$intGeneset_venn <- renderUI({
    selectInput('intGeneset_venn', 'Gene Set', gene_set_list)
  })
  
  withRNAseq_enrichment_analysit_genelist_venn <- reactive({
    if(input$withRNAseq_venn_enrich_type == "boxplot"){
      df <- RNAseq_boxplot_venn()[["genelist"]][[input$intGroup_for_RPstatus_venn]]
      colnames(df)[2] <- "Group"
      df <- df %>% dplyr::filter(Group != "Others")
      return(df)
    }
  })
  
  selected_int_group_venn <- reactive({
    group <- intGroup_venn_debounce()
    df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
    colnames(df) <- c("ENTREZID","Group")
    for(name in group){
      
      if(input$withRNAseq_venn_enrich_type == "boxplot"){
        table <- withRNAseq_enrichment_analysit_genelist_venn()
      }else{
        table <- RP_all_table_venn()
      }
      table <- table %>% dplyr::filter(Group == name)
      if(dim(table)[1] != 0){
        entrezid <- table$gene_id
        if(str_detect(table$gene_id[1], "FBgn")){
          my.symbols <- gsub("\\..*","", table$gene_id)
          gene_IDs<-AnnotationDbi::select(org(input$Species_venn),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","ENTREZID"))
          colnames(gene_IDs) <- c("gene_id","ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(ENTREZID, .keep_all = T)
          gene_IDs <- na.omit(gene_IDs)
          table <- merge(table,gene_IDs,by="gene_id")
          entrezid <- table$ENTREZID
        }
        df2 <- data.frame(ENTREZID = entrezid, Group = table$Group)
        df <- rbind(df,df2)
      }
    }
    return(df)
  })
  int_enrich_venn <- reactive({
    return(enrich_viewer_forMulti2(data3 = selected_int_group_venn(), Species = input$Species_venn, org = org_venn(),
                                   H_t2g = int_Hallmark_set_venn(),Gene_set = input$intGeneset_venn))
  })
  int_enrich_list_venn <- reactive({
    return(enrich_gene_list(data = selected_int_group_venn(),
                            Gene_set = input$intGeneset_venn, org = org_venn(), H_t2g = int_Hallmark_set_venn()))
  })
  int_enrich_plot_venn <- reactive({
    return(enrich_genelist(data = selected_int_group_venn(),
                           enrich_gene_list = int_enrich_list_venn(),type ="withRNAseq"))
  })
  output$int_enrichment1_venn <- renderPlot({
    dotplot_for_output(data = int_enrich_venn(),
                       plot_genelist = int_enrich_plot_venn(), Gene_set = input$intGeneset_venn, 
                       Species = input$Species_venn)
  })
  int_enrich_table_venn <- reactive({
    return(enrich_for_table(data = as.data.frame(int_enrich_venn()), H_t2g = int_Hallmark_set_venn(), Gene_set = input$intGeneset_venn))
  })
  output$int_enrichment_result_venn <- DT::renderDataTable({
    int_enrich_table_venn()
  })
  output$download_venn_int_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$intGeneset_venn,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = int_enrich_venn(),
                           plot_genelist = int_enrich_plot_venn(), Gene_set = input$intGeneset_venn, 
                           Species = input$Species_venn)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_venn_int_enrichment_venn_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$intGeneset_venn,".txt")
    },
    content = function(file){write.table(int_enrich_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  #heatmap--------
  ##heatmap----------
  Venn_peak_call_files_locus <- reactive({
    data <- Venn_peak_call_files_geneId()
    df <- list()
    for(name in names(data)){
      data2 <- data[[name]]
      data2_gr <- with(data2, GRanges(seqnames = seqnames, 
                                      ranges = IRanges(start,end),
                                      strand = strand))
      names(data2_gr) <- paste0(data2$ENTREZID,"_",data2$locus)
      df[[name]] <- data2_gr
    }
    return(df)
  })
  
  integrated_legend_venn <- reactive({
    lgd <- lgd(files2 = Venn_peak_call_files_locus())
    return(lgd)
  })
  
  integrated_heatlist_venn <- reactive({
    if(input$integrated_heatmapButton_venn > 0 && updateCounter_int_venn$i > 0){
      ht_list <- NULL
      if(!is.null(integrated_heatmap_add1_venn())) ht_list <- ht_list + integrated_heatmap_add1_venn()[["heatmap"]]
      if(!is.null(integrated_heatmap_add2_venn())) ht_list <- ht_list + integrated_heatmap_add2_venn()[["heatmap"]]
      if(!is.null(integrated_heatmap_add3_venn())) ht_list <- ht_list + integrated_heatmap_add3_venn()[["heatmap"]]
      if(!is.null(integrated_heatmap_add4_venn())) ht_list <- ht_list + integrated_heatmap_add4_venn()[["heatmap"]]
      if(!is.null(rnaseq_DEGs2_venn())) ht_list <- ht_list + rnaseq_DEGs_heatmap_venn()
      if(!is.null(rnaseq_count2_venn())) ht_list <- ht_list + rnaseq_count_heatmap_venn()
      return(ht_list)
    }else return(NULL)
  })
  updateCounter_int_venn <- reactiveValues(i = 0)
  
  observe({
    input$integrated_heatmapButton_venn
    isolate({
      updateCounter_int_venn$i <- updateCounter_int_venn$i + 1
    })
  })
  
  
  #Restart
  defaultvalues_venn <- observeEvent(integrated_heatlist_venn(), {
    isolate(updateCounter_int_venn$i == 0)
    updateCounter_int_venn <<- reactiveValues(i = 0)
  }) 
  integrated_heatlist_venn_plot <- reactive({
    if(!is.null(input$venn_heatmap_group)){
      if(input$integrated_heatmapButton_venn > 0 && !is.null(Venn_peak_call_files_locus()) &&
         !is.null(integrated_heatlist_venn())){
        return(draw(integrated_heatlist_venn(),annotation_legend_list = list(integrated_legend_venn()),
                    heatmap_legend_side = "bottom", ht_gap = unit(2, "mm")))
      }
    }
  })
  output$integrated_heatmap_venn <- renderPlot({
    if(!is.null(input$venn_heatmap_group)){
      if(input$integrated_heatmapButton_venn > 0 && !is.null(Venn_peak_call_files_locus()) &&
         !is.null(integrated_heatlist_venn())){
        integrated_heatlist_venn_plot()
      }
    }
  })
  output$download_integrated_heatmap_venn = downloadHandler(
    filename = "Integrated_heatmap.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(integrated_heatlist_venn_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$rnaseq_count_venn <- renderUI({
    if(input$Species_venn != "not selected"){
      fileInput("pair_rnaseq_count_venn",
                "Select RNA-seq normalized count files",
                accept = c("txt","csv","xlsx"),
                multiple = TRUE,
                width = "90%")
    }
  })
  output$rnaseq_DEGs_venn <- renderUI({
    if(input$Species_venn != "not selected"){
      fileInput("pair_rnaseq_DEGs_venn",
                "Select RNA-seq DEG result files",
                accept = c("txt","csv","xlsx"),
                multiple = TRUE,
                width = "90%")
    }
  })
  rnaseq_count_venn <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      tmp <- input$pair_rnaseq_count_venn
      upload = list()
      if(is.null(input$pair_rnaseq_count_venn) && input$goButton_venn > 0 )  {
        tmp = "data/RNAseq_count.txt"
        upload[["rna"]]<- read_df(tmp)
        return(upload)
      }else if(is.null(tmp)) {
        return(NULL)
      }else{
        return(read_dfs(tmp))
      }
    })
  })
  rnaseq_DEGs_venn <- reactive({
    withProgress(message = "Importing DEG result files, please wait",{
      tmp <- input$pair_rnaseq_DEGs_venn
      upload = list()
      if(is.null(input$pair_rnaseq_DEGs_venn) && input$goButton_venn > 0 )  {
        tmp = "data/RNAseq.txt"
        upload[["rna"]]<- read_df(tmp)
        return(upload)
      }else if(is.null(tmp)) {
        return(NULL)
      }else{
        return(read_dfs(tmp))
      }
    })
  })
  output$pre_zscoring_venn <- renderUI({
    if(length(names(rnaseq_count_venn())) > 1){
      selectInput("pre_zscoring_venn","pre_zscoring for multiple normalized count data",choices = c(TRUE,FALSE),selected = TRUE)
    }
  })
  gene_type_rnaseq_DEGs_venn <- reactive({
    files <- rnaseq_DEGs_venn()
    return(gene_type_for_integrated_heatmap(files=files,Species=input$Species_venn,org=org_venn()))
  })
  rnaseq_DEGs2_venn <- reactive({
    files <- rnaseq_DEGs_venn()
    return(rnaseqDEGs_for_integrated_heatmap(files = files,Species=input$Species_venn,gene_type = gene_type_rnaseq_DEGs_venn()))
  })
  gene_type_rnaseq_counts_venn <- reactive({
    files <- rnaseq_count_venn()
    return(gene_type_for_integrated_heatmap(files=files,Species=input$Species_venn,org=org_venn()))
  })
  rnaseq_count2_venn <- reactive({
    files <- rnaseq_count_venn()
    return(rnaseqCounts_for_integrated_heatmap(files = files,Species=input$Species_venn,gene_type = gene_type_rnaseq_counts_venn(),pre_zscoring = input$pre_zscoring_venn))
  })
  observeEvent(input$pair_rnaseq_DEGs_venn, ({
    updateCollapse(session,id =  "z-score_count_venn", open="Uploaded_DEGs_venn")
  }))
  observeEvent(input$pair_rnaseq_count_venn, ({
    updateCollapse(session,id =  "z-score_count_venn", open="z-score_multiple_count_venn_panel")
  }))
  output$rnaseq_count_output_venn <- renderDataTable({
    if(input$Species_venn != "not selected" && !is.null(rnaseq_count_venn())){
      rnaseq_count2_venn()
    }
  })
  output$rnaseq_DEGs_output_venn <- renderDataTable({
    if(input$Species_venn != "not selected" && !is.null(rnaseq_DEGs_venn())){
      rnaseq_DEGs2_venn()
    }
  })
  output$venn_heatmap_group <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      if(!is.null(venn_overlap())){
        selectInput("venn_heatmap_group", "Select intersects", choices = c(names(venn_overlap()$peaklist)),multiple = TRUE)
      }
    }
  })
  selected_grange_venn_list_for_heatmap <- reactive({
    if(!is.null(input$venn_heatmap_group)){
      Glist <- list()
      for(name in input$venn_heatmap_group){
        data <- venn_overlap()$peaklist[[name]]
        Glist[[name]] <- data
      }
      return(Glist)
    }
  })
  Venn_peak_call_files_geneId <- reactive({
    peaks <- selected_grange_venn_list_for_heatmap()
    df <- list()
    for(name in names(peaks)){
      data <- peaks[[name]]
      data2 <- as.data.frame(ChIPseeker::as.GRanges(ChIPseeker::annotatePeak(peak = data, TxDb = txdb_venn())))
      data2 <- data2 %>% dplyr::filter(!is.na(geneId))
      data2$locus <- paste0(data2$seqnames,":",data2$start,"-",data2$end)
      colnames(data2)[which(colnames(data2) == "geneId")] <- "ENTREZID"
      df[[name]] <- data2
    }
    return(df)
  })
  rnaseq_count_heatmap_venn <- reactive({
    rna <-  rnaseq_count2_venn()
    peaks <- Venn_peak_call_files_geneId()
    rna$ENTREZID <- rownames(rna)
    m_z <- data.frame(matrix(rep(NA, length(colnames(rna))), nrow=1))[numeric(0), ]
    for(name in names(peaks)){
      data <- peaks[[name]]
      data_m <- merge(rna, data, by="ENTREZID",all=T)
      data_m <- data_m %>% dplyr::filter(!is.na(locus))
      rownames(data_m) <- paste0(data_m$ENTREZID,"_",data_m$locus)
      data_m <- data_m[,2:length(colnames(rna))]
      data_m[is.na(data_m)] <- 0
      m_z <- rbind(m_z,data_m)
    }
    cond <- gsub(".+\\.", "", colnames(m_z))
    cond <- gsub("\\_.+$", "", cond)
    cond <- factor(cond, levels = unique(cond), ordered = TRUE)
    cond_color <- rainbow_hcl(length(unique(cond)),c=100)
    names(cond_color) <- unique(cond)
    if(length(names(rnaseq_count_venn())) == 1){
      file_name <- NULL
      file_name_color <- NULL
    }else{
      file_name <- gsub("\\..+$", "", colnames(m_z))
      file_name <- factor(file_name, levels = unique(file_name), ordered = TRUE)
      file_name_color <- rainbow_hcl(length(file_name))
      names(file_name_color) <- file_name
    }
    mat <- integrated_heatmap_add1_venn()[["mat"]]
    withProgress(message = "Heatmap of RNA-seq count data",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),],name = "RNA-seq\nz-score", 
                    top_annotation = HeatmapAnnotation(files = file_name, condition = cond,
                                                       col = list(files = file_name_color,
                                                                  condition = cond_color)),
                    show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "cm"),use_raster = TRUE)
      return(ht)
    })
  })
  rnaseq_DEGs_heatmap_venn <- reactive({
    rna <-  rnaseq_DEGs2_venn()
    peaks <- Venn_peak_call_files_geneId()
    rna$ENTREZID <- rownames(rna)
    m_z <- data.frame(matrix(rep(NA, length(colnames(rna))), nrow=1))[numeric(0), ]
    for(name in names(peaks)){
      data <- peaks[[name]]
      data_m <- merge(rna, data, by="ENTREZID",all=T)
      data_m <- data_m %>% dplyr::filter(!is.na(locus))
      rownames(data_m) <- paste0(data_m$ENTREZID,"_",data_m$locus)
      data_m <- data_m[,1:length(colnames(rna))]
      data_m[is.na(data_m)] <- 0
      m_z <- rbind(m_z,data_m)
    }
    for(i in 1:length(colnames(rna))){
      m_z[,i] <- -1 * as.numeric(m_z[,i])
    }
    m_z[m_z > 5]<- 5
    m_z[m_z < -5]<- -5
    colnames(m_z) <- gsub("\\.log2F.+$", "", colnames(m_z))
    mat <- integrated_heatmap_add1_venn()[["mat"]]
    withProgress(message = "Heatmap of RNA-seq log2FoldChange",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),2:length(colnames(rna))],name = "RNA-seq\nlog2FC", 
                    show_row_names = FALSE, width = unit(2.5, "cm"), column_names_gp = grid::gpar(fontsize = 9),
                    use_raster = TRUE,column_names_side = "top",show_column_dend = FALSE,
                    col = c("blue","white","gold"))
      return(ht)
    })
  })
  
  output$integrated_bw2_venn <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_2_venn",
                  "Option: Select additional bigwig files (blue)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_2_venn",
                "Option: Select additional bigwig files (blue)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw3_venn <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_3_venn",
                  "Option: Select additional bigwig files (green)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_3_venn",
                "Option: Select additional bigwig files (green)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw4_venn <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_4_venn",
                  "Option: Select additional bigwig files (purple)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_4_venn",
                "Option: Select additional bigwig files (purple)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  
  pre_integrated_additional2_venn <-reactive({
    return(pre_integrated_additional(input$integrated_bw_2_venn))
  })
  integrated_additional2_venn <- reactive({
    return(bws_ordering(bws=pre_integrated_additional2_venn(),sample_order=input$sample_order_venn_comb2))
  })
  pre_integrated_additional3_venn <-reactive({
    return(pre_integrated_additional(input$integrated_bw_3_venn))
  })
  integrated_additional3_venn <- reactive({
    return(bws_ordering(bws=pre_integrated_additional3_venn(),sample_order=input$sample_order_venn_comb3))
  })
  pre_integrated_additional4_venn <-reactive({
    return(pre_integrated_additional(input$integrated_bw_4_venn))
  })
  integrated_additional4_venn <- reactive({
    return(bws_ordering(bws=pre_integrated_additional4_venn(),sample_order=input$sample_order_venn_comb4))
  })
  integrated_heatmap_add1_venn <- reactive({
    if(!is.null(bws_venn())){
      h <- batch_heatmap(files2 = Venn_peak_call_files_locus(),files_bw = bws_venn(),
                         color = c("white","red"),signal = "red")
      return(h)
    }
  })
  integrated_heatmap_add2_venn <- reactive({
    if(!is.null(integrated_additional2_venn())){
      h <- batch_heatmap(files2 = Venn_peak_call_files_locus(),files_bw = integrated_additional2_venn(),
                         color = c("white","darkblue"),signal = "darkblue")
      return(h)
    }
  })
  integrated_heatmap_add3_venn <- reactive({
    if(!is.null(integrated_additional3_venn())){
      h <- batch_heatmap(files2 = Venn_peak_call_files_locus(),files_bw = integrated_additional3_venn(),
                         color = c("white","darkgreen"),signal = "green")
      return(h)
    }
  })
  integrated_heatmap_add4_venn <- reactive({
    if(!is.null(integrated_additional4_venn())){
      h <- batch_heatmap(files2 = Venn_peak_call_files_locus(),files_bw = integrated_additional4_venn(),
                         color = c("white","purple"),signal = "purple")
      return(h)
    }
  })
  
  #Clustering--------
  output$Spe_clustering <- renderText({
    if(input$Species_clustering == "not selected" && input$Genomic_region_clustering == "Promoter") print("Please select 'Species'")
  })
  output$Spe_clustering_track <- renderText({
    if(input$Species_clustering == "not selected") print("Please select 'Species'")
  })
  output$Spe_int_clustering <- renderText({
    if(input$Species_clustering == "not selected") print("Please select 'Species'")
  })
  observeEvent(input$goButton_clustering,({
    updateSelectInput(session,inputId = "Species_clustering","Species_clustering",species_list, selected = "Homo sapiens (hg19)")
  }))
  output$sample_order_kmeans_pattern <- renderUI({
    if(!is.null(pre_kmeans_additional())){
      selectInput(inputId = "sample_order_kmeans_pattern","Sample order:",
                  choices =  names(pre_kmeans_additional()),
                  selected = names(pre_kmeans_additional()),multiple = TRUE)
    }
  })
  output$sample_order_kmeans_track <- renderUI({
    if(!is.null(pre_track_additional_files_clustering())){
      selectInput(inputId = "sample_order_kmeans_track","Sample order:",
                  choices =  names(pre_track_additional_files_clustering()),
                  selected = names(pre_track_additional_files_clustering()),multiple = TRUE)
    }
  })
  output$sample_order_kmeans_track_int <- renderUI({
    if(!is.null(pre_int_track_additional_files_clustering())){
      selectInput(inputId = "sample_order_kmeans_track_int","Sample order:",
                  choices =  names(pre_int_track_additional_files_clustering()),
                  selected = names(pre_int_track_additional_files_clustering()),multiple = TRUE)
    }
  })
  
  ##------
  output$file1_clustering <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("file1_clustering",label=NULL,choices = bw,multiple=T)
    }else{
      fileInput("file1_clustering",NULL,
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$peak_call_file1_clustering <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bed <- c(str_subset(list,".narrowPeak"),str_subset(list,".bed"))
      selectInput("peak_call_file1_clustering",label=NULL,choices = bed,multiple=T)
    }else{
      fileInput("peak_call_file1_clustering",NULL,
                accept = c("bed","narrowPeak"),
                multiple = TRUE,
                width = "80%")
    }
  })
  
  updateCounter_createCount_clustering <- reactiveValues(i = 0)
  
  observe({
    input$createcountButton_clustering
    isolate({
      updateCounter_createCount_clustering$i <- updateCounter_createCount_clustering$i + 1
    })
  })
  
  
  #Restart
  observeEvent(bws_clustering(), {
    isolate(updateCounter_createCount_clustering$i == 0)
    updateCounter_createCount_clustering <<- reactiveValues(i = 0)
  }) 
  observeEvent(peak_call_files_clustering(), {
    isolate(updateCounter_createCount_clustering$i == 0)
    updateCounter_createCount_clustering <<- reactiveValues(i = 0)
  }) 
  # input data ------------------------------------------------------------------------------
  org1_clustering <- reactive({
    return(org(Species = input$Species_clustering))
  })
  txdb_clustering <- reactive({
    if(input$Species_clustering != "not selected"){
      return(txdb_function(Species = input$Species_clustering))
    }
  })
  promoter_region_clustering <- reactive({
    if(input$Genomic_region_clustering == "Promoter"){
      if(input$Species_clustering != "not selected"){
        return(promoter_clustering(txdb_clustering(),upstream = input$upstream_clustering, downstream = input$downstream_clustering,
                                   filter = input$pair_filter_cluster))
      }
    }else return(promoter_clustering(upstream = input$upstream_clustering, downstream = input$downstream_clustering,
                                     input_type = "Genome-wide",files =peak_call_files_clustering(),
                                     filter = input$pair_filter_cluster))
  })
  gene_position_clustering <- reactive({
    if(input$Species_clustering != "not selected"){
      return(genes(txdb_clustering()))
    }
  })
  
  bws_clustering <- reactive({
    if(input$data_file_type_clustering == "Row1" || input$data_file_type_clustering == "Row1_count"){
      if(is.null(input$file1_clustering)){
        if(input$goButton_clustering > 0 ){
          df<-list()
          df[["A_1"]] <- "data/bigwig/A_1.BigWig"
          df[["A_2"]] <- "data/bigwig/A_2.BigWig"
          df[["B_1"]] <- "data/bigwig/B_1.BigWig"
          df[["B_2"]] <- "data/bigwig/B_2.BigWig"
          return(df)
        }
        return(NULL)
      }else{
        if(length(list.files("./Volume/")) > 0){
          files <- input$file1_clustering
          names(files) <- bigwig_name(input$file1_clustering)
        }else{
          files<-c()
          name<-c()
          for(nr in 1:length(input$file1_clustering[, 1])){
            file <- input$file1_clustering[[nr, 'datapath']]
            name <- c(name, bigwig_name(input$file1_clustering[nr,]$name))
            files <- c(files,file)
          }
          names(files)<-name
        }
        return(files)
      }
    }
  })
  bws_count_clustering <- reactive({
    tmp <- input$file1_count_clustering$datapath
    if(is.null(input$file1_count_clustering) && input$goButton_clustering > 0 )  tmp = "data/bws_count.txt"
    if(!is.null(tmp)) return(read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = ""))
  })
  bws_order_clustering<-reactive({
    if(input$data_file_type_clustering == "Row1" || input$data_file_type_clustering == "Row1_count"){
      return(bws_ordering(bws=bws_clustering(),sample_order=input$sample_order_clustering,additional=FALSE))
    }
  })
  bam_clustering <- reactive({
    if(input$data_file_type_clustering == "Row2"){
      if(is.null(input$file_bam_clustering)){
        if(input$goButton_clustering > 0 ){
          df<-list()
          df[["A_1"]] <- "data/bigwig/A_1.BigWig"
          df[["A_2"]] <- "data/bigwig/A_2.BigWig"
          df[["B_1"]] <- "data/bigwig/B_1.BigWig"
          df[["B_2"]] <- "data/bigwig/B_2.BigWig"
          return(df)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file_bam_clustering[, 1])){
          file <- input$file_bam_clustering[[nr, 'datapath']]
          name <- c(name, bigwig_name(input$file_bam_clustering[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
        return(files)
      }
    }
  })
  
  peak_call_files_clustering <- reactive({
    if(input$Genomic_region_clustering == "Genome-wide"){
      if(is.null(input$peak_call_file1_clustering)){
        if(input$goButton_clustering > 0 ){
          df <- list()
          df[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
          df[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
          df[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
          df[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
          name <- c("Ctrl_1","Ctrl_2","Sen_1","Sen_2")
          files2 <- lapply(df, GetGRanges, simple = TRUE)
          names(files2)<-name
          return(files2)
        }
        return(NULL)
      }else{
        if(length(list.files("./Volume/")) > 0){
          files <- input$peak_call_file1_clustering
          files2 <- lapply(files, GetGRanges, simple = TRUE)
          name <- gsub(".+\\/","",input$peak_call_file1_clustering)
          names(files2) <- gsub("\\..+$", "", name)
        }else{
          files<-c()
          name<-c()
          for(nr in 1:length(input$peak_call_file1_clustering[, 1])){
            file <- input$peak_call_file1_clustering[[nr, 'datapath']]
            name <- c(name, gsub("\\..+$", "", input$peak_call_file1_clustering[nr,]$name))
            files <- c(files,file)
          }
          files2 <- lapply(files, GetGRanges, simple = TRUE)
          names(files2)<-name
        }
        return(files2)
      }
    }
  })
  
  bw_files_clustering <- reactive({
    return(BigWigFileList(bws_clustering()))
  })
  
  output$input_bw_files_clustering <- DT::renderDT({
    if(input$data_file_type_clustering == "Row1" || input$data_file_type_clustering == "Row1_count"){
      uploaded_files = names(bws_clustering())
      data.frame(uploaded_files = uploaded_files)
    }else{
      uploaded_files = names(bam_clustering())
      as.data.frame(uploaded_files) 
    }
  })
  
  output$input_peak_call_files_clustering <- DT::renderDataTable({
    if(input$Genomic_region_clustering == "Genome-wide"){
      uploaded_files = names(peak_call_files_clustering())
      as.data.frame(uploaded_files)
    }
  })
  genelist_clustering <- reactive({
    withProgress(message = "Importing a gene list file, please wait",{
      tmp <- input$genelist_file1_clustering$datapath
      if(is.null(input$genelist_file1_clustering) && input$goButton_clustering > 0 )  tmp = "data/RNAseq.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") {
          df2 <- readxl::read_xlsx(tmp) 
          df2 <- as.data.frame(df2)
          df <- try(data.frame(row.names = df2[,1]),silent = T)
          if(class(df) != "try-error") {
            if(dim(df2)[2] == 2){
              df <- data.frame(row.names = df2[,1],a = df2[,2])
              colnames(df)[1] <- colnames(df2)[2]
            }else{
              rownames(df2) <- df2[,1]
              df <- df2[,-1]
              colnames(df) <- gsub("-",".",colnames(df))
            }
          }
        }
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
        rownames(df) = gsub("\"", "", rownames(df))
        if(length(colnames(df)) != 0){
          if(str_detect(colnames(df)[1], "^X\\.")){
            colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
          }
        }
        return(df)
      }
    })
  })
  bw_count_clustering2 <- reactive({
    if(input$data_file_type_clustering == "Row1" && !is.null(bw_files_clustering())){
      if(input$Genomic_region_clustering == "Promoter"){
        if(input$Species_clustering != "not selected"){
          count <-Bigwig2count(bw = bw_files_clustering(),promoter_region_clustering(),
                               Species = input$Species_clustering,input_type =input$Genomic_region_clustering)
          return(count) 
        }
      }else{
        if(!is.null(peak_call_files_clustering())){
          return(Bigwig2count(bw = bw_files_clustering(),promoter_region_clustering(),
                              input_type =input$Genomic_region_clustering))
        }
      }
    }
    if(input$Species_clustering != "not selected" && input$data_file_type_clustering == "Row2" && 
       !is.null(bam_clustering())){
      withProgress(message = "converting bam to gene count data",{
        switch(input$Pair_or_single_clustering,
               "Paired-end" = seq <- TRUE,
               "Single-end" = seq <- FALSE)
        regionsToCount <- data.frame(GeneID = paste("ID", seqnames(promoter_region_clustering()), 
                                                    start(promoter_region_clustering()), end(promoter_region_clustering()), sep = "_"), 
                                     Chr = seqnames(promoter_region_clustering()), 
                                     Start = start(promoter_region_clustering()), End = end(promoter_region_clustering()), 
                                     Strand = strand(promoter_region_clustering()))
        a <- as.data.frame(promoter_region_clustering())
        Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
        if(input$Genomic_region_clustering == "Promoter"){
          count <- featureCounts(bam_clustering(),annot.ext = regionsToCount,isPairedEnd = seq,
                                 countMultiMappingReads =F,maxFragLength = 100)$counts
          rownames(count) <- Row.name
          colnames(count) <- names(bam_clustering())
          return(count)
        }else{
          if(!is.null(peak_call_files_clustering())){
            count <- featureCounts(bam_clustering(),annot.ext = regionsToCount,isPairedEnd = seq,
                                   countMultiMappingReads =F,maxFragLength = 100)$counts
            rownames(count) <- Row.name
            colnames(count) <- names(bam_clustering())
            return(count)
          }
        }
      })
    }
  })
  pre_bw_count_clustering <- reactive({
    if(input$data_file_type_clustering == "Row1" && !is.null(bw_files_clustering()) && updateCounter_createCount_clustering$i > 0){
      if(input$Genomic_region_clustering == "Promoter"){
        if(input$Species_clustering != "not selected"){
          count <-bw_count_clustering2()
          if(!is.null(genelist_clustering())){
            colnum <- length(colnames(count))
            count <- merge(count,genelist_clustering(),by=0)
            rownames(count) <- count$Row.names
            count <- count[,2:(1+colnum)]
          }
          return(count) 
        }
      }else{
        return(bw_count_clustering2())
      }
    }
    if(input$Species_clustering != "not selected" && input$data_file_type_clustering == "Row2" && 
       !is.null(bam_clustering())){
      return(bw_count_clustering2())
    }
    if(input$data_file_type_clustering == "Row1_count"){
      if(!is.null(bws_count_clustering())) return(bws_count_clustering())
    }
  })
  bw_count_clustering <- reactive({
    count <- pre_bw_count_clustering()
    order <- input$sample_order_clustering
    if(!is.null(count)){
      if(dim(count)[2] != 0){
        data <- count[,order]
        return(data)
      }
    }
  })
  observeEvent(pre_bw_count_clustering(),({
    updateSelectizeInput(session,inputId = "sample_order_clustering","Sample order:",
                         choices = colnames(pre_bw_count_clustering()),selected = colnames(pre_bw_count_clustering()))
  }))
  output$raw_count_table_clustering <- DT::renderDataTable({
    bw_count_clustering()
  })
  
  output$download_raw_count_clustering_table <- downloadHandler(
    filename = function() {
      if (input$Genomic_region_clustering=='Promoter'){
        paste0("bigwig2count-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("bigwig2count-genomeWide",".txt")
      }},
    content = function(file){write.table(bw_count_clustering(), file, row.names = T, sep = "\t", quote = F)}
  )
  observeEvent(bws_clustering(), ({
    updateCollapse(session,id =  "input_collapse_panel_clustering", open="raw_count_panel")
  }))
  observeEvent(bam_clustering(), ({
    updateCollapse(session,id =  "input_collapse_panel_clustering", open="raw_count_panel")
  }))
  observeEvent(peak_call_files_clustering(), ({
    if(input$data_file_type_clustering != "Row1_count"){
      updateCollapse(session,id =  "input_collapse_panel_clustering", open="peak_call_files_panel")
    }
  }))
  observeEvent(bw_count_clustering(), ({
    if(input$data_file_type_clustering == "Row1_count"){
      updateCollapse(session,id =  "input_collapse_panel_clustering", open="raw_count_panel")
    }
  }))
  
  # clustering limma--------
  limma_clustering <- reactive({
    if(!is.null(input$regression_mode_clustering)){
      count <- bw_count_clustering()
      collist <- gsub("\\_.+$", "", colnames(count))
      collist <- gsub(" ", ".", collist)
      colnames(count) <- gsub("\\-","\\_",colnames(count))
      count <- log(count + 1,2)
      collist <- gsub("\\-","\\_",collist)
      if(is.element(TRUE, duplicated(collist)) == TRUE){
        withProgress(message = "Detecting differential accessible region",{
          meta <- data.frame(condition = factor(collist))
          design <- model.matrix(~0+collist)
          colnames(design) <- factor(unique(collist),levels = unique(collist))
          cont <- c()
          for(i in 1:choose(n=length(unique(meta$condition)),k=2)){
            contrast = paste0(as.character(unique(meta$condition)[combn(x=length(unique(meta$condition)),m=2)[1,i]]),"-",as.character(unique(meta$condition)[combn(x=length(unique(meta$condition)),m=2)[2,i]]))
            cont <- c(cont,contrast)
          }
          cont.matrix <- makeContrasts(contrasts=cont, levels=design)
          eset = new("ExpressionSet", exprs=as.matrix(count))
          fit1 <- lmFit(eset,design)
          fit2 <- contrasts.fit(fit1, cont.matrix)
          fit3 <- eBayes(fit2,trend = TRUE ,robust = input$regression_mode_clustering)
          result <- topTable(fit3,coef=1:length(cont), number = 1e12)
          lab <- paste0("log2(",cont,")")
          lab <- gsub("-","/",lab)
          lab <- gsub("_","-",lab)
          if(length(cont) != 1) label <- c(lab,"AveExpr","F","p_value","padj") else label <- c(lab,"AveExpr","F","p_value","padj","B")
          colnames(result) <- label
          incProgress(1)
        })
        return(result)
      }
    }
  })
  output$clustering_Statistical_data <- renderDataTable({
    if(!is.null(bw_count_clustering())){
      if(!is.null(limma_clustering())){
        limma_clustering()
      }
    }
  })
  output$download_clustering_Statistical_table = downloadHandler(
    filename ="PCA_table.txt",
    content = function(file){write.table(limma_clustering(), 
                                         file, row.names = T,col.names = NA, sep = "\t", quote = F)}
  )
  output$clustering_stastics_notion <- renderText({
    if(is.null(limma_clustering()) && !is.null(bw_count_clustering())){
      "Statistical analysis is skipped due to the absence of replicates (n = 1 for each sample).\nYou can use the fold change cut-off alone for filtering."
    }
  })
  
  output$regression_mode_clustering <- renderUI({
    count <- bw_count_clustering()
    collist <- gsub("\\_.+$", "", colnames(count))
    if(is.element(TRUE, duplicated(collist)) == TRUE){
      radioButtons("regression_mode_clustering","Regression",
                   c('least squares'=FALSE,
                     'robust'=TRUE
                   ), selected = FALSE)
    }
  })
  output$filtered_fdr <- renderText({
    if(is.null(limma_clustering())){
      return(NULL)
    }else{ 
      print(paste0("The number of genomic regions after the filtration (fdr < ",input$fdr_clustering ,"): ", dim(dplyr::filter(limma_clustering(),padj < input$fdr_clustering))[1]))
    }
  })
  # clustering PCA ------------------------------------------------------------------------------
  output$clustering_pca_error <- renderText({
    p <- length(colnames(bw_count_clustering()))
    if(p < 3){
      print("PCA: The sample number must be > 2")
    }else return(NULL)
  })
  output$download_clustering_PCA = downloadHandler(
    filename = "PCA-MDS-dendrogram.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(PCAplot(data = bw_count_clustering(),plot=TRUE,legend=input$PCA_legend_clustering))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$PCA_clustering <- renderPlot({
    withProgress(message = "preparing result table",{
      if(is.null(bw_count_clustering())){
        return(NULL)
      }else{
        p <- length(colnames(bw_count_clustering()))
        if(p > 2){
          print(PCAplot(data = bw_count_clustering(),plot=TRUE,legend=input$PCA_legend_clustering))
        }
      }
    })
  })
  
  output$clustering_PCA_data <- DT::renderDataTable({
    PCAplot(data = bw_count_clustering(),plot=FALSE)
  })
  
  output$download_clustering_PCA_table = downloadHandler(
    filename ="PCA_table.txt",
    content = function(file){write.table(PCAplot(data = bw_count_clustering(),plot=FALSE), 
                                         file, row.names = T, sep = "\t", quote = F)}
  )
  #Clustering umap-------
  output$clustering_umap_n <- renderUI({
    sliderInput("clustering_n_neighbors", "n_neighbors", min = 2,
                max=100, step = 1,
                value = 15)
  })
  
  clustering_umap_plot <- reactive({
    data <- bw_count_clustering()
    if(is.null(input$clustering_n_neighbors)){
      return(NULL)
    }else{
      if(is.null(data)){
        return(NULL)
      }else{
        p<- try(umap_plot(data = data, n_neighbors = input$clustering_n_neighbors))
        return(p)
      }
    }
  })
  
  output$clustering_umap_error <- renderText({
    p <- clustering_umap_plot()
    if(length(p) == 1){
      if (class(p) == "try-error") {
        print("umap: number of neighbors must be smaller than number of items")
      }
    }else return(NULL)
  })
  
  
  output$clustering_umap <- renderPlot({
    p <- clustering_umap_plot()
    withProgress(message = "umap",{
      if(length(p) == 1){
        if (class(p) == "try-error") {
          return(NULL)
        }
      }else{
        print(p)
      }
    })
  })
  
  output$download_clustering_umap = downloadHandler(
    filename = function(){
      paste0(input$clustering_n_neighbors,"_neighbors_umap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 4.7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(clustering_umap_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  #Clustering correlation plot-----
  corrplot <-reactive({
    data <- bw_count_clustering()
    mcor <- cor(data)
    p <- ggcorrplot(corr = mcor, hc.order = TRUE, method = "square",
                    colors = c("#4b61ba", "white", "red"), lab = TRUE)
    p <- p + scale_fill_gradient2(limit = c(0,1), low = "#4b61ba", high =  "red", mid = "white", midpoint = 0.5)
    return(p)
  })
  output$correlationplot <- renderPlot({
    if(!is.null(bw_count_clustering())){
      corrplot()
    }
  })
  
  output$download_clustering_corrplot = downloadHandler(
    filename = function(){
      paste0("Corrplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(corrplot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #Clustering kmeans--------
  updateCounter_kmeans <- reactiveValues(i = 0)
  
  observe({
    input$kmeans_start
    isolate({
      updateCounter_kmeans$i <- updateCounter_kmeans$i + 1
    })
  })
  
  
  #Restart
  observeEvent(input$clustering_kmeans_number, {
    isolate(updateCounter_kmeans$i == 0)
    updateCounter_kmeans <<- reactiveValues(i = 0)
  }) 
  
  observeEvent(input$basemean_clustering, {
    isolate(updateCounter_kmeans$i == 0)
    updateCounter_kmeans <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$fc_clustering, {
    isolate(updateCounter_kmeans$i == 0)
    updateCounter_kmeans <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$fdr_clustering, {
    isolate(updateCounter_kmeans$i == 0)
    updateCounter_kmeans <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$selectFC, {
    isolate(updateCounter_kmeans$i == 0)
    updateCounter_kmeans <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$selectFC2, {
    isolate(updateCounter_kmeans$i == 0)
    updateCounter_kmeans <<- reactiveValues(i = 0)
  }) 
  output$selectFC <- renderUI({
    if(is.null(bw_count_clustering())){
      return(NULL)
    }else{
      selectizeInput("selectFC", "Option: a pair for fold change cut-off", c(unique(unique(gsub("\\_.*","", colnames(bw_count_clustering()))))),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  output$selectFC2 <- renderUI({
    if(is.null(bw_count_clustering())){
      return(NULL)
    }else{
      if(length(unique(unique(gsub("\\_.*","", colnames(bw_count_clustering()))))) > 2){
        selectizeInput("selectFC2", "Option: a pair for fold change cut-off", c(unique(unique(gsub("\\_.*","", colnames(bw_count_clustering()))))),
                       selected = "", multiple = TRUE, 
                       options = list(maxItems = 2))
      }
    }
  })
  output$kmeans_order <- renderUI({
    order <- pre_kmeans_order()
    withProgress(message = "Draw heatmap",{
      selectInput("kmeans_order","Order of clusters on heatmap",order,
                  selected = order,multiple = T)
    })
  })
  
  bw_count_clustering_anno <- reactive({
    if(!is.null(bw_count_clustering()) && !is.null(txdb_clustering())){
      withProgress(message = "preparing annotation",{
        data <- peak_kmeans_grange()
        data2 <- ChIPseeker::annotatePeak(data, TxDb= txdb_clustering())
        data <- as.data.frame(ChIPseeker::as.GRanges(data2))
        Row.name <- paste0(data$seqnames,":",data$start,"-",data$end)
        data$locus <- Row.name
        my.symbols <- data$geneId
        gene_IDs<-id_convert(my.symbols,input$Species_clustering,type="ENTREZID")
        colnames(gene_IDs) <- c("geneId","NearestGene")
        data <- merge(data, gene_IDs, by="geneId")
        data <- data %>% distinct(locus, .keep_all = T)
        rownames(data)<-data$locus
        return(data)
      })
    }
  })
  
  bw_count_clusterin_anno2 <- reactive({
    if(input$Genomic_region_clustering == "Genome-wide"){
      anno <- bw_count_clustering_anno()
      anno <- dplyr::arrange(anno, locus)
      data <- data.frame(NearestGene = anno$NearestGene,
                         locus = anno$locus)
      clustering_kmeans_pattern_extract <- clustering_kmeans_pattern_extract()
      clustering_kmeans_pattern_extract$locus <- rownames(clustering_kmeans_pattern_extract())
      clustering_kmeans_pattern_extract <- dplyr::arrange(clustering_kmeans_pattern_extract, locus)
      data2 <- merge(data,clustering_kmeans_pattern_extract,by="locus")
      rownames(data2) <- data2$locus
      data2 <- data2[, - which(colnames(data2) == "locus")]
      return(data2)
    }
  })
  
  output$clustering_kmeans_num <- renderUI({
    if(is.null(bw_count_clustering())){
      return(NULL)
    }else{
      withProgress(message = "Preparing kmeans clustering",{
        sliderInput("clustering_kmeans_number", "k-means number", min = 1,
                    max=20, step = 1,
                    value = 2)
      })
    }
  })
  
  bw_count_clustering_fc_basemean_cutoff <- reactive({
    data <- as.data.frame(bw_count_clustering_fc_basemean_cutoff2())
    if(is.null(data)){
      return(NULL)
    }else{
      if(is.null(input$selectFC2)) return(data)
      if(length(input$selectFC2) != 2) return(data)
      result <- limma_clustering()
      if(!is.null(result)){
        result_padj_cutoff <- result %>% dplyr::filter(padj < input$fdr_clustering)
        result2 <- merge(result_padj_cutoff,data,by=0)
        if(dim(data)[1] != 0){
          cond1 <- input$selectFC2[1]
          cond2 <- input$selectFC2[2]
          log2fc_1 <- paste0("log2(",cond1,"/",cond2,")")
          log2fc_2 <- paste0("log2(",cond2,"/",cond1,")")
          Log2FoldChange <- try(dplyr::select(result2, .data[[log2fc_1]]))
          if(class(Log2FoldChange) == "try-error") Log2FoldChange <- try(dplyr::select(data, .data[[log2fc_2]]))
          data$Log2FoldChange <- Log2FoldChange
          data3 <- data %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc_clustering))
          data3 <- data3[, - which(colnames(data3) == "Log2FoldChange")]
        }else data3 <- NULL
        return(data3)
      }else{
        if(dim(data)[1] != 0){
          cond1 <- input$selectFC2[1]
          cond2 <- input$selectFC2[2]
          cond1_ave <- data %>% dplyr::select(starts_with(cond1))
          cond2_ave <- data %>% dplyr::select(starts_with(cond2))
          Log2FoldChange <- log((cond1_ave + 0.01)/(cond2_ave + 0.01),2)
          data$Log2FoldChange <- Log2FoldChange
          data2 <- data %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc_clustering))
          data2 <- data2[, - which(colnames(data2) == "Log2FoldChange")]
        }else data2 <- NULL
        return(data2)
      }
    }
  })
  bw_count_clustering_fc_basemean_cutoff2 <- reactive({
    data <- as.data.frame(bw_count_clustering())
    if(is.null(data)){
      return(NULL)
    }else{
      result <- limma_clustering()
      if(!is.null(result)){
        result_padj_cutoff <- result %>% dplyr::filter(padj < input$fdr_clustering)
        data <- merge(result_padj_cutoff,data,by=0)
        rownames(data) <- data$Row.names
        data2 <- data[,-1:-(length(colnames(result))+1)]
        if(length(input$selectFC) != 2) return(data2)
        data2 <- data2 %>% dplyr::filter(apply(.,1,mean) > input$basemean_clustering)
        if(dim(data2)[1] != 0){
          cond1 <- input$selectFC[1]
          cond2 <- input$selectFC[2]
          log2fc_1 <- paste0("log2(",cond1,"/",cond2,")")
          log2fc_2 <- paste0("log2(",cond2,"/",cond1,")")
          Log2FoldChange <- try(dplyr::select(data, .data[[log2fc_1]]))
          if(class(Log2FoldChange) == "try-error") Log2FoldChange <- try(dplyr::select(data, .data[[log2fc_2]]))
          data2$Log2FoldChange <- Log2FoldChange
          data3 <- data2 %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc_clustering))
          data3 <- data3[, - which(colnames(data3) == "Log2FoldChange")]
        }else data3 <- NULL
        return(data3)
      }else{
        if(length(input$selectFC) != 2) return(data)
        if(dim(data)[1] != 0){
          cond1 <- input$selectFC[1]
          cond2 <- input$selectFC[2]
          cond1_ave <- data %>% dplyr::select(starts_with(cond1))
          cond2_ave <- data %>% dplyr::select(starts_with(cond2))
          Log2FoldChange <- log((cond1_ave + 0.01)/(cond2_ave + 0.01),2)
          data$Log2FoldChange <- Log2FoldChange
          data2 <- data %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc_clustering))
          data2 <- data2[, - which(colnames(data2) == "Log2FoldChange")]
        }else data2 <- NULL
        return(data2)
      }
    }
  })
  output$filtered_region <- renderText({
    if(is.null(bw_count_clustering_fc_basemean_cutoff2())){
      return(NULL)
    }else{ 
      if(length(input$selectFC2) != 2){
        if(length(input$selectFC) != 2){
          
        }else print(paste0("The number of genomic regions after the filtration (fdr < ",input$fdr_clustering,", |log2(", input$selectFC[1],"/", input$selectFC[2],")| > ", log2(input$fc_clustering),"): ", length(rownames(bw_count_clustering_fc_basemean_cutoff()))))
      }else{
        if(length(input$selectFC) != 2){
          print(paste0("The number of genomic regions after the filtration (fdr < ",input$fdr_clustering,", |log2(", input$selectFC2[1],"/", input$selectFC2[2],")| > ", log2(input$fc_clustering),"): ", length(rownames(bw_count_clustering_fc_basemean_cutoff()))))
        }else print(paste0("The number of genomic regions after the filtration (fdr < ",input$fdr_clustering ,
                           ", |log2(", input$selectFC[1],"/", input$selectFC[2],")| > ", log2(input$fc_clustering),", |log2(", input$selectFC2[1],"/", input$selectFC2[2],")| > ", log2(input$fc_clustering),"): ", length(rownames(bw_count_clustering_fc_basemean_cutoff()))))
      }
    }
  })
  
  
  bw_count_clustering_cutoff <- reactive({
    data <- bw_count_clustering_fc_basemean_cutoff()
    if(is.null(data)){
      return(NULL)
    }else{
      data2 <- data[order(apply(data,1,mad), decreasing = T),]
      return(data2)
    }
  })
  
  clustering_data_z <- reactive({
    data <- bw_count_clustering_cutoff()
    if(is.null(data)){
      return(NULL)
    }else{
      data.z <- genescale(data, axis = 1, method = "Z")
      data.z <- na.omit(data.z)
      return(data.z)
    }
  })
  pre_clustering_kmeans <- reactive({
    data.z <- clustering_data_z()
    if(is.null(data.z) || 
       input$kmeans_start == 0 || updateCounter_kmeans$i == 0){
      return(NULL)
    }else{
      withProgress(message = "k-means clustering",{
        set.seed(123)
        cl = consensus_kmeans(data.z, input$clustering_kmeans_number, 100)
        names(cl) <- rownames(data.z)
        incProgress(1)
      })
      return(cl)
    }
  })
  pre_kmeans_order <- reactive({
    data.z <- clustering_data_z()
    cl = pre_clustering_kmeans()
    if(is.null(cl) || length(unique(cl)) == 1  || 
       input$kmeans_start == 0 || updateCounter_kmeans$i == 0){
      return(NULL)
    }else{
      cl <- data.frame(cl)
      colnames(cl)[1] <- "cluster" 
      data2 <- merge(cl,data.z, by=0)
      rownames(data2)<-data2[,1]
      data2 <- data2[,-1]
      df <- data.frame(matrix(rep(NA, 1), nrow=1))[numeric(0), ]
      for(i in 1:input$clustering_kmeans_number){
        data3 <- data2 %>% dplyr::filter(cluster == i)
        data4 <- apply(data3[,-1],2,sum)
        df <- rbind(df,data4)
      }
      colnames(df) <- colnames(data2[,-1])
      order <- hclust(dist(df), "average")$order
      return(order)
    }
  })
  clustering_kmeans <- reactive({
    data.z <- clustering_data_z()
    if(is.null(data.z) || is.null(pre_clustering_kmeans()) || is.null(input$kmeans_order) ||
       input$kmeans_start == 0 || updateCounter_kmeans$i == 0){
      return(NULL)
    }else{
      if(length(input$kmeans_order) == length(unique(pre_clustering_kmeans()))){
        withProgress(message = "k-means clustering",{
          ht <- Heatmap(data.z, name = "z-score",
                        column_order = colnames(data.z),
                        clustering_method_columns = 'ward.D2',show_row_dend = FALSE,
                        cluster_row_slices = F, split = factor(pre_clustering_kmeans(),levels = input$kmeans_order),
                        show_row_names = F,column_names_side = "top",use_raster = TRUE)
          ht <- draw(ht)
          return(ht)
        })
      }
    }
  })
  
  clustering_kmeans_cluster <- reactive({
    ht <- clustering_kmeans()
    data.z <- clustering_data_z()
    data <- bw_count_clustering_cutoff()
    if(is.null(ht) || is.null(data.z)){
      return(NULL)
    }else{
      r.dend <- row_dend(ht)
      rcl.list <- row_order(ht)
      lapply(rcl.list, function(x) length(x))
      Cluster <- NULL
      if(!is.null(input$clustering_kmeans_number)){
        if(length(lapply(rcl.list, function(x) length(x))) != input$clustering_kmeans_number){
          return(NULL)
        }else{
          out <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
          for (i in input$kmeans_order){ 
            clu <- t(t(row.names(data.z[suppressWarnings(row_order(ht)[[i]]),])))
            clu <- cbind(clu, paste("cluster", i, sep=""))
            out <- rbind(out, clu)
          }
          colnames(out) <- c("GeneID", "Cluster")
          out <- as.data.frame(out)
          rownames(out) <- out$GeneID
          clusterCount <- merge(out, data, by=0)
          rownames(clusterCount) <- clusterCount$GeneID
          clusterCount <- clusterCount[,-1:-2]
          return(clusterCount)
        }
      }else return(NULL)
    }
  })
  output$clustering_select_kmean <- renderUI({
    withProgress(message = "preparing kmeans clustering",{
      clusters <- clustering_kmeans_cluster()
      if(is.null(clusters) || input$kmeans_start == 0){
        return(NULL)
      }else{
        selectInput("clustering_select_kmean", "Create custom cluster (combine multiple clusters)", choices = c(unique(clusters$Cluster)),multiple = T)
      } 
    })
  })
  
  clustering_kmeans_pattern_extract <- reactive({
    clusters <- clustering_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$clustering_select_kmean)){
        return(NULL)
      }else{
        cluster_name <- input$clustering_select_kmean
        clusterCount <- dplyr::filter(clusters, Cluster %in% cluster_name)
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  
  output$download_clustering_all_bed = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"all_clusters_bed",".zip")
    },
    content = function(fname){
      fs <- c()
      dir.create("all_clusters_bed/",showWarnings = FALSE)
      for(name in unique(clustering_kmeans_cluster()$Cluster)){
        file <- paste0("all_clusters_bed/",name,".bed")
        fs <- c(fs,file)
        clusterCount <- dplyr::filter(clustering_kmeans_cluster(), Cluster == name)
        if(input$Genomic_region_clustering == "Genome-wide"){
          kmeans <- range_changer(clusterCount)
          kmeans <- with(kmeans, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
        }else{
          kmeans <- symbol2gene_id(clusterCount,org1_clustering()) %>% distinct(gene_id, .keep_all = T)
          kmeans <- subset(promoter_region_clustering(), gene_id %in% kmeans$gene_id) 
        }
        write.table(apply(as.data.frame(kmeans),2,as.character),file,col.names = F,row.names = F,sep = "\t",quote = F)
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  output$clustering_kmeans_extract_table <- renderDT({
    if(!is.null(input$clustering_select_kmean)  && input$kmeans_start > 0){
      if(input$Genomic_region_clustering == "Genome-wide"){
        if(input$Species_clustering == "not selected"){
          clustering_kmeans_pattern_extract() %>%
            datatable(
              selection = "single")
        }else{
          if(!is.null(peak_kmeans_grange()) ){
            bw_count_clusterin_anno2() %>%
              datatable(
                selection = "single")
          }
        } 
      }else{
        clustering_kmeans_pattern_extract() %>%
          datatable(
            selection = "single")
      }
    }
  })
  
  output$clustering_kmeans_heatmap <- renderPlot({
    withProgress(message = "plot heatmap",{
      ht <- clustering_kmeans()
      if(is.null(ht) || input$kmeans_start == 0){
        return(NULL)
      }else{
        print(ht)
      }
    })
  })
  
  output$download_clustering_kmeans_extract_count_bed = downloadHandler(
    filename = function() {
      paste0(paste(input$clustering_select_kmean,collapse = "+"),".bed")
    },
    content = function(file){
      write.table(as.data.frame(peak_kmeans_grange()), file, row.names = F, col.names = F,sep = "\t", quote = F)
    }
  )
  output$download_clustering_kmeans_heatmap = downloadHandler(
    filename = function() {
      paste(input$clustering_kmeans_number,"kmeans_heatmap.pdf",sep = "_")
    },
    content = function(file){
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(clustering_kmeans())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_clustering_kmeans_extract_count = downloadHandler(
    filename = function() {
      paste0("kmeans_selected_table",".txt")
    },
    content = function(file){
      if(input$Species_clustering == "not selected"){
        write.table(clustering_kmeans_pattern_extract(), file, row.names = T, sep = "\t", quote = F)
      }else write.table(bw_count_clusterin_anno2(), file, row.names = T, sep = "\t", quote = F)
    }
  )
  #clustering kmeans peak--------
  output$peak_pattern_kmeans_additional <- renderUI({
    fileInput("peak_pattern_kmeans_add",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  pre_kmeans_additional <-reactive({
    if(!is.null(input$peak_pattern_kmeans_add)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$peak_pattern_kmeans_add[, 1])){
        file <- input$peak_pattern_kmeans_add[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$peak_pattern_kmeans_add[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(bigwig_breakline(files))
    }
  })
  kmeans_additional <- reactive({
    if(!is.null(input$peak_pattern_kmeans_add)){
      if(!is.null(input$sample_order_kmeans_pattern)) return(pre_kmeans_additional()[input$sample_order_kmeans_pattern]) 
    }
  })
  output$peak_pattern_kmeans_heat_range <- renderUI({
    if(!is.null(peak_kmeans_grange())){
      withProgress(message = "Preparing peak pattern",{
        rg <- kmeans_pattern_range()
        sliderInput("peak_kmeans_range","Intensity range",value=rg,min = 0,max=ceiling(rg*2),step=ceiling(rg*2)/100)
      })
    }
  })
  kmeans_pattern_range <- reactive({
    rg <- c()
    sig <- peak_pattern_function(grange=peak_kmeans_grange(), files=bws_order_clustering(),
                                 additional=kmeans_additional(),plot = FALSE)
    for(name in names(sig)){
      rg <- c(rg, mean(sig[[name]][,50]))
    }
    rg <- max(rg) + max(rg)*0.1
    return(rg)
  })
  peak_kmeans_grange <- reactive({
    if(!is.null(clustering_kmeans_pattern_extract())){
      if(input$Genomic_region_clustering == "Genome-wide"){
        kmeans <- range_changer(clustering_kmeans_pattern_extract())
        kmeans <- with(kmeans, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
      }else{
        kmeans <- symbol2gene_id(clustering_kmeans_pattern_extract(),org1_clustering()) %>% distinct(gene_id, .keep_all = T)
        kmeans <- subset(promoter_region_clustering(), gene_id %in% kmeans$gene_id) 
      }
      return(kmeans)
    }else return(NULL)
  })
  
  peak_kmeans_alinedHeatmap <- reactive({
    if(!is.null(input$peak_kmeans_range) && input$peak_kmeans_range > 0){
      bigwig <- bigwig_breakline(bws_order_clustering())
      heatmap <- peak_pattern_function(grange=peak_kmeans_grange(), files=bigwig,
                                       additional=kmeans_additional(),rg = input$peak_kmeans_range)
      return(heatmap)
    }
  })
  output$peak_pattern_kmeans_heatmap <- renderPlot({
    withProgress(message = "feature aligned heatmap",{
      if(!is.null(clustering_kmeans_pattern_extract())){
        plot_grid(peak_kmeans_alinedHeatmap()[["heat"]])
      }
    })
  })
  output$peak_pattern_kmeans_line <- renderPlot({
    withProgress(message = "feature aligned distribution",{
      if(!is.null(clustering_kmeans_pattern_extract()) && !is.null(peak_kmeans_alinedHeatmap())){
        matplot(peak_kmeans_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_kmeans_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
      }
    })
  })
  
  output$download_peak_pattern_kmeans_heatmap = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_kmeans_alinedHeatmap()[["heat"]]))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_kmeans_line = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_alignedDistibution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        matplot(peak_kmeans_alinedHeatmap()[["line"]],upstream=2000, downstream=2000,
                type="l",ylab="density",lty=1,xaxt="n")
        axis(1,at = c(0,25,50,75,100),labels = c(-2000,-1000,0,1000,2000))
        legend("topright", legend=colnames(peak_kmeans_alinedHeatmap()[["line"]]), col=1:6,
               lty=1, lwd=1)
        dev.off()
        incProgress(1)
      })
    }
  )
  #Clustering trackplot-------
  ref_clustering <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species_clustering))
    return(ref)
  })
  
  
  observeEvent(input$clustering_kmeans_extract_table_rows_selected, ({
    updateCollapse(session,id =  "clustering_kmeans_collapse_panel", open="kmeans_track_panel")
  }))
  observeEvent(input$clustering_kmeans_extract_table_rows_selected,({
    if(input$Species_clustering != "not selected"){
      if(!is.null(goi_gene_position_clustering()) && !is.null(goi_promoter_position_clustering())){
        y <- goi_promoter_position_clustering()
        gene_position <- goi_gene_position_clustering()
        start_position <- min(c(y$start,gene_position$start))
        end_position <- max(c(y$end,gene_position$end))
        updateSliderInput(session,"igv_uprange_clustering","Range (x-axis):",
                          value = c(start_position,end_position),
                          step = 100, min = start_position - 10000, max = end_position + 10000)
      }
    }
  }))
  
  output$igv_uprange_clustering <- renderUI({
    if(input$Species_clustering != "not selected"){
      if(!is.null(goi_gene_position_clustering()) && !is.null(goi_promoter_position_clustering())){
        y <- goi_promoter_position_clustering()
        gene_position <- goi_gene_position_clustering()
        start_position <- min(c(y$start,gene_position$start))-1000
        end_position <- max(c(y$end,gene_position$end))+1000
        sliderInput("igv_uprange_clustering","Range (x-axis):",value = c(start_position,end_position),
                    step = 100, min = start_position - 10000, max = end_position + 10000)
      }
    }
  })
  output$igv_ylim_clustering <- renderUI({
    numericInput("igv_ylim_clustering","Max peak intensity (y-axis):", value = 2, min = 0)
  })
  
  gtrack_clustering <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=1)
      return(gtrack)
    })
  })
  
  goi_promoter_position_clustering<- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      library(Gviz)
      if(input$Genomic_region_clustering == "Promoter"){
        label_data <- rownames(clustering_kmeans_pattern_extract()[input$clustering_kmeans_extract_table_rows_selected,])
        gene_IDs<- id_convert(label_data,input$Species_clustering,type="SYMBOL_single")
        y <- as.data.frame(subset(promoter_region_clustering(), gene_id %in% gene_IDs))
      }else{
        data <- range_changer(bw_count_clusterin_anno2()[input$clustering_kmeans_extract_table_rows_selected,])
        y <- as.data.frame(with(data, GRanges(seqnames = chr, 
                                              ranges = IRanges(start,end))))
      }
      return(y)
    }
  })
  
  goi_gene_position_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      if(input$Genomic_region_clustering == "Promoter"){
        label_data <- rownames(clustering_kmeans_pattern_extract()[input$clustering_kmeans_extract_table_rows_selected,])
      }else{
        label_data <- bw_count_clusterin_anno2()[input$clustering_kmeans_extract_table_rows_selected,]$NearestGene
      }
      gene_IDs<- id_convert(label_data,input$Species_clustering,type="SYMBOL_single")
      gene_position <- as.data.frame(subset(gene_position_clustering(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  output$trackplot_additional_clustering <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("trackplot_additional1_clustering",
                  "Select additional bigwig files",
                  choices = bw,multiple=T)
    }else{
      fileInput("trackplot_additional1_clustering",
                "Select additional bigwig files",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  pre_track_additional_files_clustering <-reactive({
    if(!is.null(input$trackplot_additional1_clustering)){
      if(length(list.files("./Volume/")) > 0){
        files <- input$trackplot_additional1_clustering
        name <- gsub(".+\\/","",input$trackplot_additional1_clustering)
        names(files) <- gsub("\\..+$", "", name)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$trackplot_additional1_clustering[, 1])){
          file <- input$trackplot_additional1_clustering[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$trackplot_additional1_clustering[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
      }
      return(files)
    }
  })
  track_additional_files_clustering <- reactive({
    if(!is.null(input$trackplot_additional1_clustering)){
      if(input$sample_order_kmeans_track) return(bws_ordering(bws=pre_track_additional_files_clustering(),sample_order=input$sample_order_kmeans_track))
    }
  })
  data_track_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      return(data_trac(y=goi_promoter_position_clustering(),gene_position=goi_gene_position_clustering(),
                       gen=ref_clustering(),txdb=txdb_clustering(),org=org1_clustering(),filetype=input$data_file_type_clustering,
                       bw_files=bws_order_clustering(),bam_files=bam_clustering(),
                       track_additional_files=track_additional_files_clustering()))
    }
  })
  
  highlight_trackplot_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      library(Gviz)
      y <- goi_promoter_position_clustering()
      gene_position <- goi_gene_position_clustering()
      chr <- gene_position$seqnames
      df <- data_track_clustering()
      if(y$start < input$igv_uprange_clustering[2] && y$end > input$igv_uprange_clustering[1]){
        withProgress(message = "Highlighting the selected region",{
          ht <- HighlightTrack(trackList = df,
                               start = y$start, width = y$width,
                               chromosome = chr, alpha=0.5)
        })
      }else ht <- NULL
      return(ht)
    }
  })
  goi_trackplot_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected) &&
       !is.null(goi_promoter_position_clustering()) && 
       !is.null(goi_gene_position_clustering()) && 
       !is.null(gtrack_clustering()) &&
       !is.null(input$igv_uprange_clustering)){
      library(Gviz)
      if(!is.null(highlight_trackplot_clustering())){
        plot<- plotTracks(list(gtrack_clustering(), highlight_trackplot_clustering()),
                          from = input$igv_uprange_clustering[1], 
                          to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                          type="hist",cex.main = 1.25)
      }else{
        df <- data_track_clustering()
        df[["gtrack"]] <- gtrack_clustering()
        plot<- plotTracks(df,
                          from = input$igv_uprange_clustering[1], 
                          to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                          type="hist",cex.main = 1.25)
      }
      return(plot)
    }
  })
  output$trackplot_goi_clustering <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$clustering_kmeans_extract_table_rows_selected)  &&
         input$Species_clustering != "not selected"){
        goi_trackplot_clustering()
      }
    })
  })
  
  output$download_clustering_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(highlight_trackplot_clustering())){
          plotTracks(list(gtrack_clustering(), highlight_trackplot_clustering()),
                     from = input$igv_uprange_clustering[1], 
                     to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                     type="hist")
        }else{
          df <- data_track_clustering()
          df[["gtrack"]] <- gtrack_clustering()
          plotTracks(df,
                     from = input$igv_uprange_clustering[1], 
                     to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                     type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  ## Enrichment viewer------
  output$Spe_rnaseq2_enrich <- renderText({
    if(input$Species_enrich == "not selected") print("Please select 'Species'")
  })
  observeEvent(input$goButton_enrich,({
    updateSelectInput(session,inputId = "Species_enrich","Species_enrich",species_list, selected = "Homo sapiens (hg19)")
  }))
  output$Spe_dist_enrich <- renderText({
    if(input$Species_enrich == "not selected") print("Please select 'Species'")
  })
  observeEvent(pre_bws_enrich(),({
    updateSelectizeInput(session,inputId = "sample_order_enrich","Sample order:",
                         choices =  names(pre_bws_enrich()),
                         selected = names(pre_bws_enrich()))
  }))
  output$sample_order_enrich_comb1 <- renderUI({
    if(!is.null(pre_integrated_additional1_enrich())){
      selectInput(inputId = "sample_order_enrich_comb1","Sample order (red):",
                  choices =  names(pre_integrated_additional1_enrich()),
                  selected = names(pre_integrated_additional1_enrich()),multiple = TRUE)
    }
  })
  output$sample_order_enrich_comb2 <- renderUI({
    if(!is.null(pre_integrated_additional2_enrich())){
      selectInput(inputId = "sample_order_enrich_comb2","Sample order (blue):",
                  choices =  names(pre_integrated_additional2_enrich()),
                  selected = names(pre_integrated_additional2_enrich()),multiple = TRUE)
    }
  })
  output$sample_order_enrich_comb3 <- renderUI({
    if(!is.null(pre_integrated_additional3_enrich())){
      selectInput(inputId = "sample_order_enrich_comb3","Sample order (green):",
                  choices =  names(pre_integrated_additional3_enrich()),
                  selected = names(pre_integrated_additional3_enrich()),multiple = TRUE)
    }
  })
  output$sample_order_enrich_comb4 <- renderUI({
    if(!is.null(pre_integrated_additional4_enrich())){
      selectInput(inputId = "sample_order_enrich_comb4","Sample order (purple):",
                  choices =  names(pre_integrated_additional4_enrich()),
                  selected = names(pre_integrated_additional4_enrich()),multiple = TRUE)
    }
  })
  
  #--------
  output$Spe_motif_enrich <- renderText({
    if(input$Species_enrich == "not selected") print("Please select 'Species'")
  })
  output$Spe_GREAT_enrich <- renderText({
    if(input$Species_enrich == "not selected") print("Please select 'Species'")
  })
  
  
  txdb_enrich <- reactive({
    if(input$Species_enrich != "not selected"){
      return(txdb_function(Species = input$Species_enrich))
    }
  })
  org_enrich <- reactive({
    return(org(Species = input$Species_enrich))
  })
  
  Enrich_peak_call_files <- reactive({
    if(is.null(input$enrich_data_file)){
      if(input$goButton_enrich > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
        files[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$enrich_data_file[, 1])){
        file <- input$enrich_data_file[[nr, 'datapath']]
        name <- c(name, bed_name(input$enrich_data_file[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  
  output$enrichment_input <- DT::renderDataTable({
    uploaded_files = names(Enrich_peak_call_files())
    as.data.frame(uploaded_files)
  })
  
  #Annotation-----------
  enrichdistribution <- reactive({
    return(ChIPpeakAnno::genomicElementDistribution(GRangesList(Enrich_peak_call_files()), 
                                                    TxDb = txdb_enrich()))
  })
  output$input_peak_distribution_enrich <- renderPlot({
    if(!is.null(Enrich_peak_call_files()) && !is.null(txdb_enrich())){
      withProgress(message = "Preparing peak distribution",{
        enrichdistribution()$plot
      })
    }
  })
  output$selected_enrich_peak_distribution <- renderPlot({
    if(!is.null(Enrich_peak_call_files()) && !is.null(txdb_enrich()) && !is.null(input$annotation_select_enrich)){
      donut_replot(dplyr::filter(enrichdistribution()$plot[[1]], source == input$annotation_select_enrich))
    }
  })
  output$download_selected_enrich_annotation_table <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"all_annotation",".zip")
    },
    content = function(fname){
      fs <- c()
      dir.create("all_annotation/",showWarnings = FALSE)
      for(name in names(selected_annoData_table_enrich())){
        file <- paste0("all_annotation/",name,".txt")
        file_pdf <- paste0("all_annotation/",name,".pdf")
        fs <- c(fs,file,file_pdf)
        df <- apply(selected_annoData_table_enrich()[[name]],2,as.character)
        pdf(file_pdf, height = 4.5, width = 6)
        print(donut_replot(dplyr::filter(enrichdistribution()$plot[[1]], source == name)))
        dev.off()
        write.table(df,file,col.names = T,row.names = F,sep = "\t",quote = F)
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  selected_annoData_table_enrich <- reactive({
    withProgress(message = "Preparing annotation",{
      peak_list <- enrichdistribution()$peaks
      df <- list()
      for(name in names(peak_list)){
        peak <- ChIPseeker::annotatePeak(peak_list[[name]],TxDb = txdb_enrich()) %>% as.data.frame()
        my.symbols <- peak$geneId
        gene_IDs<-id_convert(my.symbols,input$Species_enrich,type="ENTREZID")
        colnames(gene_IDs) <- c("geneId","NearestGene")
        data <- merge(gene_IDs,peak,by="geneId")
        data <- data[,2:10] %>% as.data.frame()
        df[[name]] <- data
      }
      return(df)
    })
  })
  output$enrich_annotation <- DT::renderDT({
    if(!is.null(Enrich_peak_call_files()) && !is.null(txdb_enrich()) && !is.null(enrichdistribution()) &&
       input$Species_enrich != "not selected" && !is.null(input$annotation_select_enrich)){
      selected_annoData_table_enrich()[[input$annotation_select_enrich]] %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  output$annotation_select_enrich <- renderUI({
    if(!is.null(Enrich_peak_call_files()) && !is.null(txdb_enrich()) && !is.null(enrichdistribution()) &&
       input$Species_enrich != "not selected" ){
      selectInput('annotation_select_enrich', 'Select bed file', names(enrichdistribution()$peaks),multiple = F)
    }
  })
  
  output$download_input_peak_distribution_enrich = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"_annotation_from_enrichment_viewer",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download, please wait",{
        fs <- c()
        peak_list <- enrichdistribution()$peaks
        dir.create("Annotation/",showWarnings = FALSE)
        for(name in names(peak_list)){
          file_name <- paste0("Annotation/",name,".txt")
          peak <- selected_annoData_table_enrich()[[name]]
          fs <- c(fs, file_name)
          write.table(peak, file_name, row.names = F, col.names=TRUE, sep = "\t", quote = F)
        }
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$enrich_pdf_width
        plot_name <- paste0("Annotation/","annotation.pdf")
        fs <- c(fs, plot_name)
        pdf(plot_name, height = pdf_height, width = pdf_width)
        print(enrichdistribution())
        dev.off()
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  
  #Enrichment viewer functional analysis------
  output$Gene_set_enrich <- renderUI({
    selectInput('Gene_set_enrich', 'Gene Set', gene_set_list)
  })
  
  Hallmark_set_enrich <- reactive({
    if(!is.null(input$Gene_set_enrich)){
      return(GeneList_for_enrichment(Species = input$Species_enrich, Gene_set = input$Gene_set_enrich, org = org_enrich()))
    }
  })
  
  enrichment_enricher_enrich <- reactive({
    data <- Enrich_peak_call_files()
    H_t2g <- Hallmark_set_enrich()
    df <- list()
    source <- ref_for_GREAT(input$Species_enrich)
    for(name in names(data)){
      res = rGREAT::great(gr = data[[name]],gene_sets = gene_list_for_enrichment_genome(H_t2g,input$Species_enrich),source)
      df[[name]] <- res
    }
    return(df)
  })
  
  enrichment_1_1_enrich <- reactive({
    df <- enrichment_enricher_enrich()
    if(!is.null(input$Gene_set_enrich) && input$Species_enrich != "not selected" && !is.null(df)){
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                          "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          if(length(group1$id) != 0){
            group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
            group1$Group <- name
            data <- rbind(data, group1)
          }else group1 <-NULL
        }else group1 <- NULL
      }
      if(length(data$id) != 0) {
        return(data)
      }else return(NULL)
    }else{return(NULL)}
  })
  
  pair_enrich1_H_enrich <- reactive({
    if(!is.null(input$Gene_set_enrich) && input$Species_enrich != "not selected"){
      df <- enrichment_enricher_enrich()
      data3 <- Enrich_peak_call_files()
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          if(length(group1$id) != 0){
            group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
            if (length(group1$p_adjust_hyper) > input$enrich_showCategory){
              group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
              group1 <- group1[1:input$enrich_showCategory,]
            }}else group1 <- NULL
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                             "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      if(length(data$Description) != 0){
        data$Group <- gsub("_", " ", data$Group)
        for(i in 1:length(data$Group)){
          data$Group[i] <- paste(strwrap(data$Group[i], width = 15),collapse = "\n")
        }
        data$Group <- gsub(" \\(", "\n\\(", data$Group)
        data$Description <- gsub("_", " ", data$Description)
        data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
        data$x <- gsub(":","", data$x)
        data <- dplyr::arrange(data, x)
        idx <- order(data[["x"]], decreasing = FALSE)
        data$Description <- factor(data$Description,
                                   levels=rev(unique(data$Description[idx])))
        p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                        geom_point() +
                        scale_color_continuous(low="red", high="blue",
                                               guide=guide_colorbar(reverse=TRUE)) +
                        scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                        scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      }else p1 <- NULL
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment_enrich <- renderPlot({
    if(input$Species_enrich != "not selected" && !is.null(Enrich_peak_call_files())){
      dotplot_for_output(data = Enrich_peak_call_files(),
                         plot_genelist = pair_enrich1_H_enrich(), Gene_set = input$Gene_set_enrich, 
                         Species = input$Species_enrich)
    }
  })
  
  
  output$download_enrich_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$Gene_set_enrich,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H_enrich()
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_enrich_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$Gene_set_enrich,".txt")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1_enrich()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  pair_enrich_table_enrich <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1_enrich()), 
                            H_t2g = Hallmark_set_enrich(), Gene_set = input$Gene_set_enrich))
  })
  
  output$enrich_enrichment_result <- DT::renderDataTable({
    if(!is.null(input$Gene_set_enrich) && 
       input$Species_enrich != "not selected" && !is.null(Enrich_peak_call_files())){
      as.data.frame(enrichment_1_1_enrich())
    }
  })
  
  
  
  output$whichGeneSet_enrich <- renderUI({
    if(input$Species_enrich != "not selected" && !is.null(Enrich_peak_call_files())){
      if(!is.null(input$intersection_enrich_fun)){
        group <- input$intersection_enrich_fun
        if(dim(as.data.frame(enrichment_1_1_enrich()))[1] != 0){
          set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1_enrich()), Group == group)$id)
          selectInput('Pathway_list_enrich', 'Pathway list', set_list)
        }
      }
    }
  })
  
  region_gene_associate_enrich <- reactive({
    if(!is.null(input$Pathway_list_enrich)){
      set_list <- input$Pathway_list_enrich
      df <- enrichment_enricher_enrich()
      data3 <- Enrich_peak_call_files()
      data <- data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))
          group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      data <- apply(data,2,as.character)
      return(data)
    }
  })
  output$region_gene_enrich_associations <- DT::renderDT({
    if(!is.null(input$Pathway_list_enrich) && 
       !is.null(input$intersection_enrich_fun) && input$Species_enrich != "not selected"){
      region_gene_associate_enrich()
    }
  })
  output$whichGenes_enrich <- renderUI({
    if(input$Species_enrich != "not selected" && !is.null(input$Gene_set_enrich)){ 
      selectInput('intersection_enrich_fun', 'Select intersection', names(enrichment_enricher_enrich()))
    }
  })
  
  region_gene_associate_plot_enrich <- reactive({
    set_list <- input$Pathway_list_enrich
    df <- enrichment_enricher_enrich()
    name <- input$intersection_enrich_fun
    if(!is.null(df[[name]])) {
      res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
    }else res <- NULL
    return(res)
  })
  output$region_gene_associations_enrich_plot <- renderPlot({
    if(!is.null(input$Pathway_list_enrich) && 
       !is.null(input$intersection_enrich_fun) && input$Species_enrich != "not selected"){
      region_gene_associate_plot_enrich()
    }
  })
  output$download_enrich_region_gene_associations_plot = downloadHandler(
    filename = function(){
      paste0("region_gene_associations-",input$intersection_enrich_fun,"-",input$Pathway_list_enrich,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$enrich_pdf_width
        set_list <- input$Pathway_list_enrich
        df <- enrichment_enricher_enrich()
        name <- input$intersection_enrich_fun
        pdf(file, height = pdf_height, width = pdf_width)
        rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_enrich_region_gene_associations = downloadHandler(
    filename = function() {
      paste0("region_gene_associations-",input$intersection_enrich_fun,"-",input$Pathway_list_enrich,".txt")
    },
    content = function(file){write.table(region_gene_associate_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  ##heatmap----------
  Enrich_peak_call_files_locus <- reactive({
    data <- Enrich_peak_call_files_geneId()
    df <- list()
    for(name in names(data)){
      data2 <- data[[name]]
      data2_gr <- with(data2, GRanges(seqnames = seqnames, 
                                      ranges = IRanges(start,end),
                                      strand = strand))
      names(data2_gr) <- paste0(data2$ENTREZID,"_",data2$locus)
      df[[name]] <- data2_gr
    }
    return(df)
  })
  
  integrated_legend_enrich <- reactive({
    lgd <- lgd(files2 = Enrich_peak_call_files_locus())
    return(lgd)
  })
  
  integrated_heatlist_enrich <- reactive({
    if(input$integrated_heatmapButton_enrich > 0 && updateCounter_int_enrich$i > 0){
      ht_list <- NULL
      if(!is.null(integrated_heatmap_add1_enrich())) ht_list <- ht_list + integrated_heatmap_add1_enrich()[["heatmap"]]
      if(!is.null(integrated_heatmap_add2_enrich())) ht_list <- ht_list + integrated_heatmap_add2_enrich()[["heatmap"]]
      if(!is.null(integrated_heatmap_add3_enrich())) ht_list <- ht_list + integrated_heatmap_add3_enrich()[["heatmap"]]
      if(!is.null(integrated_heatmap_add4_enrich())) ht_list <- ht_list + integrated_heatmap_add4_enrich()[["heatmap"]]
      if(!is.null(rnaseq_DEGs2_enrich())) ht_list <- ht_list + rnaseq_DEGs_heatmap_enrich()
      if(!is.null(rnaseq_count2_enrich())) ht_list <- ht_list + rnaseq_count_heatmap_enrich()
      return(ht_list)
    }else return(NULL)
  })
  updateCounter_int_enrich <- reactiveValues(i = 0)
  
  observe({
    input$integrated_heatmapButton_enrich
    isolate({
      updateCounter_int_enrich$i <- updateCounter_int_enrich$i + 1
    })
  })
  
  
  #Restart
  defaultvalues_enrich <- observeEvent(integrated_heatlist_enrich(), {
    isolate(updateCounter_int_enrich$i == 0)
    updateCounter_int_enrich <<- reactiveValues(i = 0)
  }) 
  integrated_heatlist_enrich_plot <- reactive({
    if(input$integrated_heatmapButton_enrich > 0 && !is.null(Enrich_peak_call_files_locus()) &&
       !is.null(integrated_heatlist_enrich())){
      return(draw(integrated_heatlist_enrich(),annotation_legend_list = list(integrated_legend_enrich()),
                  heatmap_legend_side = "bottom", ht_gap = unit(2, "mm")))
    }
  })
  output$integrated_heatmap_enrich <- renderPlot({
    if(input$integrated_heatmapButton_enrich > 0 && !is.null(Enrich_peak_call_files_locus()) &&
       !is.null(integrated_heatlist_enrich())){
      integrated_heatlist_enrich_plot()
    }
  })
  output$download_integrated_heatmap_enrich = downloadHandler(
    filename = "Integrated_heatmap.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(integrated_heatlist_enrich_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$rnaseq_count_enrich <- renderUI({
    if(input$Species_enrich != "not selected"){
      fileInput("pair_rnaseq_count_enrich",
                "Select RNA-seq normalized count files",
                accept = c("txt","csv","xlsx"),
                multiple = TRUE,
                width = "90%")
    }
  })
  output$rnaseq_DEGs_enrich <- renderUI({
    if(input$Species_enrich != "not selected"){
      fileInput("pair_rnaseq_DEGs_enrich",
                "Select RNA-seq DEG result files",
                accept = c("txt","csv","xlsx"),
                multiple = TRUE,
                width = "90%")
    }
  })
  
  rnaseq_count_enrich <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      tmp <- input$pair_rnaseq_count_enrich
      upload = list()
      if(is.null(input$pair_rnaseq_count_enrich) && input$goButton_enrich > 0 )  {
        tmp = "data/RNAseq_count.txt"
        upload[["rna"]]<- read_df(tmp)
        return(upload)
      }else if(is.null(tmp)) {
        return(NULL)
      }else{
        return(read_dfs(tmp))
      }
    })
  })
  rnaseq_DEGs_enrich <- reactive({
    withProgress(message = "Importing DEG result files, please wait",{
      tmp <- input$pair_rnaseq_DEGs_enrich
      upload = list()
      if(is.null(input$pair_rnaseq_DEGs_enrich) && input$goButton_enrich > 0 )  {
        tmp = "data/RNAseq.txt"
        upload[["rna"]]<- read_df(tmp)
        return(upload)
      }else if(is.null(tmp)) {
        return(NULL)
      }else{
        return(read_dfs(tmp))
      }
    })
  })
  output$pre_zscoring_enrich <- renderUI({
    if(length(names(rnaseq_count_enrich())) > 1){
      selectInput("pre_zscoring_enrich","pre_zscoring for multiple normalized count data",choices = c(TRUE,FALSE),selected = TRUE)
    }
  })
  gene_type_rnaseq_DEGs_enrich <- reactive({
    files <- rnaseq_DEGs_enrich()
    return(gene_type_for_integrated_heatmap(files=files,org=org_enrich(),Species=input$Species_enrich))
  })
  rnaseq_DEGs2_enrich <- reactive({
    files <- rnaseq_DEGs_enrich()
    return(rnaseqDEGs_for_integrated_heatmap(files = files,Species=input$Species_enrich,gene_type = gene_type_rnaseq_DEGs_enrich()))
  })
  gene_type_rnaseq_counts_enrich <- reactive({
    files <- rnaseq_count_enrich()
    return(gene_type_for_integrated_heatmap(files=files,org=org_enrich(),Species=input$Species_enrich))
  })
  rnaseq_count2_enrich <- reactive({
    files <- rnaseq_count_enrich()
    return(rnaseqCounts_for_integrated_heatmap(files = files,Species=input$Species_enrich,gene_type = gene_type_rnaseq_counts_enrich(),pre_zscoring = input$pre_zscoring_enrich))
  })
  observeEvent(input$pair_rnaseq_DEGs_enrich, ({
    updateCollapse(session,id =  "z-score_count_enrich", open="Uploaded_DEGs_enrich")
  }))
  observeEvent(input$pair_rnaseq_count_enrich, ({
    updateCollapse(session,id =  "z-score_count_enrich", open="z-score_multiple_count_enrich_panel")
  }))
  output$rnaseq_count_output_enrich <- renderDataTable({
    if(input$Species_enrich != "not selected" && !is.null(rnaseq_count_enrich())){
      rnaseq_count2_enrich()
    }
  })
  output$rnaseq_DEGs_output_enrich <- renderDataTable({
    if(input$Species_enrich != "not selected" && !is.null(rnaseq_DEGs_enrich())){
      rnaseq_DEGs2_enrich()
    }
  })
  Enrich_peak_call_files_geneId <- reactive({
    peaks <- Enrich_peak_call_files()
    df <- list()
    for(name in names(peaks)){
      data <- peaks[[name]]
      data2 <- as.data.frame(ChIPseeker::as.GRanges(ChIPseeker::annotatePeak(peak = data, TxDb = txdb_enrich())))
      data2 <- data2 %>% dplyr::filter(!is.na(geneId))
      data2$locus <- paste0(data2$seqnames,":",data2$start,"-",data2$end)
      colnames(data2)[which(colnames(data2) == "geneId")] <- "ENTREZID"
      df[[name]] <- data2
    }
    return(df)
  })
  rnaseq_count_heatmap_enrich <- reactive({
    rna <-  rnaseq_count2_enrich()
    peaks <- Enrich_peak_call_files_geneId()
    rna$ENTREZID <- rownames(rna)
    m_z <- data.frame(matrix(rep(NA, length(colnames(rna))), nrow=1))[numeric(0), ]
    for(name in names(peaks)){
      data <- peaks[[name]]
      data_m <- merge(rna, data, by="ENTREZID",all=T)
      data_m <- data_m %>% dplyr::filter(!is.na(locus))
      rownames(data_m) <- paste0(data_m$ENTREZID,"_",data_m$locus)
      data_m <- data_m[,2:length(colnames(rna))]
      data_m[is.na(data_m)] <- 0
      m_z <- rbind(m_z,data_m)
    }
    cond <- gsub(".+\\.", "", colnames(m_z))
    cond <- gsub("\\_.+$", "", cond)
    cond <- factor(cond, levels = unique(cond), ordered = TRUE)
    cond_color <- rainbow_hcl(length(unique(cond)),c=100)
    names(cond_color) <- unique(cond)
    if(length(names(rnaseq_count_enrich())) == 1){
      file_name <- NULL
      file_name_color <- NULL
    }else{
      file_name <- gsub("\\..+$", "", colnames(m_z))
      file_name <- factor(file_name, levels = unique(file_name), ordered = TRUE)
      file_name_color <- rainbow_hcl(length(file_name))
      names(file_name_color) <- file_name
    }
    mat <- integrated_heatmap_add1_enrich()[["mat"]]
    withProgress(message = "Heatmap of RNA-seq count data",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),],name = "RNA-seq\nz-score", 
                    top_annotation = HeatmapAnnotation(files = file_name, condition = cond,
                                                       col = list(files = file_name_color,
                                                                  condition = cond_color)),
                    show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "cm"),use_raster = TRUE)
    })
    return(ht)
  })
  rnaseq_DEGs_heatmap_enrich <- reactive({
    rna <-  rnaseq_DEGs2_enrich()
    peaks <- Enrich_peak_call_files_geneId()
    rna$ENTREZID <- rownames(rna)
    m_z <- data.frame(matrix(rep(NA, length(colnames(rna))), nrow=1))[numeric(0), ]
    for(name in names(peaks)){
      data <- peaks[[name]]
      data_m <- merge(rna, data, by="ENTREZID",all=T)
      data_m <- data_m %>% dplyr::filter(!is.na(locus))
      rownames(data_m) <- paste0(data_m$ENTREZID,"_",data_m$locus)
      data_m <- data_m[,1:length(colnames(rna))]
      data_m[is.na(data_m)] <- 0
      m_z <- rbind(m_z,data_m)
    }
    for(i in 1:length(colnames(rna))){
      m_z[,i] <- -1 * as.numeric(m_z[,i])
    }
    m_z[m_z > 5]<- 5
    m_z[m_z < -5]<- -5
    colnames(m_z) <- gsub("\\.log2F.+$", "", colnames(m_z))
    mat <- integrated_heatmap_add1_enrich()[["mat"]]
    withProgress(message = "Heatmap of RNA-seq log2FoldChange",{
      ht <- Heatmap(as.matrix(m_z)[rownames(mat),2:length(colnames(rna))],name = "RNA-seq\nlog2FC", 
                    show_row_names = FALSE, width = unit(2.5, "cm"), column_names_gp = grid::gpar(fontsize = 9),
                    use_raster = TRUE,column_names_side = "top",show_column_dend = FALSE,
                    col = c("blue","white","gold"))
      return(ht)
    })
  })
  
  output$integrated_bw1_enrich <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_1_enrich",
                  "Select bigwig files (red)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_1_enrich",
                "Select bigwig files (red)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw2_enrich <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_2_enrich",
                  "Option: Select additional bigwig files (blue)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_2_enrich",
                "Option: Select additional bigwig files (blue)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw3_enrich <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_3_enrich",
                  "Option: Select additional bigwig files (green)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_3_enrich",
                "Option: Select additional bigwig files (green)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  output$integrated_bw4_enrich <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("integrated_bw_4_enrich",
                  "Option: Select additional bigwig files (purple)",
                  choices = bw,multiple=T)
    }else{
      fileInput("integrated_bw_4_enrich",
                "Option: Select additional bigwig files (purple)",
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  pre_integrated_additional1_enrich <-reactive({
    if(!is.null(input$integrated_bw_1_enrich)){
      return(pre_integrated_additional(input$integrated_bw_1_enrich))
    }else{
      if(input$goButton_enrich > 0 ){
        file <- c("data/bigwig/A_1.BigWig",
                  "data/bigwig/B_1.BigWig")
        names(file) <- c("A_1","B_1")
        return(file)
      }
    }
  })
  integrated_additional1_enrich <- reactive({
    return(bws_ordering(bws=pre_integrated_additional1_enrich(),sample_order=input$sample_order_enrich_comb1))
  })
  pre_integrated_additional2_enrich <-reactive({
    return(pre_integrated_additional(input$integrated_bw_2_enrich))
  })
  integrated_additional2_enrich <- reactive({
    return(bws_ordering(bws=pre_integrated_additional2_enrich(),sample_order=input$sample_order_enrich_comb2))
  })
  pre_integrated_additional3_enrich <-reactive({
    return(pre_integrated_additional(input$integrated_bw_3_enrich))
  })
  integrated_additional3_enrich <- reactive({
    return(bws_ordering(bws=pre_integrated_additional3_enrich(),sample_order=input$sample_order_enrich_comb3))
  })
  pre_integrated_additional4_enrich <-reactive({
    return(pre_integrated_additional(input$integrated_bw_4_enrich))
  })
  integrated_additional4_enrich <- reactive({
    return(bws_ordering(bws=pre_integrated_additional4_enrich(),sample_order=input$sample_order_enrich_comb4))
  })
  integrated_heatmap_add1_enrich <- reactive({
    if(!is.null(integrated_additional1_enrich())){
      h <- batch_heatmap(files2 = Enrich_peak_call_files_locus(),files_bw = integrated_additional1_enrich(),
                         color = c("white","red"),signal = "red")
      return(h)
    }
  })
  integrated_heatmap_add2_enrich <- reactive({
    if(!is.null(integrated_additional2_enrich())){
      h <- batch_heatmap(files2 = Enrich_peak_call_files_locus(),files_bw = integrated_additional2_enrich(),
                         color = c("white","darkblue"),signal = "darkblue")
      return(h)
    }
  })
  integrated_heatmap_add3_enrich <- reactive({
    if(!is.null(integrated_additional3_enrich())){
      h <- batch_heatmap(files2 = Enrich_peak_call_files_locus(),files_bw = integrated_additional3_enrich(),
                         color = c("white","darkgreen"),signal = "green")
      return(h)
    }
  })
  integrated_heatmap_add4_enrich <- reactive({
    if(!is.null(integrated_additional4_enrich())){
      h <- batch_heatmap(files2 = Enrich_peak_call_files_locus(),files_bw = integrated_additional4_enrich(),
                         color = c("white","purple"),signal = "purple")
      return(h)
    }
  })
  
  #enrichment withRNAseq-------
  #enrich regulatory potential------
  output$enrich_RNA <- renderUI({
    if(!is.null(Enrich_peak_call_files())){
      selectInput("enrich_RNA", "Select Peaks",
                  c(names(Enrich_peak_call_files())),
                  selected = c(names(Enrich_peak_call_files())),multiple = F)
    }
  })
  output$enrich_select_RNA <- renderUI({
    if(!is.null(Enrich_peak_call_files())){
      enrich <- c(names(Enrich_peak_call_files()))
      selectInput("enrich_select_RNA", "Select Peak",
                  c(names(Enrich_peak_call_files())),
                  selected = enrich,multiple = T)
    }
  })
  react_enrich_select_RNA_debounce <- reactive({
    if(!is.null(input$enrich_select_RNA)) return(input$enrich_select_RNA) else return(NULL)
  })
  enrich_select_RNA_debounce <- debounce(react_enrich_select_RNA_debounce, 1000)
  
  peak_enrich_grange_RNA  <- reactive({
    if(!is.null(input$enrich_RNA)){
      return(Enrich_peak_call_files()[[input$enrich_RNA]])
    }
  })
  
  
  
  output$enrichRNAseqresult <- renderUI({
    if(input$RNAseq_data_type_enrich == "Result") {
      label <- "Select RNA-seq DEG result file"
    }else{
      if(input$RNAseq_data_type_enrich == "List") label <- "Select Gene list files" else 
        label <- "Select RNA-seq raw count file"
    }
    fileInput("enrich_DEG_result",
              label,
              accept = c("txt","csv","xlsx"),
              multiple = TRUE,
              width = "80%")
  })
  output$pre_genelist_input_choice_enrich <- renderUI({
    if(!is.null(pre_RNAseq_file_enrich())){
      list <- unique(pre_RNAseq_file_enrich()[,2])
      if(length(list) > 20 && length(input$enrich_DEG_result[, 1]) == 1){ 
        selectInput("pre_genelist_input_choice_enrich","Group name",c("File name","Second column"),selected = "File name",multiple = F)
      }else return(NULL)
    }
  })
  output$genelist_input_choice_enrich <- renderUI({
    if(!is.null(pre_RNAseq_file_enrich())){
      list <- unique(pre_RNAseq_file_enrich()[,2])
      if(!is.null(input$pre_genelist_input_choice_enrich)){
        if(input$pre_genelist_input_choice_enrich == "File name"){
          file_name <- gsub(paste0("\\.",tools::file_ext(input$enrich_DEG_result[[1, 'datapath']]),"$"), "", input$enrich_DEG_result[1,]$name)
          list <- file_name
        }
      }
      list <- sort(list)
      selectInput("genelist_input_choice_enrich","Groups (RNA)",list, selected = list, multiple = TRUE)
    }
  })
  pre_RNAseq_file_enrich <- reactive({
    withProgress(message = "Importing RNA-seq data, please wait",{
      tmp <- input$enrich_DEG_result$datapath 
      if(is.null(input$enrich_DEG_result) && input$goButton_enrich > 0 )  {
        if(input$RNAseq_data_type_enrich == "Result") tmp = "data/RNAseq.txt" else{
          if(input$RNAseq_data_type_enrich == "List") tmp = "data/genelist.txt" else tmp = "data/RNAseq_count.txt"
        } 
      }
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(input$RNAseq_data_type_enrich == "List"){
          upload = list()
          name = c()
          if(is.null(input$enrich_DEG_result) && input$goButton_enrich > 0 ){
            upload["genelist"] <- list(read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = ""))
            name <- "genelist"
          }else{
            for(nr in 1:length(input$enrich_DEG_result[, 1])){
              if(tools::file_ext(input$enrich_DEG_result[[nr, 'datapath']]) == "xlsx") {
                df2 <- read_xlsx(input$enrich_DEG_result[[nr, 'datapath']])
                df2 <- as.data.frame(df2)
                df <- try(data.frame(row.names = df2[,1]),silent = T)
                if(class(df) != "try-error") {
                  if(dim(df2)[2] == 2){
                    df <- data.frame(row.names = df2[,1],a = df2[,2])
                    colnames(df)[1] <- colnames(df2)[2]
                  }else{
                    rownames(df2) <- df2[,1]
                    df <- df2[,-1]
                    colnames(df) <- gsub("-",".",colnames(df))
                  }
                }
              }
              if(tools::file_ext(input$enrich_DEG_result[[nr, 'datapath']]) == "csv") df <- read.csv(input$enrich_DEG_result[[nr, 'datapath']], header=TRUE, sep = ",",quote = "")
              if(tools::file_ext(input$enrich_DEG_result[[nr, 'datapath']]) == "txt" || 
                 tools::file_ext(input$enrich_DEG_result[[nr, 'datapath']]) == "tsv") df <- read.table(input$enrich_DEG_result[[nr, 'datapath']], header=TRUE, sep = "\t",quote = "")
              file_name <- gsub(paste0("\\.",tools::file_ext(input$enrich_DEG_result[[nr, 'datapath']]),"$"), "", input$enrich_DEG_result[nr,]$name)
              name <- c(name, file_name)
              upload[nr] <- list(df)
            }
          }
          names(upload) <- name
          if(length(names(upload)) == 1){
            tmp <- upload[[name]]
            rownames(tmp) = gsub("\"", "", rownames(tmp))
            if(str_detect(colnames(tmp)[1], "^X\\.")){
              colnames(tmp) = str_sub(colnames(tmp), start = 3, end = -2) 
            }
            if(rownames(tmp)[1] == 1){
              if(dim(tmp)[2] >= 2){
                tmp <- data.frame(Gene = tmp[,1], Group = tmp[,2])
              }else{
                tmp <- data.frame(Gene = tmp[,1], 
                                  Group = gsub(paste0("\\.",tools::file_ext(input$enrich_DEG_result[[1, 'datapath']]),"$"), "", input$enrich_DEG_result[1,]$name))
              }
            }else{
              tmp <- data.frame(Gene = rownames(tmp), Group = tmp[,1])
            }
          }else{
            df2 <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
            for(file in names(upload)){
              df <- upload[[file]]
              if(rownames(df)[1] == 1){
                df[,1] = gsub("\"", "", df[,1])
                df <- data.frame(Gene = df[,1], Group = file)
              }else{
                rownames(df) = gsub("\"", "", rownames(df))
                df <- data.frame(Gene = rownames(df), Group = file)
              }
              df2 <- rbind(df2,df)
            }
            tmp <- df2
          }
          tmp$Group <- gsub(":","-",tmp$Group)
          return(tmp)
        }else{
          if(tools::file_ext(tmp) == "xlsx") {
            df2 <- readxl::read_xlsx(tmp) 
            df2 <- as.data.frame(df2)
            df <- try(data.frame(row.names = df2[,1]),silent = T)
            if(class(df) != "try-error") {
              if(dim(df2)[2] == 2){
                df <- data.frame(row.names = df2[,1],a = df2[,2])
                colnames(df)[1] <- colnames(df2)[2]
              }else{
                rownames(df2) <- df2[,1]
                df <- df2[,-1]
                colnames(df) <- gsub("-",".",colnames(df))
              }
            }
          }
          if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
          if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
          rownames(df) = gsub("\"", "", rownames(df))
          if(length(colnames(df)) != 0){
            if(str_detect(colnames(df)[1], "^X\\.")){
              colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
            }
          }
          return(df)
        }
      }
    })
  })
  output$enrich_RNAseq_order <- renderUI({
    if(!is.null(pre_RNAseq_file_enrich())){
      if(input$RNAseq_data_type_enrich != "Result" && input$RNAseq_data_type_enrich != "List"){
        if(length(grep("log2FoldChange", colnames(pre_RNAseq_file_enrich()))) != 0) validate("Uploaded data is a DEG_result file. Please Select 'DEG_Result' from 'Input'.")
        selectInput("enrich_RNAseq_order","Select samples:",
                    choices = colnames(pre_RNAseq_file_enrich()),selected = colnames(pre_RNAseq_file_enrich()),multiple = T)
      }
    }
  })
  RNAseq_file_enrich <- reactive({
    if(!is.null(pre_RNAseq_file_enrich())){
      if(input$RNAseq_data_type_enrich != "Result"){
        count <- pre_RNAseq_file_enrich()
        order <- input$enrich_RNAseq_order
        data <- try(count[,order])
        if(length(data) == 1){
          if(class(data) == "try-error") validate("")
        }
        return(data)
      }
    }
  })
  updateCounter_DEGanalysis_enrich <- reactiveValues(i = 0)
  
  observe({
    input$DEGanalysis_Button_enrich
    isolate({
      updateCounter_DEGanalysis_enrich$i <- updateCounter_DEGanalysis_enrich$i + 1
    })
  })
  
  #Restart
  observeEvent(RNAseq_file_enrich(), {
    isolate(updateCounter_DEGanalysis_enrich$i == 0)
    updateCounter_DEGanalysis_enrich <<- reactiveValues(i = 0)
  }) 
  pre_RNAseqDEG_enrich <- reactive({
    if(input$DEGanalysis_Button_enrich > 0 && updateCounter_DEGanalysis_enrich$i > 0){
      count <- RNAseq_file_enrich()
      if(!is.null(count)){
        collist <- gsub("\\_.+$", "", colnames(count))
        if(length(unique(collist)) == 2){
          group <- data.frame(con = factor(collist))
          dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
          dds$con <- factor(dds$con, levels = unique(collist))
          dds <- DESeq(dds)
          return(dds)
        }else validate(print(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 2 conditions.")))
      }
    }
  })
  RNAseqDEG_enrich <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if(input$RNAseq_data_type_enrich == "Result"){
        if(length(grep("log2FoldChange", colnames(pre_RNAseq_file_enrich()))) == 0) validate("Error: the uploaded file is in an unexpected format.\nThe file does not contain 'log2FoldChange' or 'padj'.")
        return(pre_RNAseq_file_enrich())
      }else{
        if(input$RNAseq_data_type_enrich != "List"){
          if(input$DEGanalysis_Button_enrich > 0 && updateCounter_DEGanalysis_enrich$i > 0){
            count <- RNAseq_file_enrich()
            collist <- gsub("\\_.+$", "", colnames(count))
            dds <- pre_RNAseqDEG_enrich()
            contrast <- c("con", unique(collist))
            res <- results(dds,  contrast = contrast)
            res <- as.data.frame(res)
            return(res)
          }
        }else {
          if(!is.null(input$genelist_input_choice_enrich)){
            RNAdata <- pre_RNAseq_file_enrich()
            if(!is.null(input$pre_genelist_input_choice_enrich)){
              if(input$pre_genelist_input_choice_enrich == "File name"){
                file_name <- gsub(paste0("\\.",tools::file_ext(input$enrich_DEG_result[[1, 'datapath']]),"$"), "", input$enrich_DEG_result[1,]$name)
                RNAdata$Group <- file_name
              }
            }
            group <- factor(input$genelist_input_choice_enrich, levels = input$genelist_input_choice_enrich)
            RNAdata <- RNAdata %>% dplyr::filter(Group %in% group)
            if(is.element(TRUE, duplicated(RNAdata$Gene)) == TRUE) validate("Duplication of gene names are not allowed.")
            rownames(RNAdata) <- RNAdata$Gene
            return(RNAdata)
          }
        }
      }
    })
  })
  RNAseqDEG_norm_enrich <- reactive({
    if(input$RNAseq_data_type_enrich != "Result" && input$RNAseq_data_type_enrich != "List" && !is.null(pre_RNAseqDEG_enrich())){
      count <- RNAseq_file_enrich()
      collist <- gsub("\\_.+$", "", colnames(count))
      dds <- pre_RNAseqDEG_enrich()
      contrast <- c("con", unique(collist))
      normalized_counts <- counts(dds, normalized=TRUE)
      return(normalized_counts)
    }
  })
  RNAseqDEG_anno_enrich <- reactive({
    if(is.null(pre_RNAseq_file_enrich())) validate ("Please upload 'RNA-seq data'")
    if(input$Species_enrich == "not selected") validate ("Please select 'Species'")
    return(RNAseqDEG_ann(RNAdata=RNAseqDEG_enrich(),Species=input$Species_enrich,gene_type = gene_type_enrich_DEG_result(),input_type=input$RNAseq_data_type_enrich))
  })
  gene_type_enrich_DEG_result <- reactive({
    if(is.null(pre_RNAseq_file_enrich())) validate ("Please upload 'RNA-seq data'")
    if(input$Species_enrich == "not selected") validate ("Please select 'Species'")
    return(gene_type(my.symbols=rownames(RNAseqDEG_enrich()),org=org_enrich(),Species=input$Species_enrich))
  })
  output$enrich_RNAseq_raw <- renderDataTable({
    if(input$RNAseq_data_type_enrich != "Result"){
      if(!is.null(RNAseq_file_enrich())){
        RNAseq_file_enrich() 
      }
    }
  })
  output$RNAseq_condition_enrich <- renderText({
    if(!is.null(pre_RNAseq_file_enrich())){
      paste0(input$RNAseq_cond1_enrich, " vs ", input$RNAseq_cond2_enrich, "\n",
             "Please confirm if the Log2FoldChange in the uploaded result file was calulated as Log2(",input$RNAseq_cond1_enrich,"/",input$RNAseq_cond2_enrich,").")
    }
  })
  output$RNAseq_RPmean_enrich <- renderText({
    paste0("RP > 0 gene is associated with the indicated peak within ",input$peak_distance_enrich," kb from its TSS.\n",
           "RP = 0 gene is not associated with the indicated peak within ",input$peak_distance_enrich," kb from its TSS.")
  })
  output$enrich_DEG_result <- renderDataTable({
    if(!is.null(pre_RNAseq_file_enrich())){
      RNAseqDEG_enrich() 
    }
  })
  RNAseq_name_enrich <- reactive({
    if(!is.null(pre_RNAseq_file_enrich())){
      if(input$RNAseq_data_type_enrich != "Result"){
        if(!is.null(RNAseq_file_enrich())){
          if(input$RNAseq_data_type_enrich == "List"){
            return(c(unique(RNAseqDEG_enrich()$Group)))
          }else{
            collist <- gsub("\\_.+$", "", colnames(RNAseq_file_enrich()))
            cond1 <- paste0(unique(collist)[1],"_high")
            cond2 <- paste0(unique(collist)[2],"_high")
            return(c(cond1, cond2))
          }
        }
      }else{
        cond1 <- paste0(input$RNAseq_cond1_enrich,"_high")
        cond2 <- paste0(input$RNAseq_cond2_enrich,"_high")
        return(c(cond1, cond2))
      }
    }
  })
  
  pre_mmAnno_enrich <- reactive({
    if(input$Species_enrich == "not selected") validate ("Please select 'Species'")
    mmAnno_list <- list()
    for(name in names(Enrich_peak_call_files())){
      peak <- Enrich_peak_call_files()[[name]]
      mAnn <- mmAnno(peak=peak,
                     genomic_region="Genome-wide",
                     txdb=txdb_enrich(),
                     peak_distance=input$peak_distance_enrich,
                     mode=input$RNAseq_mode_enrich,
                     group_name=name,
                     distribution=enrichdistribution()$peaks[[name]],
                     DAR=NULL)
      mmAnno_list[[name]] <- mAnn
    }
    return(mmAnno_list)
  })
  
  
  
  
  
  
  
  
  mmAnno_enrich <- reactive({
    return(pre_mmAnno_enrich()[[input$enrich_RNA]])
  })
  
  RP_enrich <- reactive({
    withProgress(message = "Calculating regulatory potential",{
      if(!is.null(pre_mmAnno_enrich())){
        RP_list <- list()
        for(name in names(pre_mmAnno_enrich())){
          data <- RP_f(mmAnno=pre_mmAnno_enrich()[[name]],txdb=txdb_enrich())
          RP_list[[name]] <- data
        }
        return(RP_list)
      }
      incProgress(1)
    })
  })
  regulatory_potential_enrich <- reactive({
    if(!is.null(RNAseqDEG_anno_enrich()) && !is.null(RP_enrich()) && 
       !is.null(RNAseq_name_enrich())){
      RP_list <- list()
      for(name in names(RP_enrich())){
        data <- regulatory_potential_f(species=input$Species_enrich,data=RNAseqDEG_anno_enrich(),
                                       result_geneRP= RP_enrich()[[name]],DEG_fc=input$DEG_fc_enrich,
                                       DEG_fdr=input$DEG_fdr_enrich,name=RNAseq_name_enrich())
        RP_list[[name]] <- data
      }
      return(RP_list)
    }
  })
  
  output$ks_plot_enrich <- renderPlot({    
    if(!is.null(RNAseqDEG_anno_enrich()) && !is.null(enrich_select_RNA_debounce()) &&
       !is.na(input$DEG_fdr_enrich) && !is.na(input$DEG_fc_enrich) && 
       input$Species_enrich != "not selected"  && !is.null(pre_mmAnno_enrich())){
      regulatory_potential_enrich()[[input$enrich_RNA]]
    }
  })
  ks_tables_enrich <- reactive({
    if(!is.null(RNAseqDEG_anno_enrich()) && !is.null(enrich_select_RNA_debounce()) &&
       !is.na(input$DEG_fdr_enrich) && !is.na(input$DEG_fc_enrich) && 
       input$Species_enrich != "not selected"  && !is.null(pre_mmAnno_enrich())){
      df <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
      for(name in names(regulatory_potential_enrich())){
        df2 <- regulatory_potential_enrich()[[name]]$statistics
        df2$Peak <- name
        df <- rbind(df,df2)
      }
      df <- df %>% dplyr::select(Peak,everything())
      return(df)
    }
  })
  output$ks_table_enrich <- renderDataTable({    
    if(!is.null(RNAseqDEG_anno_enrich()) && !is.null(enrich_select_RNA_debounce()) &&
       !is.na(input$DEG_fdr_enrich) && !is.na(input$DEG_fc_enrich) && 
       input$Species_enrich != "not selected"  && !is.null(pre_mmAnno_enrich())){
      ks_tables_enrich()
    }
  })
  output$DEG_fc_enrich <- renderUI({
    numericInput("DEG_fc_enrich","Fold change cutoff for RNA-seq data",
                 min=0,max=NA,value=1.5,step = 0.5)
  })
  output$DEG_fdr_enrich <- renderUI({
    numericInput("DEG_fdr_enrich","FDR cutoff for RNA-seq data",
                 min=0,max=1, value=0.05,step = 0.001)
  })
  
  output$RNAseqGroup_enrich <- renderUI({
    if(input$Species_enrich != "not selected" && 
       !is.null(regulatory_potential_enrich())){
      if(!is.null(RP_all_table_enrich())){
        selectInput("RNAseqGroup_enrich","Group (Peak:RNAseq)",
                    unique(RP_all_table_enrich()$Group),
                    multiple = FALSE)
      }
    }
  })
  
  RNAseq_boxplot_enrich <- reactive({
    RNA <- RNAseqDEG_anno_enrich()
    cond1 <- gsub("\\_.+$", "", RNAseq_name_enrich()[1])
    cond2 <- gsub("\\_.+$", "", RNAseq_name_enrich()[2])
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    genelist <- list()
    bed_list <- list()
    for(name in enrich_select_RNA_debounce()){
      data <- merge(RNA,RP_enrich()[[name]], by="gene_id",all=T)
      if(input$RNAseq_data_type_enrich != "List") data <- dplyr::filter(data, baseMean != 0)
      data$group <- "Others"
      data$group[data$sumRP > 2] <- "2 < RP"
      data$group[data$sumRP <= 2 & data$sumRP > 1] <- "1 < RP < 2"
      data$group[data$sumRP <= 1 & data$sumRP > 0.1] <- "0.1 < RP < 1"
      data$group[data$sumRP <= 0.1 & data$sumRP > 0.01] <- "0.01 < RP < 0.1"
      data$group[data$sumRP <= 0.01 & data$sumRP > 0] <- "0 < RP < 0.01"
      level <- NULL
      col <- NULL
      if(dim(dplyr::filter(data, group == "Others"))[1] != 0) {
        level <- c(level,"Others")
        col <- c(col,"gray")
      }
      if(dim(dplyr::filter(data, group == "0 < RP < 0.01"))[1] != 0) {
        level <- c(level,"0 < RP < 0.01")
        col <- c(col,"blue")
      }
      if(dim(dplyr::filter(data, group == "0.01 < RP < 0.1"))[1] != 0) {
        level <- c(level,"0.01 < RP < 0.1")
        col <- c(col,"#00BFC4")
      }
      if(dim(dplyr::filter(data, group == "0.1 < RP < 1"))[1] != 0) {
        level <- c(level,"0.1 < RP < 1")
        col <- c(col,"lightgreen")
      }
      if(dim(dplyr::filter(data, group == "1 < RP < 2"))[1] != 0) {
        level <- c(level,"1 < RP < 2")
        col <- c(col,"#F8766D")
      }
      if(dim(dplyr::filter(data, group == "2 < RP"))[1] != 0) {
        level <- c(level,"2 < RP")
        col <- c(col,"red")
      }
      data$group <- factor(data$group,levels=level,ordered=TRUE)
      data$intersection <- as.factor(name)
      df <- rbind(df,data)
      if(input$RNAseq_data_type_enrich != "List") {
        genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,RNAlog2FC = -1*data$log2FoldChange,
                                       sumRP = data$sumRP, gene_id = data$gene_id)
      }else{
        genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,
                                       sumRP = data$sumRP, gene_id = data$gene_id)
      }
      bed_list1 <- list()
      for(name2 in unique(data$group)){
        if(name2 != "Others"){
          table <- genelist[[name]] %>% dplyr::filter(group == name2)
          print(head(table))
          gene <- table$gene_id
          up_peak3 <- NULL
          if(length(gene) != 0){
            up_peak <- subset(pre_mmAnno_enrich()[[name]], gene_id %in% gene)
            up_peak2 <- as.data.frame(up_peak)
            up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
            up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
            up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
            mcols(up_peak3) <- DataFrame(Group = paste0(name,":",name2))
            print(up_peak3)
            bed_list1[[name2]] <- up_peak3
          }
        }
      }
      bed_list[[name]] <- bed_list1
    }
    data <- df 
    for(name in enrich_select_RNA_debounce()){
      check <- data %>% dplyr::filter(intersection == name) %>% 
        dplyr::filter(group != "Others") %>% summarise(n())
      if(check <= 1) data <- data %>% dplyr::filter(intersection != name)
    }
    data$intersection <- gsub("-","-\n",data$intersection)
    collist <- unique(data$group)
    
    col <-c("gray","blue","#00BFC4","lightgreen","#F8766D","red")
    if(input$RNAseq_data_type_enrich != "List"){
      data$log10FoldChange <- log10(2^data$log2FoldChange)
      data$log10FoldChange <- as.numeric(data$log10FoldChange)
      if(dim(data)[1] == 0) validate("boxplot: There are few genes with |RP| > 1")
      if (length(collist) >= 3){
        stat.test <- data %>% dplyr::group_by(intersection) %>% 
          tukey_hsd(log10FoldChange ~ group)
        stat.test <- stat.test %>% add_significance("p.adj")
        stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      }else{
        group1 <- dplyr::filter(data, group == collist[1])
        group2 <- dplyr::filter(data, group == collist[2])
        if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
          stat.test <- data %>% dplyr::group_by(intersection) %>% 
            t_test(log10FoldChange ~ group)
          stat.test <- stat.test %>% add_significance()
          stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
        }else stat.test <- NULL
      }
      p <- try(ggpubr::ggboxplot(data, x = "group", y = "log10FoldChange",
                                 fill = "group", scales = "free", 
                                 xlab = FALSE, ylab = paste0("RNAseq log10(",cond2,"/",cond1,")"))+theme_bw(base_size = 15)+
                 xlab(NULL)+scale_fill_manual(values = col) + scale_x_discrete(labels = label_wrap_gen(8)))
      if(input$Statistics_enrich !="not selected") p <- p + stat_pvalue_manual(stat.test,hide.ns = T, size = 5)         
      p <- facet(p, facet.by = "intersection",
                 panel.labs.background = list(fill = "transparent", color = "transparent"),
                 scales = "free", short.panel.labs = T, panel.labs.font = list(size=15))
      stat.test <- stat.test[,-2]
      stat.test <- stat.test[,1:9]
      colnames(stat.test)[1] <- "Peak"
      df <- list()
      df[["plot"]] <- p
      df[["statistical_test"]] <- stat.test
    }else {
      df <- list()
      table2 <- data %>% group_by(group, intersection) %>%
        summarise(count = n(), .groups = "drop")
      p <- ggplot(table2,aes(x = group,y= count, fill= group)) +
        geom_col(position=position_dodge2(preserve = "single"))+scale_fill_manual(values = col) +
        theme_bw(base_size = 15)+facet_wrap(~intersection,scales = "free") +
        xlab(NULL) + ylab("Gene number") + scale_x_discrete(labels = label_wrap_gen(8))
      df[["plot"]] <- p
    }
    df[["genelist"]] <- genelist
    df[["bedlist"]] <- bed_list
    print(bed_list)
    
    return(df)
  })
  
  
  RNAseq_popu_enrich <- reactive({
    if(!is.null(pre_RP_all_table_enrich()) && !is.null(RNAseq_name_enrich())){
      withProgress(message = "Preparing barplot",{
        up_name <- RNAseq_name_enrich()[2]
        down_name <- RNAseq_name_enrich()[1]
        df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        for(name in enrich_select_RNA_debounce()){
          epi_name <- paste0(name, "_associated")
          table <- pre_RP_all_table_enrich()[[name]] %>% dplyr::mutate(
            type =if_else(withPeakN > 0, epi_name, "not_associated")
          )
          table$Group <- gsub(".+\\:","",table$Group)
          table$Group <- paste0(table$Group," genes")
          table$type <- factor(table$type,levels = c(epi_name,"not_associated"))
          table$intersection <- name
          df <- rbind(df, table)
        }
        table <- df
        table$intersection <- gsub("-","-\n",table$intersection)
        table2 <- table %>% group_by(Group, intersection, withPeakN) %>%
          summarise(count = n(), .groups = "drop")
        p <- ggplot(table2,aes(x = withPeakN,y= count,fill = intersection)) +
          geom_col(position=position_dodge2(preserve = "single")) +
          theme_bw(base_size = 15)+facet_wrap(~Group,scales = "free",ncol = 3) +
          xlab("Number of associated peaks")+guides(fill=guide_legend(title="associated_\npeak_type"))
        incProgress(1)
      })
      return(p)
    }
  })
  output$int_box_enrich <- renderPlot({
    withProgress(message = "Boxplot",{
      if(!is.null(enrich_select_RNA_debounce()) && !is.null(pre_RP_all_table_enrich()) && 
         input$Species_enrich != "not selected" && !is.null(pre_mmAnno_enrich()) &&
         !is.null(RNAseqDEG_anno_enrich())){
        if(input$RNAseq_data_type_enrich == "List") validate("Boxplot is not available for this input type. Please use 'DEG_result' or 'Raw_count' data.")
        RNAseq_boxplot_enrich()[["plot"]]
      }
    })
  })
  output$int_boxplot_table_enrich <- DT::renderDataTable({
    if(!is.null(enrich_select_RNA_debounce()) && !is.null(pre_RP_all_table_enrich()) && 
       input$Species_enrich != "not selected" && !is.null(pre_mmAnno_enrich()) &&
       !is.null(RNAseqDEG_anno_enrich())){
      if(input$Statistics_enrich !="not selected"){
        if(input$RNAseq_data_type_enrich == "List") validate("Boxplot is not available for this input type. Please use 'DEG_result' or 'Raw_count' data.")
        RNAseq_boxplot_enrich()[["statistical_test"]]
      }
    }
  })
  output$bar_rna_enrich <- renderPlot({
    if(!is.null(enrich_select_RNA_debounce()) && !is.null(pre_RP_all_table_enrich()) && 
       !is.na(input$DEG_fdr_enrich) && !is.na(input$DEG_fc_enrich) && !is.null(RNAseq_popu_enrich()) &&
       input$Species_enrich != "not selected" && !is.null(pre_mmAnno_enrich()) &&
       !is.null(RNAseqDEG_anno_enrich())){
      RNAseq_popu_enrich()
    }
  })
  output$bar_chip_enrich <- renderPlot({
    if(!is.null(enrich_select_RNA_debounce()) && !is.null(pre_RP_all_table_enrich()) && 
       !is.na(input$DEG_fdr_enrich) && !is.na(input$DEG_fc_enrich) && !is.null(ChIPseq_popu_enrich()) &&
       input$Species_enrich != "not selected" && !is.null(pre_mmAnno_enrich()) &&
       !is.null(RNAseqDEG_anno_enrich())){
      ChIPseq_popu_enrich()
    }
  })
  ChIPseq_popu_enrich <- reactive({
    RNA <- RNAseqDEG_anno_enrich()
    if(!is.null(RNA) && !is.null(pre_mmAnno_enrich())){
      df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
      for(name in enrich_select_RNA_debounce()){
        mmano <- as.data.frame(pre_mmAnno_enrich()[[name]])
        mmano$locus <- paste0(mmano$seqnames,":",mmano$start,"-",mmano$end)
        merge <- merge(mmano,RNA, by = "gene_id")
        if(input$RNAseq_data_type_enrich != "List"){
          up_name <- paste0(RNAseq_name_enrich()[2]," gene only")
          down_name <- paste0(RNAseq_name_enrich()[1]," gene only")
          merge <- merge  %>% dplyr::mutate(Up_gene = if_else(padj < input$DEG_fdr_enrich & log2FoldChange > log(input$DEG_fc_enrich,2), 1, 0),
                                            Down_gene = if_else(padj < input$DEG_fdr_enrich & log2FoldChange < -log(input$DEG_fc_enrich,2), 1, 0),
                                            NS_gene = if_else(padj >= input$DEG_fdr_enrich | baseMean == 0 | is.na(padj) |
                                                                (padj <= input$DEG_fdr_enrich & abs(log2FoldChange) <= log(input$DEG_fc_enrich,2)), 1, 0))
          merge[is.na(merge)] <- 0
          merge2 <- merge %>% group_by(locus,Group) %>% summarise(Total_associated_gene = sum(Up_gene)+sum(Down_gene)+sum(NS_gene),Group = Group,
                                                                  Up_gene = sum(Up_gene), Down_gene = sum(Down_gene), NS_gene = sum(NS_gene))
          table <- merge2 %>% dplyr::mutate(type = if_else(Up_gene > 0 & Down_gene == 0 & NS_gene == 0, up_name,
                                                           if_else(Up_gene == 0 & Down_gene > 0 & NS_gene == 0, down_name,
                                                                   if_else(Up_gene == 0 & Down_gene == 0 & NS_gene > 0, "NS gene only", "Multiple type of genes"))))
          table$type <- factor(table$type,levels = c(down_name,up_name,"NS gene only","Multiple type of genes"))
        }else{
          merge2 <- merge %>% group_by(locus, Group.y) %>% summarise(Count = n()) 
          merge2 <- spread(merge2,Group.y,Count)
          merge2[is.na(merge2)] <- 0
          namelist <- unique(merge$Group.y)
          merge2$Total_associated_gene <- rowSums(merge2[,-1],na.rm = TRUE)
          vec <- c()
          merge3 <- merge2[,-1]
          merge3 <- merge3[,-length(colnames(merge3))]
          for(i in 1:dim(merge3)[1]){
            type <- merge3[i,]
            type <- paste(names(type[which(type != 0)]),collapse = "&")
            vec <- c(vec,type) 
          }
          merge2$type <- vec
          table <- merge2
        }
        table$Group <- paste0(table$Group," peak")
        if(max(table$Total_associated_gene) < 5) table$Total_associated_gene <- as.character(table$Total_associated_gene) else table$Total_associated_gene <- as.numeric(table$Total_associated_gene)
        table$RNA_group <- name
        df <- rbind(df, table)
      }
      table <- df
      table$RNA_group <- gsub("-","\n",table$RNA_group)
      
      p2 <- ggplot(table,aes(x = Total_associated_gene, fill = type)) + geom_bar(position = "stack") +
        theme_bw(base_size = 15)+facet_wrap(~RNA_group,scales = "free") + ylab("Number of peaks") +
        xlab("Number of associated genes")+guides(fill=guide_legend(title="associated_gene_type"))
      if(input$RNAseq_data_type_enrich != "List"){
        col <- c("#00BFC4","#F8766D","grey","black")
        p2 <- p2 +
          scale_fill_manual(values = col) 
      }
      return(p2)
    }
  })
  pre_RP_all_table_enrich <- reactive({
    if(!is.null(RP_enrich()) && !is.null(regulatory_potential_enrich())){
      table_list <- list()
      for(name in names(RP_enrich())){
        target_result <- regulatory_potential_enrich()[[name]]$data
        target_result$epigenome_category <- name
        table <- NULL
        if(str_detect(target_result$gene_id[1], "FBgn")){
          symbol <- id_convert(my.symbols = target_result$gene_id,Species = input$Species_enrich,type = "ENSEMBL")
          symbol <- symbol %>% distinct(ENSEMBL, .keep_all = TRUE)
          symbol <- symbol$SYMBOL
        }else symbol <- target_result$Symbol
        if(!is.null(pre_mmAnno_enrich()[[name]])) {
          if(input$RNAseq_data_type_enrich != "List"){
            table <- data.frame(Symbol = symbol,
                                Group = paste0(target_result$epigenome_category,":",target_result$gene_category),
                                RNA_log2FC = -target_result$log2FoldChange,
                                RNA_padj = target_result$padj,
                                regulatory_potential = target_result$sumRP,
                                withPeakN = target_result$withPeakN,
                                gene_id = target_result$gene_id)
          }else{
            my.symbols <- target_result$gene_id
            gene_IDs <- id_convert(my.symbols,Species = input$Species_enrich,type = "ENTREZID")
            colnames(gene_IDs) <- c("gene_id","SYMBOL")
            gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T)
            target_result <- merge(target_result, gene_IDs, by="gene_id")
            table <- data.frame(Symbol = target_result$SYMBOL,
                                Group = paste0(target_result$epigenome_category,":",target_result$Group),
                                regulatory_potential = target_result$sumRP,
                                withPeakN = target_result$withPeakN,
                                gene_id = target_result$gene_id)
          }
          table_list[[name]] <- table
        }
      }
      return(table_list)
    }
  })
  RP_all_table_enrich <- reactive({
    if(!is.null(pre_RP_all_table_enrich())){
      df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
      table <- pre_RP_all_table_enrich()
      for(name in names(table)){
        df <- rbind(df, table[[name]])
      }
      return(df)
    }
  })
  
  RP_selected_table_enrich <- reactive({
    table <- RP_all_table_enrich() %>% dplyr::filter(Group == input$RNAseqGroup_enrich)
    return(table)
  })
  
  output$RP_table_enrich <- renderDT({
    if(!is.null(input$RNAseqGroup_enrich) && !is.null(enrich_select_RNA_debounce()) &&
       !is.null(input$peak_distance_enrich && !is.null(pre_mmAnno_enrich())) &&
       !is.null(regulatory_potential_enrich())){
      RP_selected_table_enrich() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  output$download_enrichintbar = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        gridExtra::grid.arrange(RNAseq_popu_enrich(), ChIPseq_popu_enrich(), ncol = 1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_enrichintbox = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), "%Y%m%d_"),"boxplot",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- pdf_h(input$enrich_select_RNA)
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- pdf_w(input$enrich_select_RNA)
        }else pdf_width <- input$enrich_pdf_width
        fs <- c()
        dir.create("boxplot/",showWarnings = FALSE)
        boxplot <- paste0("boxplot/","boxplot.pdf")
        fs <- c(fs, boxplot)
        pdf(boxplot, height = pdf_height, width = pdf_width)
        print(RNAseq_boxplot_enrich()[["plot"]])
        dev.off()
        if(input$RNAseq_data_type_enrich != "List"){
          boxplot_table <- paste0("boxplot/","boxplot.txt")
          fs <- c(fs,boxplot_table)
          write.table(RNAseq_boxplot_enrich()[["statistical_test"]],boxplot_table,col.names = T,row.names = F,sep = "\t",quote = F)
        }
        dir.create("boxplot/gene_list/",showWarnings = FALSE)
        dir.create("boxplot/bed/",showWarnings = FALSE)
        for(name in enrich_select_RNA_debounce()){
          genelist <- paste0("boxplot/gene_list/",name,".txt")
          up <- RNAseq_boxplot_enrich()[["genelist"]][[name]]
          fs <- c(fs,genelist)
          write.table(up,genelist,col.names = T,row.names = F,sep = "\t",quote = F)
        }
        for(name in names(RNAseq_boxplot_enrich()[["bedlist"]])){
          dir.create(paste0("boxplot/bed/",name),showWarnings = FALSE)
          for(file in names(RNAseq_boxplot_enrich()[["bedlist"]][[name]])){
            bed <- RNAseq_boxplot_enrich()[["bedlist"]][[name]][[file]]
            bed_name <- paste0("boxplot/bed/",name,"/",file,".bed")
            bed_name <- gsub(" < ","....",bed_name)
            fs <- c(fs, bed_name)
            write.table(as.data.frame(bed), bed_name, row.names = F, col.names = F,sep = "\t", quote = F)
          }
        }
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    },
    contentType = "application/zip"
  )
  
  
  output$download_enrichKSplot = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), "%Y%m%d_"),"KSplot",".zip")
    },
    content = function(fname) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$enrich_pdf_width
        fs <- c()
        dir.create("KSplot/",showWarnings = FALSE)
        for(name in names(regulatory_potential_enrich())){
          ksplot <- paste0("KSplot/",name,".pdf")
          fs <- c(fs,ksplot)
          pdf(ksplot, height = pdf_height, width = pdf_width)
          print(regulatory_potential_enrich()[[name]])
          dev.off()
        }
        kstable <- paste0("KSplot/ks_test.txt")
        fs <- c(fs,kstable)
        write.table(ks_tables_enrich(), kstable, row.names = F, sep = "\t", quote = F)
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    }
  )
  output$download_RP_enrich_table = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), "%Y%m%d_"),"RP_table",".zip")
    },
    content = function(fname) {
      withProgress(message = "Preparing download",{
        fs <- c()
        dir.create("RP_table/",showWarnings = FALSE)
        dir.create("RP_table/selected_table(epigenome--RNA)/",showWarnings = FALSE)
        dir.create("RP_table/selected_bed(epigenome--RNA)/",showWarnings = FALSE)
        RP_summary <- paste0("RP_table/summary.txt")
        fs <- c(fs,RP_summary)
        write.table(RP_all_table_enrich(), RP_summary, row.names = F, sep = "\t", quote = F)
        for(name in unique(RP_all_table_enrich()$Group)){
          RP_selected <- paste0("RP_table/selected_table(epigenome--RNA)/",name,".txt")
          RP_selected_bed <- paste0("RP_table/selected_bed(epigenome--RNA)/",name,".bed")
          RP_selected <- gsub(":","--",RP_selected)
          RP_selected_bed <- gsub(":","--",RP_selected_bed)
          fs <- c(fs,RP_selected,RP_selected_bed)
          table <- RP_all_table_enrich() %>% dplyr::filter(Group == name)
          write.table(table, RP_selected, row.names = F, sep = "\t", quote = F)
          
          gene <- table$gene_id
          peak <- gsub("\\:.+$", "", name)
          y <- NULL
          if(!is.null(pre_mmAnno_enrich())) {
            up_peak <- subset(pre_mmAnno_enrich()[[peak]], gene_id %in% gene)
            up_peak2 <- as.data.frame(up_peak)
            up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
            up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
            up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
            mcols(up_peak3) <- DataFrame(Group = "up")
            y <- as.data.frame(up_peak3)
          }
          write.table(y, RP_selected_bed, row.names = F, col.names = F,sep = "\t", quote = F)
        }
        zip(zipfile=fname, files=fs)
        incProgress(1)
      })
    }
  )
  
  observeEvent(input$enrich_DEG_result, ({
    updateCollapse(session,id =  "input_collapse_enrich_RP", open="KS_panel")
  }))
  #enrich Integrative trackplot-------
  output$file_enrich1 <- renderUI({
    if(length(list.files("./Volume/")) > 0){
      list <- list.files("Volume",full.names = T,recursive=T)
      bw <- c(str_subset(list,".bw"),str_subset(list,".BigWig"))
      selectInput("file_enrich1",label=NULL,choices = bw,multiple=T)
    }else{
      fileInput("file_enrich1",
                NULL,
                accept = c("bw","BigWig"),
                multiple = TRUE,
                width = "80%")
    }
  })
  pre2_bws_enrich <- reactive({
    if(is.null(input$file_enrich1)){
      if(input$goButton_enrich > 0 ){
        df<-list()
        df[["A_1.bw"]] <- "data/bigwig/A_1.BigWig"
        df[["A_2.bw"]] <- "data/bigwig/A_2.BigWig"
        df[["B_1.bw"]] <- "data/bigwig/B_1.BigWig"
        df[["B_2.bw"]] <- "data/bigwig/B_2.BigWig"
        return(df)
      }
      return(NULL)
    }else{
      if(length(list.files("./Volume/")) > 0){
        files <- input$file_enrich1
        names(files) <- bigwig_name(input$file_enrich1)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file_enrich1[, 1])){
          file <- input$file_enrich1[[nr, 'datapath']]
          name <- c(name, bigwig_name(input$file_enrich1[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
      }
      return(files)
    }
  })
  pre_bws_enrich <- debounce(pre2_bws_enrich, 1000)
  bws_enrich <- reactive({
    return(bws_ordering(bws=pre_bws_enrich(),sample_order=input$sample_order_enrich,additional=FALSE))
  })
  observeEvent(input$RP_table_enrich_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_enrich_RP", open="int_Trackplot_panel")
  }))
  gene_position_enrich <- reactive({
    if(input$Species_enrich != "not selected"){
      return(genes(txdb_enrich()))
    }
  })
  output$int_igv_uprange_enrich <- renderUI({
    if(!is.null(int_goi_gene_position_enrich()) && !is.null(int_goi_promoter_position_enrich())){
      y <- int_goi_promoter_position_enrich()
      gene_position <- int_goi_gene_position_enrich()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("int_igv_uprange_enrich","Range (x-axis):",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim_enrich <- renderUI({
    numericInput("int_igv_ylim_enrich","Max peak intensity (y-axis):", value = 2, min = 0)
  })
  
  int_goi_promoter_position_enrich<- reactive({
    if(!is.null(input$RP_table_enrich_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table_enrich()[input$RP_table_enrich_rows_selected,]$gene_id
      y <- NULL
      if(!is.null(pre_mmAnno_enrich())) {
        up_peak <- subset(mmAnno_enrich(), gene_id %in% gene)
        up_peak2 <- as.data.frame(up_peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      return(y)
    }
  })
  
  int_goi_gene_position_enrich <- reactive({
    if(!is.null(input$RP_table_enrich_rows_selected)){
      gene <- RP_selected_table_enrich()[input$RP_table_enrich_rows_selected,]$gene_id
      gene_position <- as.data.frame(subset(gene_position_enrich(), gene_id %in% gene))
      return(gene_position)
    }
  })
  
  ref_enrich <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species_enrich))
    return(ref)
  })
  int_data_track_enrich <- reactive({
    if(!is.null(input$RP_table_enrich_rows_selected)){
      return(data_trac(y=int_goi_promoter_position_enrich(),gene_position=int_goi_gene_position_enrich(),
                       gen=ref_enrich(),txdb=txdb_enrich(),org=org_enrich(),
                       bw_files=bws_enrich(),
                       track_additional_files=NULL))
    }
  })
  
  int_highlight_trackplot_enrich <- reactive({
    if(!is.null(input$RP_table_enrich_rows_selected) && !is.null(input$int_igv_uprange_enrich)){
      library(Gviz)
      y <- int_goi_promoter_position_enrich()
      gene_position <- int_goi_gene_position_enrich()
      chr <- gene_position$seqnames
      df <- int_data_track_enrich()
      start <-c()
      width <- c()
      col <- c()
      fill <- c()
      for(i in 1:length(rownames(y))){
        if(y[i,]$start < input$int_igv_uprange_enrich[2] && y[i,]$end > input$int_igv_uprange_enrich[1]){
          start <- c(start, y[i,]$start)
          width <- c(width, y[i,]$width)
          if(y[i,]$Group == "up") {
            col <- c(col,"#E41A1C")
            fill <- c(fill,"#FBB4AE")
          }
          if(y[i,]$Group == "down") {
            col <- c(col,"#377EB8")
            fill <- c(fill,"#B3CDE3")
          }
        }
      }
      if(length(start) != 0){
        ht <- HighlightTrack(trackList = df,
                             start = start, width = width,col = col,fill = fill,
                             chromosome = chr, alpha=0.5) 
      }else ht <- NULL
      return(ht)
    }
  })
  gtrack_enrich <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=0.8)
      return(gtrack)
    })
  })
  int_goi_trackplot_enrich <- reactive({
    if(!is.null(input$RP_table_enrich_rows_selected) &&
       !is.null(int_goi_promoter_position_enrich()) && 
       !is.null(int_goi_gene_position_enrich()) && 
       !is.null(gtrack_enrich()) &&
       !is.null(input$int_igv_uprange_enrich)){
      library(Gviz)
      if(!is.null(int_highlight_trackplot_enrich())){
        plot<- plotTracks(list(gtrack_enrich(), int_highlight_trackplot_enrich()),
                          from = input$int_igv_uprange_enrich[1], 
                          to = input$int_igv_uprange_enrich[2],ylim=c(0,input$int_igv_ylim_enrich),
                          type="hist")
      }else{
        df <- int_data_track_enrich()
        df[["gtrack"]] <- gtrack_enrich()
        plot<- plotTracks(df,
                          from = input$int_igv_uprange_enrich[1], 
                          to = input$int_igv_uprange_enrich[2],ylim=c(0,input$int_igv_ylim_enrich),
                          type="hist")
      }
      return(plot)
    }
  })
  output$int_trackplot_goi_enrich <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$RP_table_enrich_rows_selected)){
        int_goi_trackplot_enrich()
      }
    })
  })
  output$download_enrich_int_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(int_highlight_trackplot_enrich())){
          plot<- plotTracks(list(gtrack_enrich(), int_highlight_trackplot_enrich()),
                            from = input$int_igv_uprange_enrich[1], 
                            to = input$int_igv_uprange_enrich[2],ylim=c(0,input$int_igv_ylim_enrich),
                            type="hist")
        }else{
          df <- int_data_track_enrich()
          df[["gtrack"]] <- gtrack_enrich()
          plot<- plotTracks(df,
                            from = input$int_igv_uprange_enrich[1], 
                            to = input$int_igv_uprange_enrich[2],ylim=c(0,input$int_igv_ylim_enrich),
                            type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  #enrich integrative functional enrichment-------
  int_Hallmark_set_enrich <- reactive({
    if(!is.null(input$intGeneset_enrich)){
      return(GeneList_for_enrichment(Species = input$Species_enrich, Gene_set = input$intGeneset_enrich, org = org_enrich()))
    }
  })
  
  order_for_intGroup_enrich <- reactive({
    order <- c("0 < RP < 0.01","0.01 < RP < 0.1","0.1 < RP < 1","1 < RP < 2","2 < RP")
    order <- factor(order,levels = order,ordered = TRUE)
    return(order)
  })
  output$intGroup_for_RPstatus_enrich <- renderUI({
    if(input$withRNAseq_enrich_enrich_type == "boxplot"){
      selectInput("intGroup_for_RPstatus_enrich","Peak",names(RNAseq_boxplot_enrich()[["genelist"]]),multiple = F)
    }
  })
  
  output$intGroup_enrich <- renderUI({
    if(!is.null(RP_all_table_enrich())){
      if(input$withRNAseq_enrich_enrich_type == "boxplot"){
        if(input$RNAseq_data_type_enrich == "List") validate("When using gene lists as input for RNA-seq data, RP-status-based enrichment analysis cannot be applied. \nPlease utilize the 'Relationship (Epigenome:RNAseq)' mode instead.")
        selectInput("intGroup_enrich","Group",
                    order_for_intGroup_enrich(),selected=order_for_intGroup_enrich(),multiple = T)
      }else{
        selectInput("intGroup_enrich","Group (intersection:RNAseq)",
                    unique(RP_all_table_enrich()$Group),multiple = T)
      }
    }
  })
  
  pre_intGroup_enrich_debounce <- reactive({
    if(!is.null(input$intGroup_enrich)) return(input$intGroup_enrich) else return(NULL)
  })
  intGroup_enrich_debounce <- debounce(pre_intGroup_enrich_debounce, 1000)
  
  output$intGeneset_enrich <- renderUI({
    selectInput('intGeneset_enrich', 'Gene Set', gene_set_list)
  })
  
  withRNAseq_enrichment_analysit_genelist_enrich <- reactive({
    if(input$withRNAseq_enrich_enrich_type == "boxplot"){
      df <- RNAseq_boxplot_enrich()[["genelist"]][[input$intGroup_for_RPstatus_enrich]]
      colnames(df)[2] <- "Group"
      df <- df %>% dplyr::filter(Group != "Others")
      return(df)
    }
  })
  
  selected_int_group_enrich <- reactive({
    group <- intGroup_enrich_debounce()
    df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
    colnames(df) <- c("ENTREZID","Group")
    for(name in group){
      if(input$withRNAseq_enrich_enrich_type == "boxplot"){
        table <- withRNAseq_enrichment_analysit_genelist_enrich()
      }else{
        table <- RP_all_table_enrich()
      }
      table <- table %>% dplyr::filter(Group == name)
      if(dim(table)[1] != 0){
        entrezid <- table$gene_id
        if(str_detect(table$gene_id[1], "FBgn")){
          my.symbols <- gsub("\\..*","", table$gene_id)
          gene_IDs<-AnnotationDbi::select(org(input$Species_enrich),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","ENTREZID"))
          colnames(gene_IDs) <- c("gene_id","ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(ENTREZID, .keep_all = T)
          gene_IDs <- na.omit(gene_IDs)
          table <- merge(table,gene_IDs,by="gene_id")
          entrezid <- table$ENTREZID
        }
        df2 <- data.frame(ENTREZID = entrezid, Group = table$Group)
        df <- rbind(df,df2)
      }
    }
    return(df)
  })
  int_enrich_enrich <- reactive({
    return(enrich_viewer_forMulti2(data3 = selected_int_group_enrich(), Species = input$Species_enrich, org = org_enrich(),
                                   H_t2g = int_Hallmark_set_enrich(),Gene_set = input$intGeneset_enrich))
  })
  int_enrich_list_enrich <- reactive({
    return(enrich_gene_list(data = selected_int_group_enrich(),
                            Gene_set = input$intGeneset_enrich, org = org_enrich(), H_t2g = int_Hallmark_set_enrich()))
  })
  int_enrich_plot_enrich <- reactive({
    return(enrich_genelist(data = selected_int_group_enrich(),
                           enrich_gene_list = int_enrich_list_enrich(),type ="withRNAseq"))
  })
  output$int_enrichment1_enrich <- renderPlot({
    dotplot_for_output(data = int_enrich_enrich(),
                       plot_genelist = int_enrich_plot_enrich(), Gene_set = input$intGeneset_enrich, 
                       Species = input$Species_enrich)
  })
  int_enrich_table_enrich <- reactive({
    return(enrich_for_table(data = as.data.frame(int_enrich_enrich()), H_t2g = int_Hallmark_set_enrich(), Gene_set = input$intGeneset_enrich))
  })
  output$int_enrichment_result_enrich <- DT::renderDataTable({
    int_enrich_table_enrich()
  })
  output$download_enrich_int_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$intGeneset_enrich,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = int_enrich_enrich(),
                           plot_genelist = int_enrich_plot_enrich(), Gene_set = input$intGeneset_enrich, 
                           Species = input$Species_enrich)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_enrich_int_enrichment_enrich_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$intGeneset_enrich,".txt")
    },
    content = function(file){write.table(int_enrich_table_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  #Docker only---------
  ##Motif enrichment analysis----------
  output$reference_download <- renderText({
    if(input$Species != "not selected"){
      if(input$HOMERreference_start > 0){
        withProgress(message = paste0("Install Homer reference for ", gsub(".+\\(","",gsub(")", "", input$Species)),". It takes about > 10 min."),{
          if(updateCounter_homerRef$i > 0) log <- paste0(robust.system(paste0('perl /usr/local/homer/.//configureHomer.pl -install ',gsub(".+\\(","",gsub(")", "", input$Species))))$stderr,collapse = "\n")
          paste0(gsub(".+\\(","",gsub(")", "", input$Species)), " installation has been completed.")
        })
      }else{
        paste("If the HOMER reference for your dataset species has not been installed in this container,",
              "you need to press the 'Install HOMER reference' button", sep = "\n")
      }
    }
  })
  updateCounter_homerRef <- reactiveValues(i = 0)
  
  observe({
    input$HOMERreference_start
    isolate({
      updateCounter_homerRef$i <- updateCounter_homerRef$i + 1
    })
  })
  
  
  #Restart
  observeEvent(input$reference_download, {
    isolate(updateCounter_homerRef$i == 0)
    updateCounter_homerRef <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species, {
    isolate(updateCounter_homerRef$i == 0)
    updateCounter_homerRef <<- reactiveValues(i = 0)
  }) 
  
  output$homer_bg <- renderUI({
    if(input$Genomic_region == "Genome-wide"){
      radioButtons('homer_bg','Background sequence',
                   c('random'="random",
                     'peak call files'="peakcalling"
                   ),selected = "random")
    }
  })
  
  updateCounter <- reactiveValues(i = 0)
  
  observe({
    input$motifButton
    isolate({
      updateCounter$i <- updateCounter$i + 1
    })
  })
  
  
  #Restart
  observeEvent(preMotif_list(), {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Genomic_region, {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species, {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_unknown, {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_size, {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_size2, {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_bg, {
    isolate(updateCounter$i == 0)
    updateCounter <<- reactiveValues(i = 0)
  }) 
  
  output$homer_size2 <- renderUI({
    if(!is.null(input$homer_size)){
      if(input$homer_size == "custom"){
        numericInput('homer_size2','Size of the region for motif finding',value=200, step=100)}}
  })
  
  
  preMotif_list <- reactive({
    if(input$Species != "not selected"){
      df <- list()
      collist <- collist_bw_pair()
      df[[paste0(unique(collist)[2],"_high")]] <- data_degcount_up()
      df[[paste0(unique(collist)[1],"_high")]] <- data_degcount_down()
      return(df)
    }
  })
  
  output$homer_sample_order <- renderUI({
    if(input$motifButton > 0 && !is.null(enrich_motif()) && 
       !is.null(input$homer_unknown) && input$Species != "not selected"){
      selectInput("homer_sample_order", "Order of groups", names(preMotif_list()),
                  selected=names(preMotif_list()),multiple = TRUE)
    }
  })
  enrich_motif <- reactive({
    if(updateCounter$i > 0 && input$motifButton > 0 && !is.null(preMotif_list()) 
       && input$Species != "not selected" && !is.null(input$homer_unknown)){
      if(input$homer_size == "given") size <- "given"
      if(input$homer_size == "custom") size <- input$homer_size2
      if(input$Genomic_region == "Genome-wide"){
        return(findMotif(df= preMotif_list(), anno_data = deg_result_anno2(),back = input$homer_bg,motif_length=input$homer_length,
                         Species = input$Species, motif=input$homer_unknown,size=size,bw_count=bw_count()))
      }else return(findMotif(df= data_degcount2(), anno_data = promoter_region(),
                             Species = input$Species,motif_length=input$homer_length,
                             type = "Promoter", motif=input$homer_unknown,size=size))
    }else return(NULL)
  })
  
  output$motif_plot <- renderPlot({
    if(input$motifButton > 0 && !is.null(enrich_motif()) && 
       !is.null(input$homer_unknown) && input$Species != "not selected"){
      homer_Motifplot(df = enrich_motif(),showCategory = input$homer_showCategory,group_order=input$homer_sample_order)
    }
  })
  output$download_motif_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- homer_Motifplot(df = enrich_motif(),showCategory = input$homer_showCategory)
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_homer_report = downloadHandler(
    filename = function() {
      paste0("HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- enrich_motif()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  
  output$download_motif_table = downloadHandler(
    filename = function() {
      paste0("known_motif_table",".txt")
    },
    content = function(file){write.table(motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  motif_table <- reactive({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      return(known_motif(enrich_motif()))
    }
  })
  
  output$motif_result <- DT::renderDT({
    if(input$motifButton > 0 && !is.null(enrich_motif()) && input$Species != "not selected"){
      motif_table()
    }
  })
  # with RNAseq HOMER-----------
  output$with_reference_download <- renderText({
    if(input$Species != "not selected"){
      if(input$with_HOMERreference_start > 0){
        withProgress(message = paste0("Install Homer reference for ", gsub(".+\\(","",gsub(")", "", input$Species)),". It takes about > 10 min."),{
          if(with_updateCounter_homerRef$i > 0) log <- paste0(robust.system(paste0('perl /usr/local/homer/.//configureHomer.pl -install ',gsub(".+\\(","",gsub(")", "", input$Species))))$stderr,collapse = "\n")
          paste0(gsub(".+\\(","",gsub(")", "", input$Species)), " installation has been completed.")
        })
      }else {
        paste("If the HOMER reference for your dataset species has not been installed in this container,",
              "you need to press the 'Install HOMER reference' button", sep = "\n")
      }
    }
  })
  with_updateCounter_homerRef <- reactiveValues(i = 0)
  
  observe({
    input$with_HOMERreference_start
    isolate({
      with_updateCounter_homerRef$i <- with_updateCounter_homerRef$i + 1
    })
  })
  
  
  #Restart
  observeEvent(input$with_reference_download, {
    isolate(with_updateCounter_homerRef$i == 0)
    with_updateCounter_homerRef <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species, {
    isolate(with_updateCounter_homerRef$i == 0)
    with_updateCounter_homerRef <<- reactiveValues(i = 0)
  }) 
  output$with_homer_bg <- renderUI({
    if(input$Genomic_region == "Genome-wide"){
      radioButtons('with_homer_bg','Background sequence',
                   c('random'="random",
                     'peak call files'="peakcalling"
                   ),selected = "peakcalling")
    }
  })
  
  with_updateCounter <- reactiveValues(i = 0)
  
  observe({
    input$with_motifButton
    isolate({
      with_updateCounter$i <- with_updateCounter$i + 1
    })
  })
  
  
  #Restart
  observeEvent(with_preMotif_list(), {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Group_homer, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$with_homer_unknown, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$with_homer_size, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$with_homer_size2, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$with_homer_bg, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$withRNAseq_pair_homer_type, {
    isolate(with_updateCounter$i == 0)
    with_updateCounter <<- reactiveValues(i = 0)
  }) 
  
  output$with_homer_size2 <- renderUI({
    if(!is.null(input$with_homer_size)){
      if(input$with_homer_size == "custom"){
        numericInput('with_homer_size2','Size of the region for motif finding',value=200, step=100)}}
  })
  
  
  
  withRNAseq_homer_analysit_genelist <- reactive({
    if(input$withRNAseq_pair_homer_type == "boxplot"){
      if(input$RNAseq_mode_option == "fcw_combined" || input$RNAseq_mode_option == "combined"){
        df <- pari_RNAseq_boxplot()[["bedlist"]]
      }else{
        df <- pari_RNAseq_boxplot()[["bedlist"]][[input$homerGroup_for_RPstatus]]
      }
      return(df)
    }
  })
  
  output$homerGroup_for_RPstatus <- renderUI({
    if(input$withRNAseq_pair_homer_type == "boxplot"){
      if(input$RNAseq_mode_option == "fcw_separate" || input$RNAseq_mode_option == "Classical"){
        selectInput("homerGroup_for_RPstatus","Peak",names(pari_RNAseq_boxplot()[["bedlist"]]),multiple = F)
      }
    }
  })
  
  
  output$Group_homer <- renderUI({
    if(!is.null(pre_RP_selected_table())){
      if(input$withRNAseq_pair_homer_type == "boxplot"){
        selectInput("Group_homer","Group",order_for_intGroup(),selected = order_for_intGroup(),multiple = T)
      }else{
        group <- unique(pre_RP_selected_table()$Group)
        selectInput("Group_homer","Group (Epigenome:RNAseq)",group,selected = group,multiple = T)
      }
    }
  })
  with_preMotif_list <- reactive({
    if(!is.null(RNAseqDEG()) && !is.null(input$Group_homer) && !is.null(pari_RNAseq_boxplot()) && !is.null(mmAnno_pair()) &&
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      group <- input$Group_homer
      if(input$withRNAseq_pair_homer_type == "boxplot"){
        list <- withRNAseq_homer_analysit_genelist()[input$Group_homer]
      }else{
        list <- list()
        if(is.null(group)) validate("Select groups of interest.")
        for(name in group){
          table <- pre_RP_selected_table() %>% dplyr::filter(Group == name)
          if(dim(table)[1] != 0){
            gene <- table$gene_id
            if(input$Genomic_region == "Promoter"){
              tss <- promoters(genes(txdb()),upstream = 0,downstream = 1)
              peak <- subset(tss, gene_id %in% gene)
            }else{
              peak_type <- gsub("\\:.+$","", name)
              peak <- subset(mmAnno_pair()[[peak_type]], gene_id %in% gene)
            }
            list[[gsub(":","--",name)]] <- peak
          }
        } 
      }
      return(list)
    }
  })
  
  
  output$with_homer_sample_order <- renderUI({
    if(input$with_motifButton > 0 && !is.null(with_enrich_motif()) && 
       !is.null(input$with_homer_unknown) && input$Species != "not selected"){
      selectInput("with_homer_sample_order", "Order of groups", names(with_preMotif_list()),
                  selected=names(with_preMotif_list()),multiple = TRUE)
    }
  })
  with_enrich_motif <- reactive({
    if(with_updateCounter$i > 0 && input$with_motifButton > 0 && !is.null(with_preMotif_list()) 
       && input$Species != "not selected" && !is.null(input$with_homer_unknown)){
      if(input$with_homer_size == "given") size <- "given"
      if(input$with_homer_size == "custom") size <- input$with_homer_size2
      return(findMotif(df= with_preMotif_list(), other_data = promoter_region(),back = input$with_homer_bg,
                       motif_length=input$with_homer_length,type="Other",section = "withRNAseq",
                       Species = input$Species, motif=input$with_homer_unknown,size=size,bw_count=bw_count()))
    }else return(NULL)
  })
  
  output$with_motif_plot <- renderPlot({
    if(input$with_motifButton > 0 && !is.null(with_enrich_motif()) && 
       !is.null(input$with_homer_unknown) && input$Species != "not selected"){
      homer_Motifplot(df = with_enrich_motif(),showCategory = input$with_homer_showCategory,
                      section = "withRNAseq",group_order = input$with_homer_sample_order)
    }
  })
  output$download_with_motif_plot = downloadHandler(
    filename = function() {
      paste0("withRNAseq_motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- homer_Motifplot(df = with_enrich_motif(),showCategory = input$with_homer_showCategory)
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_with_homer_report = downloadHandler(
    filename = function() {
      paste0("withRNAseq_HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- with_enrich_motif()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  
  output$download_with_motif_table = downloadHandler(
    filename = function() {
      paste0("withRNAseq_known_motif_table",".txt")
    },
    content = function(file){write.table(with_motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  with_motif_table <- reactive({
    if(input$with_motifButton > 0 && !is.null(with_enrich_motif())){
      return(known_motif(with_enrich_motif()))
    }
  })
  
  output$with_motif_result <- DT::renderDT({
    if(input$with_motifButton > 0 && !is.null(with_enrich_motif()) && input$Species != "not selected"){
      with_motif_table()
    }
  })
  ##Venn motif -----------
  output$reference_download_venn <- renderText({
    if(input$Species_venn != "not selected"){
      if(input$HOMERreference_start_venn > 0){
        withProgress(message = paste0("Install Homer reference for ", gsub(".+\\(","",gsub(")", "", input$Species_venn)),". It takes about > 10 min."),{
          if(updateCounter_homerRef_venn$i > 0) log <- paste0(robust.system(paste0('perl /usr/local/homer/.//configureHomer.pl -install ',gsub(".+\\(","",gsub(")", "", input$Species_venn))))$stderr,collapse = "\n")
          paste0(gsub(".+\\(","",gsub(")", "", input$Species_venn)), " installation has been completed.")
        })
      }else{
        paste("If the HOMER reference for your dataset species has not been installed in this container,",
              "you need to press the 'Install HOMER reference' button", sep = "\n")
      }
    }
  })
  updateCounter_homerRef_venn <- reactiveValues(i = 0)
  
  observe({
    input$HOMERreference_start_venn
    isolate({
      updateCounter_homerRef_venn$i <- updateCounter_homerRef_venn$i + 1
    })
  })
  
  
  #Restart
  observeEvent(input$reference_download_venn, {
    isolate(updateCounter_homerRef_venn$i == 0)
    updateCounter_homerRef_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species_venn, {
    isolate(updateCounter_homerRef_venn$i == 0)
    updateCounter_homerRef_venn <<- reactiveValues(i = 0)
  }) 
  output$homer_showCategory_venn <- renderUI({
    sliderInput("homer_showCategory_venn","Most significant motifs", value=5, min=1,max=20)
  })
  output$homer_size_venn <- renderUI({
    radioButtons('homer_size_venn','Type of the region for motif finding',
                 c('given (exact size)'="given",
                   'custom size'="custom"
                 ),selected = "custom")
  })
  output$homer_bg_venn <- renderUI({
    radioButtons('homer_bg_venn','Background sequence',
                 c('random'="random",
                   'bed files'="peakcalling"
                 ),selected = "random")
  })
  output$homer_bg2_venn <- renderUI({
    if(!is.null(input$homer_bg_venn)){
      if(input$homer_bg_venn == "peakcalling"){
        fileInput('homer_bg2_venn',
                  'Select bed files',
                  accept = c("bed","narrowPeak"),
                  multiple = TRUE,
                  width = "80%")
      }}
  })
  updateCounter_venn <- reactiveValues(i = 0)
  
  observe({
    input$motifButton_venn
    isolate({
      updateCounter_venn$i <- updateCounter_venn$i + 1
    })
  })
  
  
  #Restart
  observeEvent(preMotif_list_venn(), {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$venn_whichGroup1, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species_venn, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_unknown_venn, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_size_venn, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_size2_venn, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_bg_venn, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_bg2_venn, {
    isolate(updateCounter_venn$i == 0)
    updateCounter_venn <<- reactiveValues(i = 0)
  }) 
  Venn_peak_call_files_homer <- reactive({
    if(is.null(input$homer_bg2_venn)){
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$homer_bg2_venn[, 1])){
        file <- input$homer_bg2_venn[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$homer_bg2_venn[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      if(length(names(files2)) > 1){
        files2 <- soGGi:::runConsensusRegions(GRangesList(files2), "none")
      }
      return(files2)
    }
  })
  
  output$homer_size2_venn <- renderUI({
    if(!is.null(input$homer_size_venn)){
      if(input$homer_size_venn == "custom"){
        numericInput('homer_size2_venn','Size of the region for motif finding',value=200, step=100)}}
  })
  output$homer_unknown_venn <- renderUI({
    selectInput("homer_unknown_venn","Type of enrichment analysis",c("known motif","known and de novo motifs"), selected = "known motif")
  })
  observeEvent(input$homer_unknown_venn,({
    updateActionButton(session,"motifButton_venn", "Start")
  }))
  
  
  output$venn_whichGroup1 <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("venn_whichGroup1", "Select intersects", choices = c(names(venn_overlap()$peaklist)),multiple = TRUE)
    }
  })
  preMotif_list_venn <- reactive({
    df <- list()
    for(name in input$venn_whichGroup1){
      df[[name]] <- venn_overlap()$peaklist[[name]]
    }
    return(df)
  })
  
  output$homer_sample_order_venn <- renderUI({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && 
       !is.null(input$homer_unknown_venn) && input$Species_venn != "not selected"){
      selectInput("homer_sample_order_venn", "Order of groups", names(preMotif_list_venn()),
                  selected=names(preMotif_list_venn()),multiple = TRUE)
    }
  })
  enrich_motif_venn <- reactive({
    if(updateCounter_venn$i > 0 && input$motifButton_venn > 0 && !is.null(preMotif_list_venn()) 
       && input$Species_venn != "not selected" && !is.null(input$homer_unknown_venn)){
      if(input$homer_size_venn == "given") size <- "given"
      if(input$homer_size_venn == "custom") size <- input$homer_size2_venn
      if(input$homer_bg_venn == "peakcalling" && is.null(input$homer_bg2_venn)){
        return(NULL)
      }else return(findMotif(df= preMotif_list_venn(),Species = input$Species_venn,size=size,back = input$homer_bg_venn,motif_length=input$homer_length_venn,
                             motif=input$homer_unknown_venn, other_data = Venn_peak_call_files_homer(),type="Other",section="venn"))
    }
  })
  venn_motif_plot <- reactive({
    return(homer_Motifplot(df = enrich_motif_venn(),showCategory = input$homer_showCategory_venn,
                           section="venn",group_order = input$homer_sample_order_venn))
  })
  
  output$motif_venn_plot <- renderPlot({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && 
       !is.null(input$homer_unknown_venn) && input$Species_venn != "not selected"){
      venn_motif_plot()
    }
  })
  
  output$download_motif_venn_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(venn_motif_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_motif_venn_table = downloadHandler(
    filename = function() {
      paste0("known_motif_table",".txt")
    },
    content = function(file){write.table(motif_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  motif_table_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      return(known_motif(enrich_motif_venn()))
    }
  })
  denovo_motif_table_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      return(denovo_motif(enrich_motif_venn()))
    }
  })
  
  output$denovo_motif_venn_result <- DT::renderDT({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && input$Species_venn != "not selected"){
      denovo_motif_table_venn()
    }
  })
  
  output$motif_venn_result <- DT::renderDT({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && input$Species_venn != "not selected"){
      motif_table_venn()
    }
  })
  
  
  output$download_denovo_motif_venn_table = downloadHandler(
    filename = function() {
      paste0("denovo_motif_table",".txt")
    },
    content = function(file){write.table(denovo_motif_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_homer_report_venn = downloadHandler(
    filename = function() {
      paste0("HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- enrich_motif_venn()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  ##Enrich motif -----------
  output$reference_download_enrich <- renderText({
    if(input$Species_enrich != "not selected"){
      if(input$HOMERreference_start_enrich > 0){
        withProgress(message = paste0("Install Homer reference for ", gsub(".+\\(","",gsub(")", "", input$Species_enrich)),". It takes about > 10 min."),{
          if(updateCounter_homerRef_enrich$i > 0) log <- paste0(robust.system(paste0('perl /usr/local/homer/.//configureHomer.pl -install ',gsub(".+\\(","",gsub(")", "", input$Species_enrich))))$stderr,collapse = "\n")
          paste0(gsub(".+\\(","",gsub(")", "", input$Species_enrich)), " installation has been completed.")
        })
      }else{
        paste("If the HOMER reference for your dataset species has not been installed in this container,",
              "you need to press the 'Install HOMER reference' button", sep = "\n")
      }
    }
  })
  updateCounter_homerRef_enrich <- reactiveValues(i = 0)
  
  observe({
    input$HOMERreference_start_enrich
    isolate({
      updateCounter_homerRef_enrich$i <- updateCounter_homerRef_enrich$i + 1
    })
  })
  
  
  #Restart
  observeEvent(input$reference_download_enrich, {
    isolate(updateCounter_homerRef_enrich$i == 0)
    updateCounter_homerRef_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species_enrich, {
    isolate(updateCounter_homerRef_enrich$i == 0)
    updateCounter_homerRef_enrich <<- reactiveValues(i = 0)
  }) 
  output$homer_size_enrich <- renderUI({
    radioButtons('homer_size_enrich','Type of the region for motif finding',
                 c('given (exact size)'="given",
                   'custom size'="custom"
                 ),selected = "custom")
  })
  output$homer_size2_enrich <- renderUI({
    if(!is.null(input$homer_size_enrich)){
      if(input$homer_size_enrich == "custom"){
        numericInput('homer_size2_enrich','Size of the region for motif finding',value=200, step=100)
      }}
  })
  output$homer_bg_enrich <- renderUI({
    radioButtons('homer_bg_enrich','Background sequence',
                 c('random'="random",
                   'bed files'="peakcalling"
                 ),selected = "random")
  })
  output$homer_bg2_enrich <- renderUI({
    if(!is.null(input$homer_bg_enrich)){
      if(input$homer_bg_enrich == "peakcalling"){
        fileInput('homer_bg2_enrich',
                  'Select bed files',
                  accept = c("bed","narrowPeak"),
                  multiple = TRUE,
                  width = "80%")
      }}
  })
  updateCounter_enrich <- reactiveValues(i = 0)
  
  observe({
    input$motifButton_enrich
    isolate({
      updateCounter_enrich$i <- updateCounter_enrich$i + 1
    })
  })
  
  
  #Restart
  observeEvent(Enrich_peak_call_files(), {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Group_homer, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$Species_enrich, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_unknown_enrich, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_size_enrich, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_size2_enrich, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_bg_enrich, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  observeEvent(input$homer_bg2_enrich, {
    isolate(updateCounter_enrich$i == 0)
    updateCounter_enrich <<- reactiveValues(i = 0)
  }) 
  Enrich_peak_call_files_homer <- reactive({
    if(is.null(input$homer_bg2_enrich)){
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$homer_bg2_enrich[, 1])){
        file <- input$homer_bg2_enrich[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$homer_bg2_enrich[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      if(length(names(files2)) > 1){
        files2 <- soGGi:::runConsensusRegions(GRangesList(files2), "none")
      }
      return(files2)
    }
  })
  output$homer_unknown_enrich <- renderUI({
    selectInput("homer_unknown_enrich","Type of enrichment analysis",c("known motif","known and de novo motifs"), selected = "known motif")
  })
  
  output$homer_sample_order_enrich <- renderUI({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && 
       !is.null(input$homer_unknown_enrich) && input$Species_enrich != "not selected"){
      selectInput("homer_sample_order_enrich", "Order of groups", names(Enrich_peak_call_files()),
                  selected=names(Enrich_peak_call_files()),multiple = TRUE)
    }
  })
  enrich_motif_enrich <- reactive({
    if(updateCounter_enrich$i > 0 && input$motifButton_enrich > 0 && !is.null(Enrich_peak_call_files()) 
       && input$Species_enrich != "not selected" && !is.null(input$homer_unknown_enrich)){
      if(input$homer_size_enrich == "given") size <- "given"
      if(input$homer_size_enrich == "custom") size <- input$homer_size2_enrich
      if(input$homer_bg_enrich == "peakcalling" && is.null(input$homer_bg2_enrich)){
        return(NULL)
      }else return(findMotif(df= Enrich_peak_call_files(), Species = input$Species_enrich,size=size,back = input$homer_bg_enrich,motif_length=input$homer_length_enrich,
                             motif=input$homer_unknown_enrich, other_data = Enrich_peak_call_files_homer(),type="Other",section="enrich"))
    }
  })
  
  output$motif_enrich_plot <- renderPlot({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && 
       !is.null(input$homer_unknown_enrich) && input$Species_enrich != "not selected"){
      homer_Motifplot(df = enrich_motif_enrich(),showCategory = input$enrich_showCategory,section="enrich",
                      group_order = input$homer_sample_order_enrich)
    }
  })
  
  output$download_motif_enrich_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- homer_Motifplot(df = enrich_motif_enrich(),showCategory = input$enrich_showCategory,section="enrich")
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  output$download_motif_enrich_table = downloadHandler(
    filename = function() {
      paste0("known_motif_table",".txt")
    },
    content = function(file){write.table(motif_table_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  denovo_motif_table_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      return(denovo_motif(enrich_motif_enrich()))
    }
  })
  
  output$denovo_motif_enrich_result <- DT::renderDT({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && input$Species_enrich != "not selected"){
      denovo_motif_table_enrich()
    }
  })
  
  motif_table_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      return(known_motif(enrich_motif_enrich()))
    }
  })
  
  output$motif_enrich_result <- DT::renderDT({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && input$Species_enrich != "not selected"){
      motif_table_enrich()
    }
  })
  
  output$download_denovo_motif_enrich_table = downloadHandler(
    filename = function() {
      paste0("denovo_motif_table",".txt")
    },
    content = function(file){write.table(denovo_motif_table_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_homer_report_enrich = downloadHandler(
    filename = function() {
      paste0("HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- enrich_motif_enrich()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  
  ##More-------------
  #Bed tool----------
  output$bed_file1_warning <- renderText({
    if(is.null(bed_peak_call_files1())) print("Please upload bed files (A)")
  })
  output$bed_file2_warning <- renderText({
    if(input$data_file_type_bed != "type1" && 
       is.null(bed_peak_call_files2())) print("Please upload bed files (B)")
  })
  
  bed_peak_call_files1 <- reactive({
    if(is.null(input$bed_data_file)){
      if(input$goButton_bed > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
        files[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$bed_data_file[, 1])){
        file <- input$bed_data_file[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$bed_data_file[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  bed_peak_call_files2 <- reactive({
    if(is.null(input$bed_data_file2)){
      if(input$goButton_bed > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
        files[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$bed_data_file2[, 1])){
        file <- input$bed_data_file2[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$bed_data_file2[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  bed_input<-reactive({
    if(input$data_file_type_bed == "type1"){
      uploaded_files = names(bed_peak_call_files1())
    }else if(input$data_file_type_bed == "type2"){
      uploaded_files <- c()
      for(i in 1:length(names(bed_peak_call_files1()))){
        uploaded_files <- c(uploaded_files, paste0(names(bed_peak_call_files1()[i])," - ", 
                                                   names(bed_peak_call_files2()[i])))
      }
    }else {
      if(input$intersect_bed == "default") mode <- "_intersect_"
      if(input$intersect_bed == "wa") mode <- "_intersect-wa_"
      if(input$intersect_bed == "exclude") mode <- "_intersect-v_"
      uploaded_files <- c()
      for(i in 1:length(names(bed_peak_call_files1()))){
        uploaded_files <- c(uploaded_files, paste0(names(bed_peak_call_files1()[i]),mode, 
                                                   names(bed_peak_call_files2()[i])))
      }
    }
    return(as.data.frame(uploaded_files))
  })
  
  output$bed_input <- DT::renderDataTable({
    bed_input()
  })
  
  output$bed_merge_dist <- renderUI({
    numericInput("bed_merge_dist","Gap size (bp)",value = 1000,min = 1)
  })
  
  interval_merge <- reactive({
    if(input$data_file_type_bed == "type1"){
      files <- bed_peak_call_files1()
      df <- list()
      for(name in names(files)){
        file <- files[[name]]
        df[[name]] <- bedtorch::merge_bed(file, max_dist = input$bed_merge_dist)
      }
    }else if(input$data_file_type_bed == "type2"){
      files <- bed_peak_call_files1()
      files2 <- bed_peak_call_files2()
      df <- list()
      for(i in 1:length(names(files))){
        name <- paste0(names(bed_peak_call_files1()[i])," - ", 
                       names(bed_peak_call_files2()[i]))
        df[[name]] <- bedtorch::subtract_bed(files[[i]],files2[[i]])
      }
    }else{
      if(input$intersect_bed == "default") mode <- "_intersect_"
      if(input$intersect_bed == "wa") mode <- "_intersect-wa_"
      if(input$intersect_bed == "exclude") mode <- "_intersect-v_"
      files <- bed_peak_call_files1()
      files2 <- bed_peak_call_files2()
      df <- list()
      for(i in 1:length(names(files))){
        name <- paste0(names(bed_peak_call_files1()[i]),mode, 
                       names(bed_peak_call_files2()[i]))
        df[[name]] <- bedtorch::intersect_bed(files[[i]],files2[[i]],mode=input$intersect_bed)
      }
    }
    return(df)
  })
  
  output$download_bed <- downloadHandler(
    filename = function() {
      paste("bed_files", "zip", sep=".")
    },
    content = function(fname) {
      withProgress(message = "Preparing download, please wait",{
        fs <- c()
        files <- interval_merge()
        for(name in names(files)){
          file_name <- paste0(name, ".bed")
          fs <- c(fs, file_name)
          write.table(as.data.frame(files[[name]]), file_name,quote = F, row.names = F, col.names = F,sep = "\t")
        }
      })
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  observeEvent(input$goButton_bed,({
    updateSelectInput(session,inputId = "Species_bed","Species_bed",species_list, selected = "Homo sapiens (hg19)")
  }))
  
})
