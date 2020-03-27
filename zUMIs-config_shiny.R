#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(yaml)
library(shinyBS)

# Define UI for application
ui <- fluidPage(
  # Application title
  theme = shinythemes::shinytheme("simplex"),
  titlePanel("zUMIs-config: Generate yaml file"),
  navlistPanel(id = "mainNav",widths = c(2, 10),
               tabPanel("Mandatory Parameters",
                        #set basic things
                        fluidRow(
                          column(6,wellPanel(
                            textInput(inputId="runID", label="Name of the run/project:", placeholder="eg: my_zUMIs_run")
                          )),
                          shinyBS::bsTooltip(id="runID", title="This name will be used to label output files produced by zUMIs.", 
                                             placement = "bottom", trigger = "hover",options = list(container = "body")),
                          column(6,wellPanel(
                            textInput(inputId="outDir", label="Path to the output directory:", placeholder="eg: /path/to/output"),
                            shinyBS::bsTooltip(id="outDir", title="Please remember to give the full path, as relative paths may not work.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body"))
                          ))
                        ),



                        # slider input for number of fastq reads
                        fluidRow(
                          h4("Input options:",style = "padding-left: 20px;"),
                          column(4,wellPanel(
                            sliderInput("nfiles", "Number of input fastq files:", min = 1, max = 4,value = 2),
                            shinyBS::bsTooltip(id="nfiles", title="How many reads (including index reads) were obtained by your sequencing layout?", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            checkboxInput("patternsearch", "Search for sequence pattern in reads?", value = F),
                            checkboxInput("frameshiftsearch", "Correct for frameshift in barcode reads?", value = F),
                            uiOutput("patternUI"), uiOutput("patternReadUI"),
                            uiOutput("frameshiftUI"), uiOutput("frameshiftReadUI")
                          )),
                          column(8,wellPanel(
                            uiOutput("fqUI"),
                            shinyBS::bsTooltip(id="fqUI", title="Please remember to give full paths, as relative paths may not work.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body"))
                          ))
                        ),


                        fluidRow(
                          h4("In this section, fill only those fields that fit your input files:"),
                          column(3,wellPanel(
                            uiOutput("fqBCui"),
                            shinyBS::bsTooltip(id="fqBCui", title="List any barcode (BC) ranges to be extracted from the reads. You can also extract several barcode ranges from the same read using comma-separation: 1-6,11-16 ", 
                                               placement = "left", trigger = "hover",options = list(container = "body"))
                          ),offset = 1),
                          column(3,wellPanel(
                            uiOutput("fqUMIui"),
                            shinyBS::bsTooltip(id="fqUMIui", title="List any unique molecular identifier (UMI) ranges to be extracted from the reads.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body"))
                          )),
                          column(3,wellPanel(
                            uiOutput("fqCDNAui"),
                            shinyBS::bsTooltip(id="fqCDNAui", title="List the cDNA range(s) in the reads to be mapped to the genome. May be one range (single-end) or two ranges (paired-end)", 
                                               placement = "right", trigger = "hover",options = list(container = "body"))
                          ))
                        ),

                        fluidRow(
                          h4("Mapping/Reference options:",style = "padding-left: 20px;"),
                          column(6,wellPanel(
                            textInput(inputId="STARpath", label="Full path to the STAR index directory:", placeholder="eg: /path/to/output"),
                            shinyBS::bsTooltip(id="STARpath", title="The STAR index should be generated without splice junction database. Please remember to give the full path, as relative paths may not work.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            textInput(inputId="GTFpath", label="Full path to the annotation GTF file:", placeholder="eg: /path/to/output"),
                            shinyBS::bsTooltip(id="GTFpath", title="Make sure the gene annotation matches the genome. Please remember to give the full path, as relative paths may not work.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            textInput(inputId="STARparams", label="Optional: Additional mapping parameters for STAR:", value=""),
                            shinyBS::bsTooltip(id="STARparams", title="You may list additional STAR mapping parameters. For instance, try trimming adapter sequenes using --clip3pAdapterSeq", 
                                               placement = "top", trigger = "hover",options = list(container = "body")),
                            numericInput(inputId = "NUMadditionalFA",label = "Optional: Number of additional reference sequences:",value = 0,min = 0,step = 1),
                            shinyBS::bsTooltip(id="NUMadditionalFA", title="Here you can give additional reference sequences zUMIs should map to but are not integrated in the STAR index. For instance, you could add the ERCC spike-in reference on the fly (eg. ERCC.fa).", 
                                               placement = "top", trigger = "hover",options = list(container = "body"))
                          )),
                          column(6,wellPanel(
                            p(strong("Additional references:")),
                            uiOutput("refUI"),
                            shinyBS::bsTooltip(id="refUI", title="Please remember to give the full paths, as relative paths may not work.", 
                                               placement = "top", trigger = "hover",options = list(container = "body"))
                          ))
                        )


               ),
               tabPanel("Optional Parameters",
                        fluidRow(
                          column(6,wellPanel(
                            h4("General options:"),
                            numericInput("numThreads","Number of CPU Threads:",value=8,min=1,step=1),
                            shinyBS::bsTooltip(id="numThreads", title="More CPU Threads increase processing speeds, but speed usually becomes I/O limited above 32 Threads.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            numericInput("memLimit","Memory limit (Gb) 0 = no limit:",value=0,step=1),
                            shinyBS::bsTooltip(id="memLimit", title="More memory usually speeds up processing. Note that this setting cannot prevent STAR from using as much RAM as necessary for the specified reference genome.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            checkboxInput("makeStats",label="Produce summary statistics",value = T),
                            shinyBS::bsTooltip(id="makeStats", title="Plots and statistics files can be found in zUMIs_output/stats/.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            selectInput("whichStage",label = "Start zUMIs from following stage:",choices = c("Filtering", "Mapping", "Counting", "Summarising"),selected = "Filtering",multiple = F),
                            shinyBS::bsTooltip(id="whichStage", title="(Re)-Run the pipeline from the specified stage. Default: Filtering.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body"))
                          )),


                          column(6,wellPanel(
                            h4("Barcode & UMI filtering options:"),
                            numericInput("BCbases","Number of barcode bases below quality:",value=1,min=1,step=1),
                            numericInput("BCphred","Barcode Phred quality threshold:",value=20,min=1,max = 40,step=1),
                            numericInput("UMIbases","Number of UMI bases below quality:",value=1,min=1,step=1),
                            numericInput("UMIphred","UMI Phred quality threshold:",value=20,min=1,max = 40,step=1)
                          )
                          )),
                        fluidRow(
                          column(6,wellPanel(
                            h4("Counting options:"),
                            checkboxInput("countIntrons",label="Generate intron & exon+intron count tables?",value = T),
                            shinyBS::bsTooltip(id="countIntrons", title="Exonic read count tables will by output in any case. Untick this box if you want to count exonic reads only.", 
                                               placement = "top", trigger = "hover",options = list(container = "body")),
                            textInput("downsamp",label = "Number of reads to downsample to:",value = "0"),
                            shinyBS::bsTooltip(id="downsamp", title="This value can be a fixed number of reads (e.g. 10000) or a desired range (e.g. 10000-20000). You can also use several depths (e.g. 10000,20000,40000-50000). 0 invokes adaptive downsampling using median-absolute deviation of observed read counts.", 
                                               placement = "bottom", trigger = "hover",options = list(container = "body")),
                            selectInput("strand",label = "Is the library stranded?",choices = c("unstranded" = 0, "positively stranded" = 1, "negatively stranded" = 2),selected = "unstranded",multiple = F),
                            numericInput("HamDist","Hamming distance collapsing of UMI sequences:",value=0,min=0,max=5,step=1),
                            shinyBS::bsTooltip(id="HamDist", title="Note: Using this will considerably slow down the processing.", 
                                               placement = "top", trigger = "hover",options = list(container = "body")),
                            checkboxInput("writeHam", label="Write hamming distance mappings to file?", value = F),
                            checkboxInput("doVelocity",label="Generate RNA velocity counting of intron-exon spanning reads? Assumes velocyto is installed in path.",value = F),
                            checkboxInput("countPrimary",label="Count the primary Hits of multimapping reads towards gene expression levels?",value = T),
                            shinyBS::bsTooltip(id="countPrimary", title="Untick this box if you want to count uniquely aligned reads only.", 
                                               placement = "top", trigger = "hover",options = list(container = "body")),
                            checkboxInput("twoPass",label="Perform STAR twoPass mapping?",value = T),
                            shinyBS::bsTooltip(id="twoPass", title="Two-Pass mapping can improve the mapping by identifying novel splice-junctions.", 
                                               placement = "top", trigger = "hover",options = list(container = "body"))
                          )),
                          column(6,wellPanel(
                            h4("Barcode options:"),
                            radioButtons(inputId = "barcodeChoice",label = "Type of barcode selection:",choices = c("Automatic","Number of top Barcodes","Barcode whitelist")),
                            uiOutput("barcodeUI"),
                            numericInput("HamBC","Hamming distance collapsing of close cell barcode sequences.",value=1,min=0,max=5,step=1),
                            numericInput("nReadsBC","Keep only the cell barcodes with atleast n number of reads",value=100,min=1,max=5,step=1),
                            textInput("sharedBC",label = "Optional: Barcode Sharing (path to file):",value = NULL),
                            checkboxInput("demux", label = "Demultiplex into per-cell bam files?", value = F),
                            shinyBS::bsTooltip(id="demux", title = "Output files will be stored in zUMIs_output/demultiplexed/ .", 
                                               placement = "top", trigger = "hover",options = list(container = "body"))

                          )),
                          column(6,wellPanel(
                            h4("Dependencies:"),
                            textInput("r_exec","Rscript executable:", value = "Rscript"),
                            textInput("samtools_exec","samtools executable:", value = "samtools"),
                            textInput("pigz_exec","pigz executable:", value = "pigz"),
                            textInput("star_exec","STAR executable:", value = "STAR")
                          ))
                        )
               ),
               tabPanel("Generate YAML!",
                        p(strong("Your YAML is nearly ready...")),
                        uiOutput("saveUI"),
                        actionButton(inputId = "SaveYAML",label = "Save YAML!"),
                        p(),
                        p(strong("...or download your YAML file by clicking here:")),
                        downloadButton("downloadData", "Download YAML!")
               ),
               tabPanel("Load YAML",
                        p(strong("You can load an existing YAML for zUMIs here")),
                        textInput(inputId="inYAML", label="Full path to YAML file", placeholder="eg: /path/to.yaml"),
                        fileInput('file1', '...or choose YAML file to upload',
                                  accept = c(
                                    'text/tab-separated-values',
                                    'text/plain',
                                    '.yml',
                                    '.yaml'
                                  )),
                        actionButton(inputId = "LoadYAML",label = "Load YAML!"),
                        p(strong("Please click the Load button twice for all options to be properly set in shiny! Sorry for the inconvenience."))
               )
  ))

# Define server logic
server <- function(input, output, session) {

  output$barcodeUI <- renderUI({
    switch(input$barcodeChoice,
           #"Automatic" = p(em("Intact barcodes will be detected automatically.")),
           "Automatic" = textInput(inputId = "BCfile",label = "Optional: File to barcode whitelist to use for guiding automatic detection", value = NULL),
           "Number of top Barcodes" = numericInput(inputId = "BCnum",label = "Number of barcodes to consider:",value = 100, min = 10, step = 1),
           "Barcode whitelist" = textInput(inputId = "BCfile",label = "File to barcode whitelist to use:", value = "/fullpath/to/file.txt")
    )
  })
  
  output$patternUI <- renderUI({
    if(input$patternsearch==T){
           updateCheckboxInput(session = session, inputId = "frameshiftsearch",value = F)
           textInput(inputId = "pattern", label = "Search for the following sequence:" , placeholder = "ACTGCTGC")
    }
  })
  
  output$patternReadUI <- renderUI({
    if(input$patternsearch==T){
      selectInput(inputId = "patternRead", label = "Search for the pattern in this read:" ,choices = paste("Read",1:input$nfiles),selected = "Read 1")
    }
  })
  
  output$frameshiftUI <- renderUI({
    if(input$frameshiftsearch==T){
      updateCheckboxInput(session = session, inputId = "patternsearch",value = F)
      textInput(inputId = "pattern", label = "Correct frameshifts using the following sequence:" , placeholder = "ACTGCTGC")
    }
  })
  
  output$frameshiftReadUI <- renderUI({
    if(input$frameshiftsearch==T){
      selectInput(inputId = "patternRead", label = "Correct for frameshift with pattern in this read:" ,choices = paste("Read",1:input$nfiles),selected = "Read 1")
    }
  })

  output$refUI <- renderUI({
    if(input$NUMadditionalFA>0){
      lapply(1:input$NUMadditionalFA, function(i) {
        textInput(inputId = paste0("FA_",i),label = paste("Full path to additional sequence",i))
      })
    }
  })

  output$fqUI <- renderUI({
    if(input$nfiles>0){
      lapply(1:input$nfiles, function(i) {
        textInput(inputId = paste0("fqpath_",i),label = paste("Full path to fastq file",i,":"),placeholder = "/path/to/read.fq.gz")
      })
    }
  })

  output$fqBCui <- renderUI({
    # if(input$nfiles>0){
    #   lapply(1:input$nfiles, function(i) {
    #     checkboxGroupInput(inputId = paste0("fqptype_",i),choices = c("cDNA","BC","UMI"),label = paste("Content of read",i),inline = T)
    #   })
    # }
    lapply(1:input$nfiles, function(i) {
      textInput(inputId = paste0("BC_",i),label = paste("Read",i,"BC",":"),placeholder = "eg: 1-6")
    })

  })

  output$fqUMIui <- renderUI({
    lapply(1:input$nfiles, function(i) {
      textInput(inputId = paste0("UMI_",i),label = paste("Read",i,"UMI",":"),placeholder = "eg: 7-16")
    })
  })

  output$fqCDNAui <- renderUI({
    lapply(1:input$nfiles, function(i) {
      textInput(inputId = paste0("cDNA_",i),label = paste("Read",i,"cDNA",":"),placeholder = "eg: 1-50")
    })
  })


  output$saveUI <- renderUI({
    textInput(inputId = "savePath",label = "Save YAML file in this path:",value = paste0(input$outDir,"/",input$runID,".yaml"))
  })

  output$layoutUI <- renderUI({
    selectInput(inputId = "layout", label="cDNA Read Layout", choices = c("SE","PE"),selected = "SE")
  })

  # output$basedefui <- renderUI({
  #   strong(paste("found"))
  #   lapply(1:input$nfiles, function(i){
  #     #strong(paste0('Hi, this is output B#', i))
  #
  #     #lapply(input$paste0("fqptype_",i), function(j){
  #       #textInput(inputId = tmp,label = paste("read",i,"input",j))
  #     #})
  #   })
  # })


  observeEvent(input$SaveYAML, {
    write_yaml(
      x = makeYAML(input),
      file = input$savePath
    )

  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$runID, ".yaml", sep = "")
    },
    content = function(file) {
      write_yaml(
        makeYAML(input),
        file 
      )
      #write.csv(datasetInput(), file, row.names = FALSE)
    }
  )

  makeYAML<-function(input){
    #collect fastq file information
    seqf<-list()
    for(i in 1:input$nfiles){
      bc_struc<-c(input[[paste0("cDNA_",i)]],input[[paste0("BC_",i)]],input[[paste0("UMI_",i)]])
      names(bc_struc)<-c("cDNA","BC","UMI")
      bc_struc<-bc_struc[which( bc_struc != "" )]
      
      
      if(input$patternsearch==F & input$frameshiftsearch==F){
        seqf[[i]] <- list(
          "name" = input[[paste0("fqpath_",i)]],
          "base_definition" = paste0(names(bc_struc),"(",bc_struc,")")
        )
      }else{
        if(input$patternsearch==T & substr(input$patternRead,6,6)==i){
          seqf[[i]] <- list(
            "name" = input[[paste0("fqpath_",i)]],
            "base_definition" = paste0(names(bc_struc),"(",bc_struc,")"),
            "find_pattern" = input$pattern
          )  
        }else{
          if(input$frameshiftsearch==T & substr(input$patternRead,6,6)==i){
            seqf[[i]] <- list(
              "name" = input[[paste0("fqpath_",i)]],
              "base_definition" = paste0(names(bc_struc),"(",bc_struc,")"),
              "correct_frameshift" = input$pattern
            ) 
          }else{
          seqf[[i]] <- list(
            "name" = input[[paste0("fqpath_",i)]],
            "base_definition" = paste0(names(bc_struc),"(",bc_struc,")")
          )
          }
        }
      }
      
      
    }
    names(seqf) <- paste0("file",1:input$nfiles)

    #collect potential additional fasta files
    if(input$NUMadditionalFA>0){
      fa_list <- c()
      for(i in 1:input$NUMadditionalFA){
        fa_list <- c(fa_list,input[[paste0("FA_",i)]])
      }
    }else{
      fa_list<-NULL
    }

    #decide on barcode param
    #if(input$barcodeChoice)

    #make yaml list
    y <- list(
      "project" = input$runID,
      "sequence_files" = seqf,
      "reference" = list(
        "STAR_index" = input$STARpath,
        "GTF_file" = input$GTFpath,
        "additional_STAR_params" = input$STARparams,
        "additional_files" = fa_list
      ),
      "out_dir" = input$outDir,
      "num_threads" = input$numThreads,
      "mem_limit" = input$memLimit,
      "filter_cutoffs" = list(
        "BC_filter" = list(
          "num_bases" = input$BCbases,
          "phred" = input$BCphred
        ),
        "UMI_filter" = list(
          "num_bases" = input$UMIbases,
          "phred" = input$UMIphred
        )
      ),
      "barcodes" = list(
        "barcode_num" = input$BCnum,
        "barcode_file" = input$BCfile,
        "barcode_sharing" = input$sharedBC,
        "automatic" = ifelse(input$barcodeChoice=="Automatic", TRUE, FALSE),
        "BarcodeBinning" = input$HamBC,
        "nReadsperCell" = input$nReadsBC,
        "demultiplex" = input$demux
      ),
      "counting_opts" = list(
        "introns" = input$countIntrons,
        "downsampling" = input$downsamp,
        "strand" = as.integer(input$strand),
        "Ham_Dist" = input$HamDist,
        "write_ham" = input$writeHam,
        "velocyto" = input$doVelocity,
        "primaryHit" = input$countPrimary,
        "twoPass" = input$twoPass
      ),
      "make_stats" = input$makeStats,
      "which_Stage" = input$whichStage,
      "Rscript_exec" = input$r_exec,
      "STAR_exec" = input$star_exec,
      "pigz_exec" = input$pigz_exec,
      "samtools_exec" = input$samtools_exec
    )
    return(y)
  }

  observeEvent(input$LoadYAML, {
    
    if(file.exists(input$inYAML)){
      print(paste("Loading",input$inYAML))
      loading_func(input$inYAML)
      updateNavlistPanel(session = session, inputId = "mainNav",selected = "Mandatory Parameters")
    }else{
      if (is.null(input$file1)){
        print("File doesn't exist!")
      }else{
        print(paste("Loading",input$file1$datapath))
        loading_func(input$file1$datapath)
        updateNavlistPanel(session = session, inputId = "mainNav",selected = "Mandatory Parameters")
      }
      
      
    }
    
  })

  loading_func <- function(myYAML){

      ya <- read_yaml(file = myYAML)
      updateTextInput(session = session,inputId = "runID",value = ya$project)
      updateTextInput(session = session,inputId = "outDir",value = ya$out_dir)
      updateTextInput(session = session,inputId = "STARpath",value = ya$reference$STAR_index)
      updateTextInput(session = session,inputId = "GTFpath",value = ya$reference$GTF_file)
      updateTextInput(session = session,inputId = "STARparams",value = ya$reference$additional_STAR_params)
      updateNumericInput(session = session, inputId = "NUMadditionalFA", value = length(ya$reference$additional_files))
      if(length(ya$reference$additional_files)>0){
        for(i in 1:length(ya$reference$additional_files)){
          updateTextInput(session = session,inputId = paste0("FA_",i), value = ya$reference$additional_files[i])
        }
      }
      updateSliderInput(session = session, inputId = "nfiles", value = length(ya$sequence_files))
      for(i in 1:length(ya$sequence_files)){
        updateTextInput(session = session,inputId = paste0("fqpath_",i), value = ya$sequence_files[[i]]$name)

        if(any(grepl("BC",ya$sequence_files[[i]]$base_definition))){
          bc_string <- grep("BC",ya$sequence_files[[i]]$base_definition,value = T)
          bc_string <- substr(bc_string,start = 4, stop = nchar(bc_string)-1)
          updateTextInput(session = session,inputId = paste0("BC_",i), value = bc_string)
        }
        if(any(grepl("UMI",ya$sequence_files[[i]]$base_definition))){
          umi_string <- grep("UMI",ya$sequence_files[[i]]$base_definition,value = T)
          umi_string <- substr(umi_string,start = 5, stop = nchar(umi_string)-1)
          updateTextInput(session = session,inputId = paste0("UMI_",i), value = umi_string)
        }
        if(any(grepl("cDNA",ya$sequence_files[[i]]$base_definition))){
          cdna_string <- grep("cDNA",ya$sequence_files[[i]]$base_definition,value = T)
          cdna_string <- substr(cdna_string,start = 6, stop = nchar(cdna_string)-1)
          updateTextInput(session = session,inputId = paste0("cDNA_",i), value = cdna_string)
        }
        if(length(ya$sequence_files[[i]]$find_pattern)==1){
          updateCheckboxInput(session = session, inputId = "patternsearch", value = T)
          updateSelectInput(session = session, inputId = "patternRead", choices = paste("Read",1:length(ya$sequence_files)), selected = paste("Read",i))
          updateTextInput(session = session, inputId = "pattern", value = ya$sequence_files[[i]]$find_pattern)
        }
        if(length(ya$sequence_files[[i]]$correct_frameshift)==1){
          updateCheckboxInput(session = session, inputId = "frameshiftsearch", value = T)
          updateSelectInput(session = session, inputId = "patternRead", choices = paste("Read",1:length(ya$sequence_files)), selected = paste("Read",i))
          updateTextInput(session = session, inputId = "pattern", value = ya$sequence_files[[i]]$correct_frameshift)
        }
      }

      updateNumericInput(session = session, inputId = "BCbases", value = ya$filter_cutoffs$BC_filter$num_bases)
      updateNumericInput(session = session, inputId = "BCphred", value = ya$filter_cutoffs$BC_filter$phred)
      updateNumericInput(session = session, inputId = "UMIbases", value = ya$filter_cutoffs$UMI_filter$num_bases)
      updateNumericInput(session = session, inputId = "UMIphred", value = ya$filter_cutoffs$UMI_filter$phred)
      updateNumericInput(session = session, inputId = "numThreads", value = ya$num_threads)
      updateNumericInput(session = session, inputId = "memLimit", value = ya$mem_limit)
      updateCheckboxInput(session = session, inputId = "makeStats", value = ya$make_stats)
      updateSelectInput(session = session, inputId = "whichStage", selected = ya$which_Stage)
      updateCheckboxInput(session = session, inputId = "countIntrons", value = ya$counting_opts$introns)
      updateTextInput(session = session, inputId = "downsamp", value = ya$counting_opts$downsampling)
      updateSelectInput(session = session, inputId = "strand", selected = ya$counting_opts$strand)
      updateNumericInput(session = session, inputId = "HamDist", value = ya$counting_opts$Ham_Dist)
      updateCheckboxInput(session = session, inputId = "writeHam", value = ya$counting_opts$write_ham)
      updateCheckboxInput(session = session, inputId = "doVelocity", value = ya$counting_opts$velocyto)
      updateCheckboxInput(session = session, inputId = "countPrimary", value = ya$counting_opts$primaryHit)
      updateCheckboxInput(session = session, inputId = "twoPass", value = ya$counting_opts$twoPass)
      if (is.null(ya$barcodes$barcode_num) & ya$barcodes$automatic == TRUE) {
        updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Automatic")
        updateTextInput(session = session, inputId = "BCfile", value = ya$barcodes$barcode_file)
      }
      if (!is.null(ya$barcodes$barcode_file) & ya$barcodes$automatic == FALSE){
        updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Barcode whitelist")
        updateTextInput(session = session, inputId = "BCfile", value = ya$barcodes$barcode_file)
      }
      if (!is.null(ya$barcodes$barcode_num)){
        updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Number of top Barcodes")
        updateNumericInput(session = session, inputId = "BCnum", value = ya$barcodes$barcode_num)
      }

      updateNumericInput(session = session, inputId = "HamBC", value = ya$barcodes$BarcodeBinning)
      updateNumericInput(session = session, inputId = "nReadsBC", value = ya$barcodes$nReadsperCell)
      updateTextInput(session = session, inputId = "sharedBC", value = ya$barcodes$barcode_sharing)
      updateCheckboxInput(session = session, inputId = "demux", value = ya$barcodes$demultiplex)
      
      if(!is.null(ya$read_layout)){
        updateSelectInput(session = session, inputId = "layout", selected = ya$read_layout)
      }
      
      updateTextInput(session = session, inputId = "r_exec", value = ya$Rscript_exec)
      updateTextInput(session = session, inputId = "samtools_exec", value = ya$samtools_exec)
      updateTextInput(session = session, inputId = "pigz_exec", value = ya$pigz_exec)
      updateTextInput(session = session, inputId = "star_exec", value = ya$STAR_exec)
  }

}

# Run the application
shinyApp(ui = ui, server = server)
