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

# Define UI for application
ui <- fluidPage(
  # Application title
  titlePanel("zUMIs-config: Generate yaml file"),
  navlistPanel(id = "mainNav",widths = c(2, 10),
               tabPanel("Mandatory Parameters",
                        #set basic things
                        fluidRow(
                          column(6,wellPanel(
                            textInput(inputId="runID", label="Name of the run/project", placeholder="eg: my_zUMIs_run")
                          )),
                          column(6,wellPanel(
                            textInput(inputId="outDir", label="Full path to the output directory", placeholder="eg: /path/to/output")
                          ))
                        ),
                        
                        #STAR options
                        
                        
                        # slider input for number of fastq reads
                        fluidRow(
                          h4("Input options:",style = "padding-left: 20px;"),
                          column(4,wellPanel(
                            sliderInput("nfiles", "Number of Input fastq files:", min = 1, max = 4,value = 2)
                          )),
                          column(8,wellPanel(
                            uiOutput("fqUI")
                          ))
                        ),
                        
                        
                        fluidRow(
                          
                          column(3,wellPanel(
                            uiOutput("fqBCui")
                          ),offset = 1),
                          column(3,wellPanel(
                            uiOutput("fqUMIui")
                          )),
                          column(3,wellPanel(
                            uiOutput("fqCDNAui")
                          ))
                        ),
                        
                        fluidRow(
                          h4("Mapping/Reference options:",style = "padding-left: 20px;"),
                          column(6,wellPanel(
                            textInput(inputId="STARpath", label="Full path to the STAR index directory", placeholder="eg: /path/to/output"),
                            textInput(inputId="GTFpath", label="Full path to the annotation GTF file", placeholder="eg: /path/to/output"),
                            textInput(inputId="STARparams", label="Optional: Additional mapping parameters for STAR", value=""),
                            numericInput(inputId = "NUMadditionalFA",label = "Optional: Number of additional reference sequences (eg. ERCC.fa)",value = 0,min = 0,step = 1)
                          )),
                          column(6,wellPanel(
                            p(strong("Additional references:")),
                            uiOutput("refUI")
                          ))
                        )
                        
                        
               ),
               tabPanel("Optional Parameters",
                        fluidRow(
                          column(6,wellPanel(
                            h4("General options:"),
                            numericInput("numThreads","Number of CPU Threads:",value=8,min=1,step=1),
                            numericInput("memLimit","Memory limit (Gb) 0 = no limit:",value=0,step=1),
                            checkboxInput("makeStats",label="Produce summary statistics",value = T),
                            checkboxInput("useSLURM",label="Use SLURM workload manager",value = F),
                            selectInput("whichStage",label = "Start zUMIs from following stage:",choices = c("Filtering", "Mapping", "Counting", "Summarising"),selected = "Filtering",multiple = F)
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
                            checkboxInput("countIntrons",label="Also count introns & exon+intron",value = T),
                            textInput("downsamp",label = "Number of reads to downsample to. This value can be a fixed number of reads (e.g. 10000) or a desired range (e.g. 10000-20000). 0 invokes adaptive downsampling.",value = "0"),
                            selectInput("strand",label = "Is the library stranded?",choices = c("unstranded" = 0, "positively stranded" = 1, "negatively stranded" = 2),selected = "unstranded",multiple = F),
                            numericInput("HamDist","Hamming distance collapsing of UMI sequences:",value=0,min=0,max=5,step=1),
                            checkboxInput("doVelocity",label="Generate velocyto-compatible counting of intron-exon spanning reads. ATTENTION! This option is currently not implemented!",value = F),
                            checkboxInput("countPrimary",label="Count the primary Hits of multimapping reads towards gene expression levels?",value = T),
                            checkboxInput("twoPass",label="Perform basic STAR twoPass mapping?",value = T)
                          )),
                          column(6,wellPanel(
                            h4("Barcode options:"),
                            radioButtons(inputId = "barcodeChoice",label = "Type of barcode selection:",choices = c("Automatic","Number of top Barcodes","Barcod whitelist")),
                            uiOutput("barcodeUI"),
                            numericInput("HamBC","Hamming distance collapsing of close cell barcode sequences. ATTENTION! This option is currently not implemented!",value=0,min=0,max=5,step=1),
                            numericInput("nReadsBC","Keep only the cell barcodes with atleast n number of reads",value=100,min=1,max=5,step=1)
                            
                          ))
                        )
               ),
               tabPanel("Generate YAML!",
                        p(strong("Your YAML is nearly ready...")),
                        uiOutput("saveUI"),
                        actionButton(inputId = "SaveYAML",label = "Save YAML!")
                         ),
                         tabPanel("Load YAML",
                           p(strong("You can load an existing YAML for zUMIs here")),
                           textInput(inputId="inYAML", label="Full path to YAML file", placeholder="eg: /path/to.yaml"),
                           actionButton(inputId = "LoadYAML",label = "Load YAML!"),
                           p(strong("Please click the Load button twice for all options to be properly set in shiny! Sorry for the inconvenience."))
               )
  ))

# Define server logic
server <- function(input, output, session) {
  
  output$barcodeUI <- renderUI({
    switch(input$barcodeChoice,
           "Automatic" = p(em("Intact barcodes will be detected automatically.")),
           "Number of top Barcodes" = numericInput(inputId = "BCnum",label = "Number of barcodes to consider:",value = 100, min = 10, step = 1),
           "Barcod whitelist" = textInput(inputId = "BCfile",label = "File to barcode whitelist to use:", value = "/fullpath/to/file.txt")
    )
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
        textInput(inputId = paste0("fqpath_",i),label = paste("Full path to fastq file",i),placeholder = "/path/to/read.fq.gz")
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
      textInput(inputId = paste0("BC_",i),label = paste("Read",i,"BC"),placeholder = "eg: 1-6")
    })
    
  })
  
  output$fqUMIui <- renderUI({
    lapply(1:input$nfiles, function(i) {
      textInput(inputId = paste0("UMI_",i),label = paste("Read",i,"UMI"),placeholder = "eg: 7-16")
    })
  })
  
  output$fqCDNAui <- renderUI({
    lapply(1:input$nfiles, function(i) {
      textInput(inputId = paste0("cDNA_",i),label = paste("Read",i,"cDNA"),placeholder = "eg: 1-50")
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
  
  makeYAML<-function(input){
    #collect fastq file information
    seqf<-list()
    for(i in 1:input$nfiles){
      bc_struc<-c(input[[paste0("cDNA_",i)]],input[[paste0("BC_",i)]],input[[paste0("UMI_",i)]])
      names(bc_struc)<-c("cDNA","BC","UMI")
      bc_struc<-bc_struc[which( bc_struc != "" )]
      
      seqf[[i]] <- list(
        "name" = input[[paste0("fqpath_",i)]],
        "base_definition" = paste0(names(bc_struc),"(",bc_struc,")")
      )
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
        "BarcodeBinning" = input$HamBC,
        "nReadsperCell" = input$nReadsBC
      ),
      "counting_opts" = list(
        "introns" = input$countIntrons,
        "downsampling" = input$downsamp,
        "strand" = as.integer(input$strand),
        "Ham_Dist" = input$HamDist,
        "velocyto" = input$doVelocity,
        "primaryHit" = input$countPrimary,
        "twoPass" = input$twoPass
      ),
      "make_stats" = input$makeStats,
      "use_SLURM" = input$useSLURM,
      "which_Stage" = input$whichStage
    )
    return(y)
  }
  
  observeEvent(input$LoadYAML, {
    
    if(file.exists(input$inYAML)){
      print(paste("Loading",input$inYAML))
      loading_func(input$inYAML)
      updateNavlistPanel(session = session, inputId = "mainNav",selected = "Mandatory Parameters")
      }else{
        print("File doesn't exist!")
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
        
      }
      
      updateNumericInput(session = session, inputId = "BCbases", value = ya$filter_cutoffs$BC_filter$num_bases)
      updateNumericInput(session = session, inputId = "BCphred", value = ya$filter_cutoffs$BC_filter$phred)
      updateNumericInput(session = session, inputId = "UMIbases", value = ya$filter_cutoffs$UMI_filter$num_bases)
      updateNumericInput(session = session, inputId = "UMIphred", value = ya$filter_cutoffs$UMI_filter$phred)
      updateNumericInput(session = session, inputId = "numThreads", value = ya$num_threads)
      updateNumericInput(session = session, inputId = "memLimit", value = ya$mem_limit)
      updateCheckboxInput(session = session, inputId = "makeStats", value = ya$make_stats)
      updateCheckboxInput(session = session, inputId = "useSLURM", value = ya$use_SLURM)
      updateSelectInput(session = session, inputId = "whichStage", selected = ya$which_Stage)
      updateCheckboxInput(session = session, inputId = "countIntrons", value = ya$counting_opts$introns)
      updateTextInput(session = session, inputId = "downsamp", value = ya$counting_opts$downsampling)
      updateSelectInput(session = session, inputId = "strand", selected = ya$counting_opts$strand)
      updateNumericInput(session = session, inputId = "HamDist", value = ya$counting_opts$Ham_Dist)
      updateCheckboxInput(session = session, inputId = "doVelocity", value = ya$counting_opts$velocyto)
      updateCheckboxInput(session = session, inputId = "countPrimary", value = ya$counting_opts$primaryHit)
      updateCheckboxInput(session = session, inputId = "twoPass", value = ya$counting_opts$twoPass)
      if (is.null(ya$barcodes$barcode_num) & is.null(ya$barcodes$barcode_file)) {
        updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Automatic")
      }
      if (!is.null(ya$barcodes$barcode_file)){
        updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Barcod whitelist")
        updateTextInput(session = session, inputId = "BCfile", value = ya$barcodes$barcode_file)
      }
      if (!is.null(ya$barcodes$barcode_num)){
        updateRadioButtons(session = session, inputId = "barcodeChoice", selected = "Number of top Barcodes")
        updateNumericInput(session = session, inputId = "BCnum", value = ya$barcodes$barcode_num)
      }
      
      updateNumericInput(session = session, inputId = "HamBC", value = ya$barcodes$BarcodeBinning)
      updateNumericInput(session = session, inputId = "nReadsBC", value = ya$barcodes$nReadsperCell)
      
      if(!is.null(ya$read_layout)){
        updateSelectInput(session = session, inputId = "layout", selected = ya$read_layout)
      }
      
  }
  
}

# Run the application
shinyApp(ui = ui, server = server)