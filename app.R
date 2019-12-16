# 3 abas
## 1 Instructions + upload file
## 2 Tabel
## 3 Plot

# load packages
library(shiny)
library(DT)
library(tidyverse)
library(ggpubr)

# set up environment
ui = fluidPage(
  
  titlePanel("RNA-seq analyzer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "regulation",
        label = "Regulation",
        choices = c("All", "Up", "Down", "Not"),
        selected = "All"
      ),
      selectInput(
        inputId = "significance",
        label = "Significance",
        choices = c("All", "Significant", "Not Significant"),
        selected = "All"
      ),
      downloadButton("downloadPlot", "Download Plot"),
      hr(),
      downloadButton("downloadData", "Download Table")
    ),
    
    
    mainPanel(
      tabsetPanel(
        
        tabPanel("Start", verbatimTextOutput("start"),
                 strong('Welcome,'),
                 
                 p('with this app, you will be able to analyze and organize RNA-seq raw data, create a graph and an interactive table, which will help you to visualize your reuslts quickly and more efficiently.'),
                 p( "After analyzing the raw RNA-seq data on Bowtie (http://bowtie-bio.sourceforge.net/index.shtml) and Htseq (https://htseq.readthedocs.io/en/release_0.11.1/), you should have an .csv table with your data. Please, before uploading your data on the space avaiable below, read the 'Instructions' tab."),
                 p("After uploading your data, the app will automaticly create an MA plot based on the gene expression fold change and the gene mean expression. You can vizualize the plot on the 'Plot' tab and download it."),
                 p("The app will also generate an interactive table. On the 'Regulation' option you can select genes that are up regulated (fold change >= 2), down regulated (<= -2), or neither up or down regulated. You can also use the 'Significance' option to select for genes that present fold change that is statistically significant (p<=0.05)."),
                 p("You can download the table by clicking on the 'Download Table' buton."),
                 p("To search for specific information, you can use the 'Search' butom at the 'Tabble' tab."),
                 p("The summary table is a easy way to vizualize how many genes are up or dawn regulated."),
                 
                 
                 fileInput(inputId = 'fileup', label = 'Upload Data')),
        
        tabPanel("Instructions", verbatimTextOutput("instructions"),
                 strong ('Please, follow the instructions below before uploading your data'),
                 
                 p(' Once you have received the raw RNA-seq data as a FASTQ-format file containing the reads sequenced from an Next Generation Sequencing platform, you will align these reads to an annotated reference genome, and quantify the expression of genes. The first step is to creat an index that will align the raw data with the reference genome of the bacteria of interest. To do so, you will use Bowtie which is a short read aligner. Please, visit web site  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml.'),
                 p ('Command to build your index: bowtie2-build <reference_in> <bt2_base>'),
                 p ('`<reference _in>` is a comma-separated list of FASTA files containing the reference sequences to be aligned to, and `<bt2_base>` is the basename of the index files to write.'),
                 p ('On the next step you will align the data from the FASTAQ files with the genome which is the index that you built.'),
                 p ('Use the following command to align: bowtie2 -x ~/Directory/index -1 .ForwardSequence.fastq -2 ./ReverseSequence.fastq -S OutputName.sam -p 2'),
                 p ('After, you will count the genes that are being expressed using htseq-count and generate a table with your data.'),
                 
                 p ('Use the following command to gene count: htseq-count -s no -t CDS -a 1 -m intersection-nonempty -i ID OutputName.sam HP_GCF_000210895.1_ASM21089v1_genomic.gff > OutputName.txt'),
                 p('IMPORTANT: do not rename the columns in the table that you generated.')
                 
        ),
        
        tabPanel("Plot", plotOutput("plot")), 
        
        tabPanel("Table", dataTableOutput("table")),
        
        tabPanel( "Summary", dataTableOutput("table2") )
        
      )
    ),
    position = 'right'
  )
)

# Define server logic ----
server <- function(input, output) {
  
  filtered_regulation = reactive({
    
    df = read.csv(input$fileup$datapath, stringsAsFactors = F)
    df = df %>% mutate(fold.change = exp(log2FoldChange))
    df = df  %>% filter(fold.change != "#VALUE!") %>% mutate(fold.change = as.numeric(fold.change))
    df = df %>% filter(!is.na(padj)) %>% mutate(obs_id = row_number())
    df = df %>% mutate(fold_change_scaled = ifelse(fold.change < 1, -1/fold.change, fold.change),
                       regulation = ifelse(fold_change_scaled < -2, "down", ifelse(fold_change_scaled > 2, "up", "not")),
                       fold_change_scaled = round(fold_change_scaled,3),
                       log2FoldChange = round(log2FoldChange,3),
                       baseMean = round(baseMean,3),
                       padj = round(padj,3))
    
    if(input$regulation == "All") {
      df %>% select(obs_id, gene, gene_description, fold_change_scaled, log2FoldChange, baseMean, padj)
    }
    
    else {
      if(input$regulation == "Up") {
        df %>% filter(regulation == 'up') %>% select(obs_id, gene, gene_description, fold_change_scaled, log2FoldChange, baseMean, padj) }
      else {
        if(input$regulation == 'Down') {
          df %>% filter(regulation =='down') %>% select(obs_id, gene, gene_description, fold_change_scaled, log2FoldChange, baseMean, padj)
        }
        else {
          df %>% filter(regulation == 'not') %>% select(obs_id, gene, gene_description, fold_change_scaled, log2FoldChange, baseMean, padj)
        }
        
      }
    }
  })
  
  filtered_significance = reactive({
    if(input$significance == 'All') {
      filtered_regulation()
    }
    else{
      if(input$significance == 'Significant') {
        filtered_regulation() %>% filter(padj < 0.05)
      }
      else{
        filtered_regulation() %>% filter(padj >= 0.05)
      }
    }
  })
  
  plotInput = function() {
    req(input$fileup)
    df <- read.csv(input$fileup$datapath)
    df = df %>% filter(!is.na(padj))
    df = df %>% select(gene, baseMean, log2FoldChange, padj)
    maplot = ggmaplot(df,
                      fdr = 0.05, fc = 2, size = 2,
                      palette = c("#B31B21", "#1465AC", "darkgray"),
                      genenames = as.vector(df$gene),
                      legend = "top", top = 20,
                      font.label = c("bold", 11),
                      font.legend = "bold",
                      font.main = "bold",
                      ggtheme = ggplot2::theme_minimal())
    
    return(maplot)
  }
  
  output$plot = renderPlot({plotInput()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, plot = plotInput(), device = device)
    })
  
  output$table = renderDataTable({
    DT::datatable(data = filtered_significance(), options = list(pageLength = 10),
                  rownames = FALSE, class = 'display', escape = FALSE)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "filtered_table.csv"
    },
    content = function(con) {
      write.csv(filtered_significance(), con, row.names = FALSE)
    }
  )
  
  output$table2 = renderDataTable({
    req(input$fileup)
    
    
    df = read.csv(input$fileup$datapath, stringsAsFactors = F)
    df = df  %>% filter(fold.change != "#VALUE!") %>% mutate(fold.change = as.numeric(fold.change))
    df = df %>% filter(!is.na(padj)) %>% mutate(obs_id = row_number())
    df = df %>% mutate(fold_change_scaled = ifelse(fold.change < 1, -1/fold.change, fold.change),
                       regulation = ifelse(fold_change_scaled < -2, "down", ifelse(fold_change_scaled > 2, "up", "not")),
                       significance = ifelse(padj < 0.05, "Significant", "Not significant"))
    
    table = df %>% group_by(regulation, significance) %>% tally()
    return(table)
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)
# runApp(shinyApp(ui = ui, server = server), launch.browser = T)