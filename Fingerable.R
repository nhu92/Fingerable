library(shiny)
library(shinyFiles)

# Source your original functions
source("chlgene_fingerprintable.R")

ui <- fluidPage(
  titlePanel("\U0001F9EC Gene Fingerprint Extractor"),
  sidebarLayout(
    sidebarPanel(
      textInput("gname1", "Gene Name 1 (required):", value = "", placeholder = "e.g., ycf1"),
      textInput("gname2", "Gene Name 2 (optional):", value = "", placeholder = "e.g., psbA"),
      shinyDirButton("dir", "Choose GenBank Files Directory", "Please select a directory"),
      textOutput("dirpath"),
      br(),
      actionButton("run", "\U0001F680 Run Analysis", class = "btn btn-success"),
      br(), br(),
      actionButton("combine", "\U0001F4CA Combine & Evaluate Classifications", class = "btn btn-primary")
    ),
    mainPanel(
      h4("\U0001F4AC Status Indicator"),
      verbatimTextOutput("status", placeholder = TRUE),
      br(),
      h4("\U0001F5C3 Final Classification Table"),
      tableOutput("final_table"),
      br(),
      h4("\U0001F4D1 Combined Summary Table"),
      tableOutput("combined_summary")
    )
  )
)

server <- function(input, output, session) {
  volumes <- c(Home = "~", getVolumes()())
  shinyDirChoose(input, "dir", roots = volumes)
  
  observe({
    if (!is.null(input$dir)) {
      dir_path <- parseDirPath(volumes, input$dir)
      output$dirpath <- renderText({ paste("Selected Directory:", dir_path) })
    }
  })
  
  observeEvent(input$run, {
    req(input$gname1)
    dir_path <- parseDirPath(volumes, input$dir)
    
    output$status <- renderText({ "Step 1: Extracting sequences..." })
    g1 <- input$gname1
    g2 <- input$gname2
    
    if (g2 == "") {
      out_fasta <- file.path(tempdir(), paste0(g1, ".fasta"))
      out_dist <- file.path(tempdir(), paste0(g1, "_distance.tsv"))
      out_class <- paste0(g1, "_classification.tsv")
      
      extract_shortest_interval_between_genes(g1, input_dir = dir_path, output_fasta = out_fasta)
      
      output$status <- renderText({ "Step 2: Aligning sequences and calculating distances..." })
      align_and_calc_distance(out_fasta, out_dist)
      
      output$status <- renderText({ "Step 3: Classifying species..." })
      result <- classify_species_by_distance(out_dist)
      write.table(result, out_class, sep = "\t", row.names = FALSE, quote = FALSE)
      
      output$status <- renderText({ paste("✅ Analysis completed for", g1) })
    } else {
      combo_label <- paste0(g1, "_", g2)
      out_fasta <- file.path(tempdir(), paste0(combo_label, ".fasta"))
      out_dist <- file.path(tempdir(), paste0(combo_label, "_distance.tsv"))
      out_class <- paste0(combo_label, "_classification.tsv")
      
      extract_shortest_interval_between_genes(g1, g2, input_dir = dir_path, output_fasta = out_fasta)
      
      output$status <- renderText({ "Step 2: Aligning sequences and calculating distances..." })
      align_and_calc_distance(out_fasta, out_dist)
      
      output$status <- renderText({ "Step 3: Classifying species..." })
      result <- classify_species_by_distance(out_dist)
      write.table(result, out_class, sep = "\t", row.names = FALSE, quote = FALSE)
      
      output$status <- renderText({ paste("✅ Analysis completed for", g1, "and", g2) })
    }
    
    output$final_table <- renderTable({
      if (input$gname2 == "") {
        read.delim(paste0(input$gname1, "_classification.tsv"))
      } else {
        read.delim(paste0(input$gname1, "_", input$gname2, "_classification.tsv"))
      }
    })
  })
  
  observeEvent(input$combine, {
    dir_path <- parseDirPath(volumes, input$dir)
    output$status <- renderText({ "\U0001F50D Analyzing gene combinations..." })
    summary_result <- combine_and_evaluate_classification(dir_path)
    write.table(summary_result, file = "gene_combination_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    output$combined_summary <- renderTable(summary_result)
    output$status <- renderText({ "\U0001F389 Gene combination analysis complete!" })
  })
}

shinyApp(ui, server)