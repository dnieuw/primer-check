library(shiny)
library(data.table)
library(stringr)
library(shinycssloaders)
library(shinyWidgets)
library(shinydashboard)
library(plotly)
library(lubridate)
library(reticulate)
library(rhandsontable)
library(Biostrings)

options(shiny.maxRequestSize=100*1024^2)
options(repos = BiocManager::repositories())

UPDATE_DATE<-"13th March 2020"

ui <- dashboardPage(skin = "red",
        dashboardHeader(title = "Primer-Check"),
        dashboardSidebar(
            h5(paste0("Updated on ", UPDATE_DATE)),
            tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon-96x96.png"),
                      tags$title("SARS-CoV-2 primer check")),
            a(href = "https://www6.erasmusmc.nl/viroscience/#", target = "_blank", img(src="virosciencelab_stamp-03.png", align = "middle", width="100%")),
            br(),
            sidebarMenuOutput("sidebarmenu"),
            br(),
            h5("The SARS-CoV-2 whole genomes were kindly shared via", a(href="https://www.gisaid.org/", target = "_blank", img(src="gisaid.png", align = "middle", width="50%"))),
            br()
        ), 
     dashboardBody(
         tabItems(
             tabItem(tabName = "input",
                     fluidPage(
                       tabBox(width = 8,
                              tabPanel(title = "Upload sequences",
                                   dropdownButton(h4("Upload a reference file containing whole genome sequences in fasta format"),
                                                  size = "xs",
                                                  circle = TRUE,
                                                  status = "info",
                                                  icon = icon("info"),
                                                  width = "100%",
                                                  tooltip = tooltipOptions(title = "Info")
                                  ),
                                  fileInput("ref_input","Reference file"),
                                  rHandsontableOutput("refs_table")
                               ),
                              tabPanel(title = "Upload oligos",
                                dropdownButton(h4("Upload a primer file for analysis in ';' separated format"),
                                               rHandsontableOutput("example_hot"),
                                               size = "xs",
                                               circle = TRUE,
                                               status = "info",
                                               icon = icon("info"),
                                               width = "100%",
                                               tooltip = tooltipOptions(title = "Info")
                                ),
                                fileInput("primer_input","Primer file"),
                                rHandsontableOutput("primer_table")
                              )
                       ),
                         column(width = 4,
                                box(title = "Add to existing database", width = NULL, status = "warning", solidHeader = TRUE,
                                    dropdownButton(h4("Select an existing database to load or to add primer/probe-reference alignments to or create a new database by typing in a database name"),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    uiOutput("availabledb"),
                                    actionButton("update_button", "UPDATE")
                                ),
                                box(title = "Create new database", width = NULL, status = "warning", solidHeader = TRUE,
                                    dropdownButton(h4("Select an existing database to load or to add primer/probe-reference alignments to or create a new database by typing in a database name"),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    textInput("newdb", "New database name:", value = ""),
                                    actionButton("newdb_button", "RUN")
                                ),
                                box(title = "Load results", width = NULL, status = "warning", solidHeader = TRUE,
                                    dropdownButton(h4("Select an existing database to load or to add primer/probe-reference alignments to or create a new database by typing in a database name"),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    textInput("newdb", "New database name:", value = ""),
                                    actionButton("newdb_button", "RUN")
                                ),
                                box(title = "Download results", width = NULL, status = "warning", solidHeader = TRUE,
                                    dropdownButton(h4("Click the LOAD/RUN button to either load old results or match the primer input with the reference input."),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    withSpinner(verbatimTextOutput("script_output")),
                                    downloadButton('downloadData', 'Download')
                                )
                         )
                     )
             ),
             tabItem(tabName = "results",
                     fluidPage(
                         column(width = 4,
                                box(width = NULL,
                                    dropdownButton(h4("Click a primer or probe to select it (red indicates one or more mismatches)"),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    uiOutput("origin"),
                                    # div(style = 'height: 800px; overflow: auto', withSpinner(plotlyOutput("query_select_heatmap")))
                                    withSpinner(plotlyOutput("query_select_heatmap", height = "800px"))
                                )
                         ),
                         column(width = 8,
                                box(width = NULL,
                                    dropdownButton(h4("Click a mismatch (number) to show it as an alignment below"),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    uiOutput("aln_heatmap_title"),
                                    plotlyOutput("aln_heatmap", height = "400px"),
                                ),
                                box(width = NULL,
                                    div(style = 'height: 400px; overflow: auto', htmlOutput("alignment"))
                                )
                         )
                     )
             ),
             tabItem(tabName = "source",
                 fluidPage(
                     box(width=12,
                         h4("We gratefully acknowledge the Authors, the Originating and Submitting Laboratories for their sequence and metadata shared through GISAID, on which this research is based"),
                         h4("All submitters of data may be contacted directly via ", a("www.gisaid.org",href="https://www.gisaid.org/")),
                         div(style = 'height: 800px; overflow:auto',tableOutput("source_table"))
                     )
                 )
             )
         )
     )
)

server <- function(input, output, session) {
  #Pretend that we are at shinyapps.io
  Sys.setenv(R_CONFIG_ACTIVE="shinyapps")
  
  #########################
  #Information window
  #########################
  
  showModal(modalDialog(
    title = "INFO",
    HTML(
      'This app visualizes the in silico performance of primers and probes designed for the detection of SARS-CoV-2<br><br>
      The centre of origin of the primer/probe sets can be selected in the dropdown menu<br><br>
      Primers and probes and their match/mismatch profile can be visualized by clicking on the interactive diagrams<br><br>
      For more information please click the (i) icons'
    ),
    easyClose = TRUE,
    footer = NULL
  ))
  
  #########################
  #Sidebarmenu server logic
  #########################  

  output$sidebarmenu <- renderMenu({
      #Remove input if run on shinyapps.io
      if (Sys.getenv("R_CONFIG_ACTIVE") == "shinyapps"){
          sidebarMenu(
              menuItem("Results", tabName = "results", icon = icon("th")),
              menuItem("Acknowledgements", tabName = "source", icon = icon("thumbs-up"))
          )
      } else {
          sidebarMenu(
              menuItem("Input", tabName = "input", icon = icon("dashboard")),
              menuItem("Results", tabName = "results", icon = icon("th")),
              menuItem("Acknowledgements", tabName = "source", icon = icon("thumbs-up"))
          )
      }
  })
  
  #Set the conda environment if not on shinyapps.io
  if (Sys.getenv("R_CONFIG_ACTIVE") != "shinyapps"){
      use_condaenv("base")
      source_python("primer_check_reticulate.py")
  }
  
  #########################
  #Input menu server logic
  #########################
  
  output$example_hot <- renderRHandsontable({
      data <- fread("name	sequence	type	Institute	gene	set
Target 1 	CCCTGTGGGTTTTACACTTAA	FWD	Institute X	Orf1b	1
Target 1	ACGATTGTGCATCAGCTGA	PROBE	Institute X	Orf1b	1
Target 1 	CCCTGTGGGTTTTACACTTAA	RE	Institute X	Orf1b	1
Target 2	ACGATTGTGCATCAGCTGA	FWD	Institute X	Orf1b	2
Target 2 	CCCTGTGGGTTTTACACTTAA	PROBE	Institute X	Orf1b	2
Target 2	ACGATTGTGCATCAGCTGA	RE	Institute X	Orf1b	2",header=T)
      rhandsontable(data, readOnly = T)
  })
  
  output$refs_table <- renderRHandsontable({
      req(input$ref_input)
      refs <- names(readDNAStringSet(input$ref_input$datapath))
      refs <- data.table(fasta_header=refs)
      rhandsontable(refs, readOnly = TRUE, height = 600)
  })
  
  output$primer_table <- renderRHandsontable({
      req(input$primer_input)
      primers <- fread(input$primer_input$datapath)
      rhandsontable(primers, readOnly = TRUE, height = 600)
  })

  result_file <- reactive({
      return(tempfile())
  })
  
  output$availabledb <- renderUI({
      pickerInput("availabledb",
                  label = "Available databases",
                  choices = list.files(".", pattern = "*.db"),
                  selected = "defaultDB.db",
                  options = list(
                      `live-search` = TRUE))
  })
  
  database_file <- reactive({
      if (isTruthy(input$newdb)){
          return(paste0(input$newdb,".db"))
      } else if (isTruthy(input$availabledb)) {
          return(input$availabledb)
      } else {
          return("defaultDB.db")
      }
  })
  
  primer_check_output <- eventReactive(input$run_button,{

      if (!isTruthy(input$ref_input) | !isTruthy(input$primer_input)){
          #If no primer or ref file is loaded load the old results from the database
          sys_out <- py_capture_output(run_primer_check(probes=NULL, 
                                       fasta=NULL, 
                                       out=result_file(), 
                                       database=database_file(), 
                                       cores=1, 
                                       loaddb=TRUE))
      } else {
          #Otherwise check primers versus references
          sys_out <- py_capture_output(run_primer_check(probes=input$primer_input$datapath,
                            fasta=input$ref_input$datapath, 
                            out=result_file(), 
                            database=database_file(), 
                            cores=32, 
                            loaddb=FALSE))
      }

      return(paste(sys_out, collapse = '\n'))
  })
  
  output$script_output <- renderText({
      req(primer_check_output())
      return(primer_check_output())
  })
  
  #Retrieve alignment data from script output or from disk depending on shinyapps.io or not
  if (Sys.getenv("R_CONFIG_ACTIVE") == "shinyapps"){
      aln_data <- reactive({
          data <- fread("current_results.csv", sep='\t')
          names(data) <- c("qname","pseq","type","origin","gene","set","refname","score","start","stop","refseq","alnstr","qseq")
          return(data)
      })
  } else {
      aln_data <- reactive({
          req(primer_check_output())
          data <- fread(result_file(), sep='\t')
          names(data) <- c("qname","pseq","type","origin","gene","set","refname","score","start","stop","refseq","alnstr","qseq")
          return(data)
      })
  }
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(aln_data(), con, row.names = F, quote = F, sep = "\t")
    }
  )
  
  #########################
  #Result menu server logic
  #########################
  
  check_mut <- function(aln_str, ref){
      #If the query is the same as the ref return a '.' otherwise return the difference in the ref
      return(ifelse(aln_str=='|','*',ref))
  }
  
  output$origin <- renderUI({
      data <- req(aln_data())
      origins <- data[,.N,origin][,origin]
      pickerInput("origin",
                  label = "Institute",
                  choices = origins)
  })
  
  output$query_select_heatmap <- renderPlotly({
      data <- req(aln_data())
      
      data[,gene:=factor(gene)]
      
      data <- data[origin==req(input$origin)]
      
      #Find all primers with one or more errors
      data[str_detect(alnstr,pattern = "\\.|X"), mm:=1]
      
      #Sum mismatches
      data <- data[,.("tot_mm"=sum(mm, na.rm = T)),.(qname, origin, gene, type, set)]
      
      data[,mm:=ifelse(tot_mm>0, "yes","no")]
      
      data[,yfactor:=factor(paste(qname,type,sep='|'), levels = rev(paste(qname,type,sep='|')), ordered=T)]
      
      # if(nrow(data)<25){
      #   p <- ggplot(data, aes(y=yfactor, x=gene, fill=mm, text=qname, key=qname)) +
      #     geom_tile() +
      #     scale_x_discrete(position = "top", drop=F) +
      #     scale_fill_manual(values=c("no"="lightgray","yes"="#e72807")) +
      #     facet_grid(gene + set ~ ., scales = "free_y") +
      #     theme_minimal() +
      #     theme(legend.position = "none",
      #           strip.text = element_blank(),
      #           text = element_text(size=14),
      #           axis.text.x = element_text(angle=45, vjust=1, hjust=0)) +
      #     labs(x=NULL, y=NULL)
      # } else {
      #   p <- ggplot(data, aes(y=yfactor, x=gene, fill=mm, text=qname, key=qname)) +
      #     geom_tile() +
      #     scale_x_discrete(position = "top", drop=F) +
      #     scale_fill_manual(values=c("no"="lightgray","yes"="#e72807")) +
      #     theme_minimal() +
      #     theme(legend.position = "none",
      #           strip.text = element_blank(),
      #           text = element_text(size=14),
      #           axis.text.x = element_text(angle=45, vjust=1, hjust=0)) +
      #     labs(x=NULL, y=NULL)
      # }
      
      p <- ggplot(data, aes(y=yfactor, x=gene, fill=mm, text=qname, key=qname)) +
        geom_tile() +
        scale_x_discrete(position = "top", drop=F) +
        scale_fill_manual(values=c("no"="lightgray","yes"="#e72807")) +
        facet_grid(gene + set ~ ., scales = "free_y") +
        theme_minimal() +
        theme(legend.position = "none",
              strip.text = element_blank(),
              text = element_text(size=14),
              axis.text.x = element_text(angle=45, vjust=1, hjust=0)) +
        labs(x=NULL, y=NULL)

      # # p <- ggplotly(p, tooltip = "text", height = 700 + nrow(data)*10, source = "plotly_query_select") %>% 

      p <- ggplotly(p, tooltip = "text", source = "plotly_query_select") %>%
        event_register('plotly_click') %>%
        config(displayModeBar = F) %>%
        layout(dragmode = F,
               xaxis=list(
                   tickfont=list(
                       size = 14,
                       color = "black")
               ),
               yaxis=list(
                   tickfont=list(
                       size = 14,
                       color = "black")
               )
        )
      
      return(p)
  })
  
  output$aln_heatmap_title <- renderUI({
      clickData <- event_data("plotly_click", source = "plotly_query_select")
      req(clickData)
      
      h3(clickData$key)
  })
  
  output$aln_heatmap <- renderPlotly({
      
      data <- req(aln_data())
      
      clickData <- event_data("plotly_click", source = "plotly_query_select")
      if (is.null(clickData)) return(NULL)
      
      data <- data[qname==clickData$key]
      
      #Split the aligned strings and check every position for a mutation
      mm_data <- data[,.(check_mut(str_split(alnstr,'', simplify = T), str_split(refseq,'', simplify = T)))]
      #Change colnames to query nucleotides
      colnames(mm_data) <- data[1, str_split(qseq,'', simplify = T)]
      #Add the refnames
      mm_data[,ref:=data[,refname]]
      #Melt the table to further process
      mm_melt <- melt(mm_data, id.vars = 'ref', variable.name = "query", value.name = "match")
      #Add query nucleotide positions
      mm_melt[,pos:=rep(seq(ncol(mm_data)-1),each=nrow(mm_data))]
      #Replace degenerate positions with "other"
      mm_melt[!match%in%c("A","C","T","G","*"),match:="other"]
      #Sum the matches by query position and match
      mm_melt <- mm_melt[,.(count=.N),.(query, pos, match)]
      #Include all possible matches
      mm_melt[,match:=factor(match, levels=c("*","A","C","T","G","other"), ordered = T)]
      
      nuc_labels <- as.character(mm_melt[,.N,.(query,pos)][,query])
      nuc_labels <- paste0('<b>',nuc_labels)
      match_labels <- paste0('<b>',c("*","A","C","T","G","other"))
      
      p <- ggplot(mm_melt) + 
          geom_text(aes(x=factor(pos),
                        y=match,
                        text = ifelse(match=='*','no mismatches',sprintf("%s mismatches<br>%s>%s", count, query, match)),
                        label=count), vjust=0.5) +
          scale_y_discrete(drop=F, label=match_labels)+
          scale_x_discrete(label=nuc_labels)+
          theme_minimal() +
          theme(legend.position = "none",
                text = element_text(size = 14))+
          labs(x=NULL, y=NULL)
      
      p <- ggplotly(p, tooltip = 'text', source = "plotly_muttable") %>% 
          event_register('plotly_click') %>% 
          config(displayModeBar = F) %>%
          layout(dragmode = F,
                 xaxis=list(
                     tickfont=list(
                         family = "Courier",
                         size = 20,
                         color = "black")
                 ),
                 yaxis=list(
                     tickfont=list(
                         family = "Courier",
                         size = 20,
                         color = "black")
                 )
          )
      
      return(p)
  })
  
  alignment_data <- reactive({
      req(aln_data())
      clickData <- event_data("plotly_click", source = "plotly_query_select")
      if (is.null(clickData)) print("NULL")
      
      data <- req(aln_data())
      data <- data[qname==clickData$key]
      
      return(data)
  })
  
  output$alignment <- renderUI({
      
      clickData <- event_data("plotly_click", source = "plotly_muttable")
      if (is.null(clickData)) return(NULL)
      #Degenerate positions or missing "*","A","C","T","G","R","Y","W","S","M","K","H","B","D","V","X","N","."
      
      mutations <- c("*","A","C","T","G","other")
      
      data <- isolate(alignment_data())
      
      if (mutations[clickData$y] == "other"){
          #Find which reference contains selected mutation at the selected position
          aln_data <- data[str_sub(refseq, clickData$x, clickData$x)%in%c("R","Y","W","S","M","K","H","B","D","V","X","N",".")]
      } else if (mutations[clickData$y] != "*") {
          aln_data <- data[str_sub(refseq, clickData$x, clickData$x)==mutations[clickData$y]]
      } else {
          aln_data <- copy(data)
      }
      
      add_color <- function(n){
          if (n=="A")
              return(paste0('<span style="background:#e72807; color:black">',n,'</span>'))
          if (n=="C")
              return(paste0('<span style="background:#009aff; color:black">',n,'</span>'))
          if (n=="G")
              return(paste0('<span style="background:#ffb507; color:black">',n,'</span>'))
          if (n=="T")
              return(paste0('<span style="background:#00bc0d; color:black">',n,'</span>'))
          else {
              return(paste0('<span style="background:#dfa3ff; color:black">',n,'</span>'))
          }
      }
      
      parse_alignment <- function(aln_str, ref){
          #If the query is the same as the ref return a '.' otherwise return the difference in the ref
          return(ifelse(aln_str=='|','.',sapply(ref, add_color)))
      }
      
      mm_data <- aln_data[,parse_alignment(str_split(alnstr,'', simplify = T), str_split(refseq,'', simplify = T))]
      
      aln_data[,alnstr:=apply(mm_data,1, paste, collapse='')]
      
      #Need to refactor this
      HTML(
          paste0(
              '<pre style="font-size:15pt">',
              paste(
                  #Paste primer
                  aln_data[1,paste(
                      paste('',
                      paste(sapply(str_split(qseq,'', simplify = T), add_color),
                          collapse=''),
                      '', sep='\t'),
                      qname, sep="\t")],
                  #Paste references
                  aln_data[,paste(
                      paste(
                          paste(`start`, alnstr, `stop`, sep='\t'), refname,
                          sep='\t'),
                      collapse = '<br>')], 
                  sep='<br>'))
          )
  })
  
  output$source_table <- renderTable({
      hot <- fread("gisaid_cov2020_acknowledgement_table.csv", header=T)
      return(hot)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

#Demo site
#rsconnect::deployApp('.', account = "dnieuw", appName = "primer-check", appTitle = "SARS-CoV-2 primer check visualization")
#Deploy site
#rsconnect::deployApp('.', account = "viroscience-emc", appName = "primer-check", appTitle = "SARS-CoV-2 primer check visualization")
