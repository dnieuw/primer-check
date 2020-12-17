library(shiny)
library(data.table)
library(stringr)
library(shinycssloaders)
library(shinyWidgets)
library(shinydashboard)
library(plotly)
library(lubridate)
library(rhandsontable)
library(Biostrings)
library(ape)

options(shiny.maxRequestSize=100*1024^2)
options(repos = BiocManager::repositories())
options(rsconnect.max.bundle.size=1000*1024^2)

ui <- dashboardPage(skin = "red",
        dashboardHeader(title = "Primer-Check"),
        dashboardSidebar(
            tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon-96x96.png"),
                      tags$title("SARS-CoV-19 primer check")),
            a(href = "https://www6.erasmusmc.nl/viroscience/#", target = "_blank", img(src="virosciencelab_stamp-03.png", align = "middle", width="100%")),
            br(),
            sidebarMenuOutput("sidebarmenu"),
            br(),
            h5("The SARS-CoV-19 whole genomes were kindly shared via", a(href="https://www.gisaid.org/", target = "_blank", img(src="gisaid.png", align = "middle", width="50%"))),
            br()
        ),
     dashboardBody(
         tabItems(
             tabItem(tabName = "input",
                     fluidPage(
                       fluidRow(
                         column(width = 5,
                                materialSwitch(
                                  inputId = "example_data",
                                  label = "Use example data?", 
                                  value = FALSE,
                                  status = "warning"
                                ),
                                box(width = NULL,
                                    dropdownButton(h4("Upload an alignment file containing aligned sequences in fasta format. The first sequence is assumed to be the reference, which should perfectly fit all primers/probes.", strong("Large gaps in the primer/probe regions will mess up the primer alignment!")),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    fileInput("ref_input","Alignment file"),
                                    withSpinner(rHandsontableOutput("refs_table"))
                                )

                         ),
                         column(width = 6,
                                downloadButton("download_primer",label = "Download example primertable"),
                                box(width = NULL,
                                    dropdownButton(h4("Upload a primer file for analysis in comma separated format"),
                                                   rHandsontableOutput("example_hot"),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    fileInput("primer_input","Primer file"),
                                    withSpinner(rHandsontableOutput("primer_table"))
                                )
                          )
                         )
                     )
             ),
             tabItem(tabName = "results",
                     fluidPage(
                         column(width = 4,
                                box(width = NULL,
                                    dropdownButton(h4("Click a primer or probe to select it. Number indicates the amount of mismatched positions. A red tile indicates one or more mismatches."),
                                                   size = "xs",
                                                   circle = TRUE, 
                                                   status = "info",
                                                   icon = icon("info"),
                                                   width = "100%",
                                                   tooltip = tooltipOptions(title = "Info")
                                    ),
                                    uiOutput("origin"),
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
                                    withSpinner(plotlyOutput("aln_heatmap", height = "400px")),
                                ),
                                box(width = NULL,
                                    div(style = 'height: 400px; overflow: auto', withSpinner(htmlOutput("alignment")))
                                )
                         )
                     )
             ),
             tabItem(tabName = "source",
                 fluidPage(
                     box(width=12,
                         h4("We gratefully acknowledge the Authors, the Originating and Submitting Laboratories for their sequence and metadata shared through GISAID, on which this research is based"),
                         h4("All submitters of data may be contacted directly via ", a("www.gisaid.org",href="https://www.gisaid.org/")),
                         downloadButton("source_table","Download Acknowledgement table")
                     )
                 )
             )
         )
     )
)

server <- function(input, output, session) {
  
    f2si<-function (number, rounding=F, digits=ifelse(rounding, NA, 6)) {
      lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06, 
               0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21, 
               1e+24, 1e+27)
      pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "k", 
               "M", "G", "T", "P", "E", "Z", "Y", NA)
      ix <- findInterval(number, lut)
      if (ix>0 && ix<length(lut) && lut[ix]!=1) {
        if (rounding==T && !is.numeric(digits)) {
          sistring <- paste0(round(number/lut[ix]), pre[ix])
        } else if (rounding == T || is.numeric(digits)) {
          sistring <- paste0(signif(number/lut[ix], digits), pre[ix])
        } else {
          sistring <- paste0(number/lut[ix], pre[ix])
        } 
      } else {
        sistring <- as.character(number)
      }
      return(sistring)
    }
  
    output$sidebarmenu <- renderMenu({
      sidebarMenu(
          menuItem("Input", tabName = "input", icon = icon("dashboard")),
          menuItem("Results", tabName = "results", icon = icon("th")),
          menuItem("Acknowledgements", tabName = "source", icon = icon("thumbs-up"))
      )
    })
    
    output$download_primer <- downloadHandler(
      filename = function() {
        paste("example_data.csv")
      },
      content = function(file) {
        fwrite(fread("example_data.csv"), file, row.names = F)
      }
    )
    
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
    
    aln_data <- reactive({
      
      if (input$example_data & !isTruthy(input$ref_input)) {
        refs <- read.dna("example_data.fasta", format = "fasta")
      } else {
        req(input$ref_input)
        refs <- read.dna(input$ref_input$datapath, format = "fasta")
      }
      
      return(refs)
    })
    
    output$refs_table <- renderRHandsontable({
        refnames <- data.table(fasta_header=rownames(aln_data()))
        rhandsontable(refnames, readOnly = TRUE, height = 600)
    })
    
    output$primer_table <- renderRHandsontable({
      
      if (input$example_data & !isTruthy(input$primer_input)) {
        primers <- fread("example_data.csv")
      } else {
        req(input$primer_input)
        primers <- fread(input$primer_input$datapath)
      }
      
      rhandsontable(primers, readOnly = TRUE, height = 600)
    })

    output$origin <- renderUI({
      req(input$primer_table)
      primers <- hot_to_r(input$primer_table)
      origins <- primers[,.N,Origin][,Origin]
      pickerInput("origin",
                  label = "Institute",
                  choices = origins)
    })
    
    ambiguous_nt_map <- lapply(IUPAC_CODE_MAP, function(x) as.vector(str_split(x, '', simplify = T)))
    ambiguous_nt_map[["-"]] <- "-"
    IUPAC_NT <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N", "-")
    
    primer_positions <- reactive({
      
      #Assume that the first sequence in the file is the reference sequence
      if (input$example_data & !isTruthy(input$ref_input)) {
        ref <- readDNAStringSet("example_data.fasta", nrec = 1)
      } else {
        req(input$ref_input)
        ref <- readDNAStringSet(input$ref_input$datapath, nrec = 1)
      }
      
      req(input$primer_table)
      
      primer_data <- hot_to_r(input$primer_table)
      patterns <- DNAStringSet(primer_data$Sequence)
      #Reverse complement those primers that are in RE orientation
      patterns[which(primer_data[,Type=="RE"])] <- reverseComplement(patterns[which(primer_data[,Type=="RE"])])
      
      #Generate substitutionmatrix for pairwise alignment
      submat <- matrix(nrow = 16, ncol = 16)
      submat[1:15,1:15] <- nucleotideSubstitutionMatrix()
      submat[,16] <- 0
      submat[16,] <- 0
      submat[submat!=0] <- 1
      diag(submat) <- 1
      colnames(submat) <- IUPAC_NT
      rownames(submat) <- IUPAC_NT
      
      #Reduce score for matching N
      submat[15,1:14] <- 0.5
      submat[1:14,15] <- 0.5
      
      #Generate a pairwise alignment of all primers with the reference sequence
      result <- pairwiseAlignment(patterns, ref, substitutionMatrix = submat,  gapOpening = -1, gapExtension = -1, type = "local")
      result_score <- score(result)
      
      #Determine if there are any misalignments with the reference
      misalign_text <- unlist(sapply(seq(nrow(primer_data)), function(n){
        primer_name <- primer_data$Name[n]
        primer_seq <- primer_data$Sequence[n]
        if (str_length(primer_seq)!=result_score[n]) {
          if (primer_data$Type[n]=="RE") {
            return(paste0("<b>Primer: ", primer_name, '</b><br><pre style="font-size:12pt">', primer_seq, " - Primer<br>", reverseComplement(DNAStringSet(subject(result[n]))), " - Reference</pre><br>"))
          } else {
            return(paste0("<b>Primer: ", primer_name, '</b><br><pre style="font-size:12pt">', primer_seq, " - Primer<br>", subject(result[n]), " - Reference</pre><br>"))
          }
        } else {
          return(NULL)
        }
      }))
      
      #If there are any misalignments with the reference, inform the user
      if (!is.null(misalign_text)) {
        showModal(modalDialog(
          title = "Primer(s) did not align perfectly",
          div(style = 'height: 400px; overflow: auto',
            HTML(paste(c("<b>These primers did not align perfectly with your reference:</b><br>",names(ref),"<br><br>",misalign_text), collapse = ''))
          ),
          easyClose = TRUE,
          footer = "This may influence the results for these primers"
        ))
      }

      return(result)
    })
    
    calculate_basefrequencies <- reactive({
      req(aln_data())
      
      basefreq <- t(sapply(seq(ncol(aln_data())), function(n) base.freq(aln_data()[,n], freq=T, all = T)))
      
      return(basefreq)
    })
    
    compare_nt <- function(x,y){
      sapply(seq_along(x), function(n){
        if (x[n]=="N"|y[n]=="N")
          return(FALSE)
        any(ambiguous_nt_map[[x[n]]] %in% ambiguous_nt_map[[y[n]]])
      })
    }
    
    find_discrepancies <- reactive({
      primer_ranges <- ranges(Views(primer_positions()))
      primer_seq <- pattern(primer_positions())
      
      basefreq <- req(calculate_basefrequencies())
      
      result <- lapply(seq(length(primer_ranges)), function(n) {
        #Extract the location of the primer from the basefrequency table ignore N's
        primer_basefreq <- basefreq[start(primer_ranges[n]):end(primer_ranges[n]),-15]
        
        #Check each nucleotide in the primer for mismatches
        mismatch <- sapply(seq(width(primer_seq[n])), function(nt) {
          primer_nt <- str_sub(primer_seq[n],nt,nt)
          basefreq_nt <- str_to_upper(colnames(primer_basefreq))
          
          #Check which nucleotides mismatch
          mismatches <- !compare_nt(basefreq_nt, rep(primer_nt,ncol(primer_basefreq)))
          
          #Sum up each column of mismatches in the basefrequency table and return true if there is any mismatch
          return(sum(primer_basefreq[nt,mismatches])>0)
        })
        
        return(sum(mismatch))
      })
      
      return(result)
    })
    
    output$query_select_heatmap <- renderPlotly({

        discrep <- req(find_discrepancies())
        
        primer_data <- hot_to_r(input$primer_table)
        primer_data[,Gene:=factor(Gene)]
        primer_data[,tot_mm:=unlist(discrep)]
        
        primer_data <- primer_data[Origin==req(input$origin)]
        
        primer_data[,mm:=ifelse(tot_mm>0, "yes","no")]
        
        p <- ggplot(primer_data, aes(y=paste(Name,Type,sep='|'), x=Gene, fill=mm, key=Name)) + 
            geom_tile(aes(text=Name), width=0.95, height=0.95) + 
            geom_text(aes(label=tot_mm), vjust=0.5, color="white") +
            scale_x_discrete(position = "top", drop=F) +
            scale_fill_manual(values=c("no"="lightgray","yes"="#e72807")) + 
            facet_grid(Gene + Set ~ ., scales = "free_y") +
            theme_minimal() + 
            theme(legend.position = "none",
                  strip.text = element_blank(),
                  text = element_text(size=14), 
                  axis.text.x = element_text(angle=45, vjust=1, hjust=0)) +
            labs(x=NULL, y=NULL)
        
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
        
        #Remove hover text for the text elements:
        text_elements <- which(sapply(p$x$data, function(x) x$mode=="text"))
        p <- style(p, hoverinfo = "skip", traces = text_elements)
        
        return(p)
    })
    
    output$aln_heatmap_title <- renderUI({
        clickData <- event_data("plotly_click", source = "plotly_query_select")
        req(clickData)
        
        h3(clickData$key)
    })
    
    #Create reactive values to be able to hide figures when selecting different primers/primer origins
    select_data <- reactiveValues(alignment_click=NULL, muttable_click=NULL)
    
    observe({
      select_data$muttable_click <- event_data("plotly_click", source = "plotly_muttable")
    })
    
    observe({
      select_data$alignment_click <- event_data("plotly_click", source = "plotly_query_select")
    })
    
    observeEvent(event_data("plotly_click", source = "plotly_query_select"),{
      select_data$muttable_click <- NULL
    })
    
    observeEvent(input$origin,{
      select_data$alignment_click <- NULL
      select_data$muttable_click <- NULL
    })
    
    output$aln_heatmap <- renderPlotly({
        
      clickData <- select_data$alignment_click
      if (is.null(clickData)) return(NULL)
      
      primer_selected <- which(hot_to_r(input$primer_table)[,Name==unlist(clickData$key)])
      
      primer_data <- hot_to_r(input$primer_table)[primer_selected,]
      
      primer_ranges <- ranges(Views(primer_positions()))[primer_selected]
      primer_seq <- pattern(primer_positions())[primer_selected]
      
      basefreq <- req(calculate_basefrequencies())

      primer_basefreq <- data.table(basefreq[start(primer_ranges):end(primer_ranges),])
      
      if (primer_data$Type=="RE") {
        primer_seq <- reverseComplement(DNAStringSet(primer_seq))
        colnames(primer_basefreq)[1:14] <- as.character(reverseComplement(DNAStringSet(colnames(primer_basefreq)[1:14])))
        primer_basefreq <- primer_basefreq[rev(seq(nrow(primer_basefreq))),]
      }
      
      primer_basefreq[,pos:=seq(nrow(primer_basefreq))]
      primer_basefreq[,ref:=str_split(primer_seq, pattern="")[[1]]]
      
      primer_basefreq <- melt(primer_basefreq, id.vars =c("pos","ref"), variable.name = "char", value.name = "count")
      primer_basefreq[,char:=str_to_upper(char)]
      
      primer_basefreq[compare_nt(ref,char),char:="*"]
      
      primer_basefreq[!char%in%c("A","C","T","G","*","-","N"), char:="Other"]
      primer_basefreq[,char:=factor(char, levels=c("*","A","C","T","G","Other","-","N"), ordered = T)]
      primer_basefreq <- primer_basefreq[,.(count=sum(count)),.(ref,pos,char)]
      primer_basefreq <- primer_basefreq[count!=0]
      primer_basefreq <- primer_basefreq[,count:=sapply(count, f2si, digits = 2)]
      
      nuc_labels <- str_split(primer_seq, pattern = '', simplify = T)
      nuc_labels <- paste0('<b>',nuc_labels)
      match_labels <- paste0('<b>',c("A","C","T","G","*","-","Other","N"))
      names(match_labels) <- c("A","C","T","G","*","-","Other","N")
      
      p <- ggplot(primer_basefreq) +
        geom_text(aes(x=factor(pos),
                      y=char,
                      text = ifelse(char=='*','no mismatches',sprintf("%s mismatch(es)<br>%s>%s", count, ref, char)),
                      label=count), vjust=0.5) +
        scale_y_discrete(label=match_labels, drop=F)+
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

    output$alignment <- renderUI({
        
        mutations <- c("*","A","C","T","G","Other","-","N")

        clickData <- select_data$alignment_click
        if (is.null(clickData)) return(NULL)
        
        primer_selected <- which(hot_to_r(input$primer_table)[,Name==unlist(clickData$key)])
        
        primer_data <- hot_to_r(input$primer_table)[primer_selected,]
        
        primer_ranges <- ranges(Views(primer_positions()))[primer_selected]
        primer_seq <- pattern(primer_positions())[primer_selected]
        
        if (primer_data$Type=="RE") {
          primer_seq <- reverseComplement(DNAStringSet(primer_seq))
        }
        
        clickData <- select_data$muttable_click
        if (is.null(clickData)) return(NULL)
        
        aln_seqs <- aln_data()[,start(primer_ranges):end(primer_ranges)]
        
        if (primer_data$Type=="RE") {
          aln_seqs <- ape::complement(aln_seqs)
        }
        
        aln_nt <- str_to_upper(aln_seqs[,clickData$x])
        query_nt <- mutations[clickData$y]
        
        if (query_nt == "Other") {
          #Find which reference contains selected mutation at the selected position
          aln_selected <- which(aln_nt%in%c("R","Y","W","S","M","K","H","B","D","V","X","."))
        } else {
          aln_selected <- which(aln_nt==query_nt)
        }
        
        aln_seqs <- apply(as.character(aln_seqs[aln_selected,]), c(1,2), str_to_upper)
        
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
        
        aln_col <- matrix(nrow=nrow(aln_seqs),ncol=ncol(aln_seqs))
        
        for (n in seq(ncol(aln_seqs))) {
          query <- aln_seqs[,n]
          primer_nt <- rep_len(str_split(primer_seq, '', simplify = T)[,n], length(query))
          #If the query is the same as the ref return a '.' otherwise return the colored nucleotide
          aln_col[,n] <- ifelse(compare_nt(query,primer_nt),'.', sapply(query, add_color))
        }
        
        aln_data <- apply(aln_col, 1, paste0, collapse='')
        
        primer_text <- paste(sapply(str_split(primer_seq,'', simplify = T), add_color), collapse='')
        
        aln_text <- sapply(aln_data, function(x) paste(start(primer_ranges), x, end(primer_ranges), sep="\t"))
        aln_text <- paste(aln_text, rownames(aln_seqs), sep="\t")
        
        HTML(
            paste0(
                '<pre style="font-size:15pt">',
                  paste(
                    #Paste primer text
                    paste('\t', primer_text, '\t\t', primer_data$Name, sep=''),
                    #Paste references
                    paste(aln_text, collapse = '<br>'), 
                  sep='<br>')
            )
        )
    })
    
    output$source_table <- downloadHandler(
      filename = function() {
        paste("gisaid_cov2020_acknowledgement_table.pdf")
      },
      content = function(file) {
        file.copy("www/gisaid_hcov-19_acknowledgement_table.pdf", file)
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)