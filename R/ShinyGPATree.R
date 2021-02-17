#' Run ShinyGPATree App
#'
#' This function will initialize the ShinyGPATree App for dynamic and interactive visualization of GPA-Tree model results.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class GPATree.
#' @return Output of GPA-Tree model.
#' @examples
#' \dontrun{
#' library(GPATree)
#'
#' # load GPATree example data
#' data(GPATreeExampleData)
#'
#' #fitting the GPATree model
#' fit <- GPATree(GPATreeExampleData$gwasPval, GPATreeExampleData$annMat)
#'
#' # initialize the ShinyGPATree app
#' ShinyGPATree(fit)
#' }
#' @export

ShinyGPATree <- function(object){

    ui <- shiny::fluidPage(
        shiny::tags$head(shiny::tags$style(".checkbox-inline {margin: 0 !important;}")),
        # App title ----
        shiny::titlePanel("ShinyGPATree"),
        shiny::tags$head(shiny::tags$style(shiny::HTML(".shiny-notification {position:fixed; width: 200px; top: 0%; right: 0%;}"))),

        # Sidebar layout with input and output definitions ----
        shiny::sidebarLayout(

            # Sidebar panel for inputs ----
            shiny::sidebarPanel(

                # Input: Slider for the value of log cp ----
                shiny::sliderInput(inputId = "log10cp",
                                   label = "Complexity parameter (cp), log10 scale",
                                   value = -5,
                                   min = -5,
                                   max = 0,
                                   step = 0.002,
                                   width= "95%"),
                shiny::helpText("Change cp to further prune the GPA-Tree model fit. To plot the original GPA-Tree model set pointer at -5."),
                shiny::hr(),
                shiny::h4("Plot options"),
                shiny::sliderInput(inputId = "plotHeight",
                                   label="Plot height",
                                   value = 350,
                                   min=100,
                                   max=2100,
                                   step = 50,
                                   width= "95%"),
                shiny::sliderInput(inputId = "plotWidth",
                                   label="Plot width",
                                   value = 600,
                                   min=400,
                                   max=4000,
                                   step = 50,
                                   width= "95%"),
                shiny::helpText("Change plot dimension using the options above.")

                ),

            # Main panel for displaying outputs ----
            shiny::mainPanel( #main panel ####

                shiny::tabsetPanel( # tabsetpanel ####
                    shiny::tabPanel("Plot",
                                    shiny::hr(),
                                    shiny::downloadButton("downloadPlot", label = "Download Plot"),
                                    shiny::htmlOutput('cpconversion'),
                                    shiny::plotOutput(outputId = "distPlot", height = "auto", width = "100%"),
                                    shiny::hr(),
                                    shiny::htmlOutput("anntext"),
                                    shiny::div(DT::dataTableOutput(outputId = "leafInfo", width = "100%"), style = "font-size:80%")
                                    ),
                    shiny::tabPanel("Info",
                                    shiny::wellPanel(
                                        shiny::fluidRow(
                                            shiny::column(5,
                                               shiny::numericInput(inputId = "FDR",
                                                                    label="False Discovery Rate (FDR) level",
                                                                    value = 0.01,
                                                                    min=0,
                                                                    max=1,
                                                                    step = 0.01),
                                               shiny::helpText("FDR must be between 0 and 1.")),
                                            shiny::column(4, offset = 2,
                                               shiny::radioButtons(inputId = 'fdrControl',
                                                                   label = 'Type of FDR',
                                                                   choices = list('global' = 'global',
                                                                                  'local' = 'local'),
                                                                   selected = 'global'))
                                        ),
                                        shiny::hr(),
                                        shiny::checkboxGroupInput(inputId="associd",
                                                                  label="Choose association status of SNPs",
                                                                  choices = list("Non-null" = 1, "Null" = 0),
                                                                  inline=TRUE,
                                                                  selected = 1,
                                                                  width = "80%"),
                                        shiny::checkboxGroupInput(inputId="leafindex",
                                                                  label="Choose Leaf",
                                                                  choices = list("LEAF 1" = "LEAF 1",
                                                                                 "LEAF 2" = "LEAF 2",
                                                                                 "LEAF 3" = "LEAF 3"),
                                                                  inline=TRUE,
                                                                  width = "80%",
                                                                  selected = c("LEAF 1")),
                                        shiny::downloadButton("downloadTable", label = "Download SNP Table")
                                        ),
                                        shiny::h5(shiny::tags$b("SNP Table")),
                                        shiny::div(shiny::tableOutput(outputId = 'tabletext'), style = "font-size:75%"),
                                        shiny::div(DT::dataTableOutput(outputId = "table", width = "100%"), style = "font-size:80%")
                                        )
                )
                )
            )
    )

    # column description for SNP Table ####
    nt <- function() {

        if (ncol(object@fit$Zmarg) == 1){
            col1 <- c("SNPID", "local FDR", 'p-value', 'Leaf', "others")
            col2 <- c("SNP ID.",
                      "1 - posterior probability of being non-null for the phenotype based on GPA-Tree model.",
                      "GWAS association p-value provided to the GPATree() function parameter 'gwasPval'.",
                      "Leaves in which SNPs fall.",
                      "Functional annotation information provided to the GPATree() function parameter 'annMat'.")
            df_out <- as.data.frame(cbind(col1, col2))
            colnames(df_out) <- c(" ","Column description")
        }
        return(df_out)
    }

    # defined function assocAll ####
    assocAll <- function(FDR, fdrControl, log10cp) {

        if (log10cp == -5){
            assoc_out <- assoc(object, FDR = FDR, fdrControl = fdrControl)
        } else {
            GPATreemodel_pruned <- prune(object, cp = 10^log10cp)
            assoc_out <- assoc(GPATreemodel_pruned, FDR = FDR, fdrControl = fdrControl)
        }
        assoc_out <- cbind(assoc_out, 1 - object@fit$Zmarg)

        if (ncol(object@fit$Zmarg) == 1){
            colnames(assoc_out)[3] <- 'local FDR'
        }

        if (is.null(rownames(assoc_out))) {
            assoc_out <- cbind('SNPID' = paste('SNP_', 1:nrow(assoc_out), sep = ''), assoc_out)
        } else {
            assoc_out <- cbind('SNPID' = rownames(assoc_out), assoc_out)
        }

        assoc_out <- cbind(assoc_out, object@gwasPval)
        ncol_start_ann <- ncol(assoc_out)

        if (ncol(object@fit$Zmarg) == 1){
            colnames(assoc_out)[5] <- 'p-value'
        }

        assoc_out <- cbind(assoc_out, object@annMat)
        colnames(assoc_out)[(ncol_start_ann+1):ncol(assoc_out)] <- colnames(object@annMat)
        return(assoc_out)
    }

    # defined function plotOutput ####
    plotOutput <- function(log10cp){

        if (log10cp == -5){
            plot(object)
        } else {
            GPATreemodel_pruned <- prune(object, cp = 10^log10cp)
            plot(GPATreemodel_pruned)
        }
    }

    #defined function leafinfo_output ####
    leafinfo_output <- function(object, log10cp){

        if (log10cp == -5){
            out <- leaf(object)
        } else {
            GPATreemodel_pruned <- prune(object, cp = 10^log10cp)
            out <- leaf(GPATreemodel_pruned)
        }

        out <- as.data.frame(cbind('Leaf' = rownames(out), out))

        # if (nrow(out) == 1){ out$Note <- 'No annotations selected' }
        return(out)
    }

    # server ####
    server <- function(input, output, session) {

        output$cpconversion <- shiny::renderText({

            if (input$log10cp == -5){
                paste(" ", "<b> GPA-Tree model<b>", sep="<br/>")
            } else {
                paste(" ", paste("<b>GPA-Tree model pruned using cp =", round(10^input$log10cp, 5),"<b>"), sep = '<br/>')
            }

            })

        pHeight <- shiny::reactive({

            new_height <- input$plotHeight
            return(new_height)

        })

        pWidth <-  shiny::reactive({

            new_width <- input$plotWidth
            return(new_width)

        })

        n0 <- shiny::reactive({

            leafinfo_out <- leafinfo_output(object,
                                            log10cp = input$log10cp)
            for (i in 1:ncol(object@fit$Zmarg)) {

                leafinfo_out[, i+1] <- formatC(as.numeric(leafinfo_out[, i+1]), format = "E", digits = 2)

            }

            return(leafinfo_out)

        })

        output$distPlot <- shiny::renderPlot(height = function()pHeight(), width = function()pWidth(), {

            # update of progress bar dependent on number of leaves in the tree
            n_index <- n0()
            n <- nrow(n_index)

            # Create a Progress object
            progress <- shiny::Progress$new(session, min = 1, max = 1.1*n)

            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress$close())
            progress$set(message = "Updating plot . . .", value = 0)

            for (i in 1:n) {

                # Update the progress bar
                progress$set(value = i)

                # Pause for 0.1 seconds to simulate a long computation.
                Sys.sleep(0.1)
            }

            cex_new <- (pWidth()*pHeight())/(pWidth()*pHeight())
            plotOutput(log10cp = input$log10cp)


        })

        output$anntext <- shiny::renderText({

            paste("<b><h5>Leaf Description<h5><b>")

        })

        output$leafInfo <- DT::renderDataTable({

            leafDat <- n0() %>%
                DT::datatable(rownames = FALSE,
                              class = 'cell-border stripe',
                              options = list(scrollY = "140px",
                                             searching = FALSE,
                                             lengthMenu = c(100, 1000),
                                             scipen = 4,
                                             scrollX = TRUE,
                                             paging = FALSE))
        })

        n1 <- shiny::reactive({

            asso_ans <- assocAll(FDR = input$FDR, fdrControl = input$fdrControl, log10cp = input$log10cp)
            return(asso_ans)

        })

        n2 <- shiny::reactive({

            assocOutnew <- n1()
            leaf_check <- c()
            for (i in 1:length(input$leafindex)) {
                leaf_out <- paste(stringr::str_split(input$leafindex[i], " ")[[1]][1],
                                  stringr::str_split(input$leafindex[i], " ")[[1]][2])
                leaf_check <- c(leaf_check, leaf_out)
            }

            select_rows <- c()
            for (i in 1:ncol(object@fit$Zmarg)) {
                select_rows <- c(select_rows, which(assocOutnew[, i+1] %in% input$associd))
            }

            assocOutnew <- assocOutnew[select_rows, ]
            assocOutnew <- subset(assocOutnew, assocOutnew$leaf %in% leaf_check)

            if (ncol(object@fit$Zmarg) == 1){

                assocOutnew <- as.data.frame(assocOutnew[, c('SNPID', 'local FDR', 'p-value',
                                                       'leaf', colnames(object@annMat))])
                assocOutnew <- assocOutnew[order(assocOutnew$`local FDR`), ]
            }

            for (i in 2:(2*ncol(object@fit$Zmarg)+1)) {
                assocOutnew[, i] <- formatC(as.numeric(assocOutnew[, i]), format = "E", digits = 2)
            }
            colnames(assocOutnew)[2*ncol(object@fit$Zmarg)+2] <- 'Leaf'

            return(assocOutnew)

        })

        output$tabletext <- shiny::renderTable({

            tout <- nt()
            return(tout)

        })

        output$table <-  DT::renderDT({

            # update of progress bar dependent on number of leaves in the tree
            n <- 20

            # Create a Progress object
            progress <- shiny::Progress$new(session, min = 1, max = 1.1*n)

            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress$close())
            progress$set(message = "Updating table . . .", value = 0)

            for (i in 1:n) {
                # Update the progress bar
                progress$set(value = i)
                # Pause for 0.1 seconds to simulate a long computation.
                Sys.sleep(0.1)
            }

            dtab <- n2() %>%
                        DT::datatable(rownames = FALSE,
                                      class = 'cell-border stripe',
                                      options = list(scrollY = "200px",
                                                     searching = TRUE,
                                                     # dom = 'b',
                                                     lengthMenu = c(100, 1000),
                                                     scipen = 4,
                                                     scrollX = TRUE,
                                                     paging = FALSE))
            return(dtab)

        })

        inputVar <- shiny::reactive({

            n1_out <- n0()
            if (ncol(object@fit$Zmarg) == 1){
                n1out_leaf <- paste(n1_out$Leaf,
                                    ' (local FDR = ',
                                    formatC(round(as.numeric(n1_out[, 2]), 3), 3, format="f"),
                                    ")", sep = '')
            }

            return(n1out_leaf)

        })

        shiny::observeEvent(inputVar(),  {shiny::updateCheckboxGroupInput(session,
                                                                          inputId="leafindex",
                                                                          label = "Choose Leaf",
                                                                          choices = inputVar(),
                                                                          selected = inputVar(),
                                                                          inline = TRUE)} )

        output$downloadTable <- shiny::downloadHandler(
            filename = function() {
                paste('ShinyGPATree-SNPtable-', Sys.time(), '.csv', sep='')
            },
            content = function(file) {
                write.csv(n2(), file, row.names = FALSE)
            })

        output$downloadPlot <- shiny::downloadHandler( filename = function() {
                paste('ShinyGPATreePlot-cp', round(10^input$log10cp, 3), '.png', sep='')
                },
            content = function(file) {
                grDevices::png(file,
                               height = pHeight(),
                               width = pWidth(),
                               # res = 300,
                               # bg = "transparent",
                               units = 'px')
                print(plotOutput(log10cp = input$log10cp))
                # print(p1())
                grDevices::dev.off()
            })
    }

    shiny::shinyApp(ui = ui, server = server)

}
