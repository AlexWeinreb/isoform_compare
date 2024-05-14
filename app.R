


# Initializations ----

#~ Packages ----
library(shiny)
library(dplyr)
library(ggplot2)


library(wbData)
gids <- wb_load_gene_ids(281)

#~ Data ----

t_exp_db <- pool::dbPool(RSQLite::SQLite(),
             dbname = "data/t_exp.sqlite.db",
             read_only = TRUE)

onStop(function() {
  pool::poolClose(t_exp_db)
})

tx_long <- tbl(t_exp_db, "t_exp")


measured_neurons <- qs::qread("data/measured_neurons.qs")


neurons_table <- readr::read_csv("data/neuron_properties.csv",
                                 col_types = "cccc") %>%
  filter(Neuron_type %in% measured_neurons) %>%
  mutate(across(Modality:Neurotransmitter, stringr::str_to_lower))
readr::stop_for_problems(neurons_table)


#~~ lookup tables ----
tsingle_color_scales <- tribble(
  ~shortcode,              ~fill_scale,                                    ~color_scale,
  'ggsci D3 (category20)', ggsci::scale_fill_d3("category20"),             ggsci::scale_color_d3("category20"),
  'viridis' ,              viridis::scale_fill_viridis(discrete = TRUE),   viridis::scale_color_viridis(discrete = TRUE),
  'ggsci npg',             ggsci::scale_fill_npg(),                        ggsci::scale_color_npg(),
  'iwanthue',              hues::scale_fill_iwanthue(),                    hues::scale_color_iwanthue())

theatmap_color_scales <- tribble(
  ~shortcode,           ~fill_scale,
  'MetBrewer OKeefe2',  scale_fill_gradientn(colors =
                                               MetBrewer::met.brewer("OKeeffe2",
                                                                     direction = -1)),
  'viridis magma',      viridis::scale_fill_viridis(direction = 1, option = "magma"),
  'Red-White-Blue',     scale_fill_gradient2())

scale_callback <- tribble(
  ~type,       ~callback,                              ~legend,
  "None",      function(x) x,                          "TPM",
  "Log2",      function(x) log2(x+1),                  "log2(TPM + 1)",
  "Z-score",   function(x) (x-mean(x))/sd(x),          "TPM Z-score",
  "Min-Max",   function(x) (x-min(x))/(max(x)-min(x)), "scaled TPM"
)


#~~ help text ----
explainNeurons_text <- 'Use ALL for all neurons, individual neuron names (e.g. "AWA", "ASEL", or "OLL"), or keywords such as "ACh", "motor", "sensory", ... The combination of all neurons corresponding to these keywords will be displayed.'

#~~ functions ----
source("R/utils.R")



# UI ----
ui <- fluidPage(
  
  titlePanel("Transcript-level quantification"),
  # Layout: 2 tabs, each with a sidebar Layout
  tabsetPanel(
    #~ Single gene ----
    tabPanel("Single gene",
             sidebarLayout(
               sidebarPanel(
                 # Inputs for tab Single gene (tsingle)
                 textInput(inputId = "tsingle_genes",
                           label = "Gene:",
                           value = "ric-4"),
                 fluidRow(
                   column(width = 10,
                          textInput(inputId = "tsingle_neurons",
                                    label = "Neurons:",
                                    value = "ALL")),
                   column(width = 1,
                          actionButton("tsingle_explainNeurons",
                                       label = "",
                                       icon = icon("question")))
                 ),
                 checkboxInput("plotIndividualSamples", "Plot individual samples"),
                 checkboxInput("tsingle_log_scale", "Display TPM on log scale"),
                 selectInput("tsingle_useColorScale",
                             label = "Color scale",
                             choices = tsingle_color_scales$shortcode),
                 actionButton("submit_tsingle",
                              "Plot"),
                 width = 2
               ),
               mainPanel(
                 fluidRow(
                   column(width = 2,
                          uiOutput("wormbase_browser_url"),
                          style = "margin-left:10px;margin-top:40px"),
                   column(width = 8,
                          plotly::plotlyOutput("single_gene_average",
                                               height = "150px",
                                               width = "90%"),
                          offset = 1,
                          align = "center"
                          ),
                 ),
                 plotly::plotlyOutput("single_gene_proportions",
                                      height = "600px"),
                 plotly::plotlyOutput("single_gene_tpm",
                                      height = "600px"),
                 width = 8
               )
             )
    ),
    tabPanel("Heatmap of transcript expression",
             #~ Heatmap ----
             sidebarLayout(
               sidebarPanel(
                 # Inputs for tab Heatmap (theatmap)
                 textInput(inputId = "theatmap_genes",
                           label = "Gene:",
                           value = "ric-4"),
                 fluidRow(
                   column(width = 10,
                          textInput(inputId = "theatmap_neurons",
                                    label = "Neurons:",
                                    value = "ALL")),
                   column(width = 1,
                          actionButton("theatmap_explainNeurons",
                                       label = "",
                                       icon = icon("question")))
                 ),
                 fluidRow(
                   column(width = 5,
                          selectInput("theatmap_scale_type",
                                      label = "Normalization",
                                      choices = scale_callback$type)),
                   column(width = 6,
                          selectInput("theatmap_scale_on",
                                      label = "Scale on",
                                      choices = c("Transcript", "Neuron")))
                 ),
                 selectInput("theatmap_useColorScale",
                             label =  "Color scale",
                             choices = theatmap_color_scales$shortcode),
                 actionButton("submit_theatmap",
                              "Plot"),
                 width = 2
               ),
               mainPanel(
                 plotly::plotlyOutput("heatmap"),
                 downloadButton("downloadData",
                                label = "Download table"),
                 downloadButton("downloadSvg",
                                label = "Download SVG"),
                 width = 9
               )
             )
    )
  )
)

# SERVER ----
server <- function(input, output, session) {
  
  
  # Reactives ----
  
  # For the tab Single gene
  r_tsingle_gene_id <- reactive({
    input$tsingle_genes %>%
      split_text_to_vector() %>%
      .[[1]] %>%
      convert_to_wb_id(gids)
  })
  
  r_tsingle_gene_name <- reactive({
    r_tsingle_gene_id() %>%
      i2s(gids)
  })
  
  r_tsingle_neurons <- reactive({
    input$tsingle_neurons %>%
      stringr::str_to_upper() %>%
      split_text_to_vector() %>%
      validate_neurons(neurons_table)
  })
  
  r_tsingle_colorChoice <- reactive(input$tsingle_useColorScale)
  
  r_tsingle_scaleFill <- reactive(tsingle_color_scales$fill_scale[
    tsingle_color_scales$shortcode == r_tsingle_colorChoice()
  ])
  r_tsingle_scaleColor <- reactive(tsingle_color_scales$color_scale[
    tsingle_color_scales$shortcode == r_tsingle_colorChoice()
  ])
  
  # data
  tsingle_data_from_db <- reactive({
    tx_long %>%
      filter(gene_id == !!r_tsingle_gene_id(),
             neuron_id %in% !!r_tsingle_neurons()) %>%
      collect()
  })
  
  
  
  
  # for the tab Heatmap
  r_theatmap_genes_id <- reactive({
    input$theatmap_genes %>%
      split_text_to_vector() %>%
      convert_to_wb_id(gids)
  })
  
  r_theatmap_neurons <- reactive({
    input$theatmap_neurons %>%
      stringr::str_to_upper() %>%
      split_text_to_vector() %>%
      validate_neurons(neurons_table)
  })
  
  r_theatmap_scaleFill <- reactive(theatmap_color_scales$fill_scale[
    theatmap_color_scales$shortcode==input$theatmap_useColorScale
  ])
  r_theatmap_scaleColor <- reactive(theatmap_color_scales$color_scale[
    theatmap_color_scales$shortcode==input$theatmap_useColorScale
  ])
  
  # Normalization
  r_theatmap_scale_callback <- reactive(scale_callback$callback[[
    which(scale_callback$type == input$theatmap_scale_type)
  ]])
  r_theatmap_scale_legend <- reactive(scale_callback$legend[[
    which(scale_callback$type == input$theatmap_scale_type)
  ]])
  
  heatmap_data <- reactive({
    tx_long %>%
      filter(gene_id %in% !!r_theatmap_genes_id(),
             neuron_id %in% !!r_theatmap_neurons()) %>%
      collect()
  })
  
  
  
  
  ## help bubbles ----
  observeEvent(input$tsingle_explainNeurons, {
    showModal(modalDialog(
      title = "Neurons selection",
      explainNeurons_text,
      easyClose = TRUE
    ))
  })
  observeEvent(input$theatmap_explainNeurons, {
    showModal(modalDialog(
      title = "Neurons selection",
      explainNeurons_text,
      easyClose = TRUE
    ))
  })
  
  
  
  # Single gene plots ----
  
  ##~ Wormbase URL ----
  r_wormbase_link <- eventReactive(
    eventExpr = input$submit_tsingle,
    valueExpr = tsingle_data_from_db()$browser_url[[1]])
  
  output$wormbase_browser_url <- renderUI(a("Open in Wormbase browser",
                                            href = r_wormbase_link(),
                                            class = "btn btn-default action-button", 
                                            style = "fontweight:800"))
  
  
  ##~ prepare average ----
  plot_tsingle_average <- eventReactive(
    eventExpr = input$submit_tsingle,
    valueExpr = {
      
      tsingle_data_from_db() %>%
        group_by(transcript_id, neuron_id) %>%
        summarize(mean_cnt = mean(TPM),
                  .groups = "drop_last") %>%
        summarize(mean_cnt = mean(mean_cnt),
                  .groups = "drop") %>%
        select(Transcript = transcript_id,
               `Transcript Usage` = mean_cnt) %>%
        ggplot(aes(x=1,
                   y=`Transcript Usage`,
                   fill = Transcript)) +
        ylab("Mean transcript usage (TPM)") +
        geom_col(position = position_stack()) +
        theme_classic() +
        coord_flip() +
        labs(x = NULL,
             title = paste0(r_tsingle_gene_name(), "/", r_tsingle_gene_id(),
                            " — average across neurons")) +
        r_tsingle_scaleFill() +
        scale_x_discrete(NULL)
      
      
    }
  )
  
  ##~ prepare proportions ----
  plot_tsingle_props <- eventReactive(
    eventExpr = input$submit_tsingle,
    valueExpr = {
      if(input$plotIndividualSamples){
        
        # by sample
        gg <- tsingle_data_from_db() %>%
          group_by(sample_id, gene_id) %>%
          mutate(sample_prop = round(100*TPM/sum(TPM), 1)) %>%
          rename(Neuron = neuron_id,
                 Transcript = transcript_id,
                 `Transcript Proportion` = sample_prop) %>%
          ggplot(aes(x=sample_id,
                     y=`Transcript Proportion`,
                     fill = Transcript)) +
          ylab("Transcript usage (%)")
        
      } else{
        
        # by neuron: first mean per neuron class, then prop of each tx
        gg <- tsingle_data_from_db() %>%
          group_by(transcript_id, neuron_id) %>%
          summarize(mean_cnt = mean(TPM),
                    sd_cnt = sd(TPM),
                    .groups = "drop") %>%
          group_by(neuron_id) %>%
          mutate(prop = round(100*mean_cnt/sum(mean_cnt), 1),
                 sd = sd_cnt/sum(mean_cnt)) %>%
          rename(Neuron = neuron_id,
                 Transcript = transcript_id,
                 `Transcript Usage` = prop) %>%
          ggplot(aes(x=Neuron,
                     y=`Transcript Usage`,
                     fill = Transcript)) +
          ylab("Mean transcript usage (%)")
      }
      
      gg <- gg +
        geom_col(position = position_stack()) +
        theme_classic() +
        coord_flip() +
        labs(x = NULL,
             title = paste0(r_tsingle_gene_name(), "/", r_tsingle_gene_id())) +
        scale_y_continuous(labels = pct) +
        r_tsingle_scaleFill()
      
      gg
    }
  )
  
  ##~ prepare tpm ----
  plot_tsingle_tpms <- eventReactive(
    eventExpr = input$submit_tsingle,
    valueExpr = {
      if(input$plotIndividualSamples){
        
        gg <- tsingle_data_from_db() %>%
          rename(Neuron = neuron_id,
                 Transcript = transcript_id,
                 `Transcript Usage` = TPM) %>%
          ggplot(aes(x = sample_id,
                     y = `Transcript Usage`,
                     color = Transcript)) +
          ylab("Transcript TPM")
        
      } else{
        
        gg <- tsingle_data_from_db() %>%
          group_by(transcript_id, neuron_id) %>%
          summarize(`Mean TPM` = mean(TPM),
                    sd = sd(TPM),
                    .groups = "drop") %>%
          rename(`Neuron` = neuron_id,
                 `Transcript` = transcript_id) %>%
          ggplot(aes(x=`Neuron`,
                     y=`Mean TPM`,
                     color = `Transcript`)) +
          geom_errorbar(aes(ymin = `Mean TPM` - sd,
                            ymax = `Mean TPM` + sd),
                        width = 0,
                        position = position_dodge(.2)) +
          ylab("Mean transcript TPM (+/-sd)")
        
      }
      
      gg <- gg +
        theme_classic() +
        geom_point(position = position_dodge(.2),
                   pch=15,
                   size=1) +
        coord_flip() +
        {if(input$tsingle_log_scale) scale_y_log10()} +
        labs(x = NULL,
             title = paste0(r_tsingle_gene_name(), "/", r_tsingle_gene_id())) +
        r_tsingle_scaleColor()
      
      gg
    }
  )
  
  #~ actual plots ----
  output$single_gene_proportions <- plotly::renderPlotly({
    plotly::ggplotly(plot_tsingle_props())
  })
  
  output$single_gene_tpm <- plotly::renderPlotly({
    plotly::ggplotly(plot_tsingle_tpms())
  })
  
  output$single_gene_average <- plotly::renderPlotly({
    plotly::ggplotly(plot_tsingle_average())
  })
  
  
  
  
  
  # Heatmap ----
  
  #~ prepare heatmap plot ----
  plot_theatmap <- eventReactive(
    eventExpr = input$submit_theatmap,
    valueExpr = {
      heatmap_data() %>%
        mutate(gene_id = i2s(gene_id, gids)) %>%
        group_by(gene_id, transcript_id, neuron_id) %>%
        summarize(TPM = mean(TPM),
                  .groups = "drop") %>%
        rename(Transcript = transcript_id,
               Neuron = neuron_id,
               Gene = gene_id) %>%
        arrange(Gene, Transcript, Neuron) %>%
        mutate(Transcript = forcats::fct_inorder(Transcript)) %>%
        group_by(.data[[input$theatmap_scale_on]]) %>%
        mutate(TPM = r_theatmap_scale_callback()(TPM)) %>%
        ungroup() %>%
        ggplot() +
        theme_minimal() +
        geom_tile(aes(x = Neuron, y =  Transcript, fill = TPM)) +
        geom_text(aes(x = 'Gene',
                      y = Transcript,
                      color = Gene,
                      label = Gene)) +
        scale_x_discrete(limits = c(unique(heatmap_data()$neuron_id), '','Gene',' ')) +
        r_theatmap_scaleFill() +
        theme(axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5)) +
        labs(title = paste0("Transcript expression level (",
                            r_theatmap_scale_legend(),
                            " — mean per neuron)"),
             fill = r_theatmap_scale_legend(),
             x = NULL,
             y = NULL)
    }
  )
  
  #~ render heatmap ----
  output$heatmap <- plotly::renderPlotly({
    plotly::ggplotly(plot_theatmap())
  })
  
  #~ download table ----
  output$downloadData <- downloadHandler(
    filename = function(){
      paste0(format(Sys.time(),"%y%m%d-%H%M"),
             "_",
             r_theatmap_genes_id() %>%
               i2s(gids) %>%
               paste0(collapse = "") %>%
               substr(start = 0,
                      stop = 20),
             "_",
             r_theatmap_neurons() %>%
               paste0(collapse = "") %>%
               substr(start = 0,
                      stop = 20)) %>%
        paste0(".tsv")
    },
    content = function(file){
      heatmap_data() %>%
        mutate(gene_id = i2s(gene_id, gids)) %>%
        group_by(gene_id, transcript_id, neuron_id) %>%
        summarize(TPM = mean(TPM),
                  .groups = "drop") %>%
        rename(Transcript = transcript_id,
               Neuron = neuron_id,
               Gene = gene_id) %>%
        arrange(Gene, Transcript, Neuron, TPM) %>%
        readr::write_tsv(file)
    }
  )
  
  #~ download SVG ----
  output$downloadSvg <- downloadHandler(
    filename = function(){
      paste0(format(Sys.time(),"%y%m%d-%H%M"),
             "_",
             r_theatmap_genes_id() %>%
               i2s(gids) %>%
               paste0(collapse = "") %>%
               substr(start = 0,
                      stop = 20),
             "_",
             r_theatmap_neurons() %>%
               paste0(collapse = "") %>%
               substr(start = 0,
                      stop = 20)) %>%
        paste0(".svg")
    },
    content = function(file){
      ggsave(file,
             plot = plot_theatmap(),
             width = 25,
             height = 10,
             units = "cm")
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)




















