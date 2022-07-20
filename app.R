


# Initializations ----

#~ Packages ----
library(shiny)
library(dplyr)
library(ggplot2)


library(wbData)
gids <- wb_load_gene_ids(281)

#~ Data ----

t_exp_db <- DBI::dbConnect(RSQLite::SQLite(),
                           "data/t_exp.sqlite.db",
                           read_only = TRUE)

tx_long <- tbl(t_exp_db, "t_exp")


measured_neurons <- tx_long %>%
  select(neuron_id) %>%
  distinct() %>%
  pull(neuron_id)


neurons_table <- readr::read_csv("data/neuron_properties.csv",
                                 col_types = "cccc") %>%
  filter(Neuron_type %in% measured_neurons) %>%
  mutate(across(Modality:Neurotransmitter, stringr::str_to_lower))
readr::stop_for_problems(neurons_table)


accepted_color_scales <- tribble(
  ~shortcode,      ~fill_scale,                                     ~color_scale,
  'viridis' ,       viridis::scale_fill_viridis(discrete = TRUE),   viridis::scale_color_viridis(discrete = TRUE),
  'ggsci_d3_cat20', ggsci::scale_fill_d3("category20"),             ggsci::scale_color_d3("category20"),
  'ggsci_npg',      ggsci::scale_fill_npg(),                        ggsci::scale_color_npg(),
  'iwanthue',       hues::scale_fill_iwanthue(),                    hues::scale_color_iwanthue())


explainNeurons_text <- 'Use ALL for all neurons, individual neuron names (e.g. "AWA", "ASEL", or "OLL"), or keywords such as "ACh", "motor", "sensory", ... The combination of all neurons corresponding to these keywords will be displayed.'

#~ Functions ----
source("R/utils.R")



# UI ----
ui <- fluidPage(
  
  titlePanel("Transcript-level quantification"),
  # Layout: 2 tabs, each with a sidebar Layout
  tabsetPanel(
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
                 selectInput("tsingle_useColorScale", "Color scale", accepted_color_scales$shortcode),
                 width = 2
               ),
               mainPanel(
                 
                 plotly::plotlyOutput("single_gene_proportions"),
                 plotly::plotlyOutput("single_gene_tpm"),
                 width = 8
               )
             )
    ),
    tabPanel("Heatmap of transcript expression",
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
                 checkboxInput("theatmap_log_scale", "Display TPM on a log scale"),
                 selectInput("theatmap_useColorScale", "Color scale", accepted_color_scales$shortcode),
                 width = 2
               ),
               mainPanel(
                 plotly::plotlyOutput("heatmap"),
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
  
  r_tsingle_neurons <- input$tsingle_neurons %>%
    stringr::str_to_upper() %>%
    split_text_to_vector() %>%
    validate_neurons(neurons_table) %>%
    reactive()
  
  r_tsingle_colorChoice <- reactive(input$tsingle_useColorScale)
  
  r_tsingle_scaleFill <- reactive(accepted_color_scales$fill_scale[
    accepted_color_scales$shortcode == r_tsingle_colorChoice()
  ])
  r_tsingle_scaleColor <- reactive(accepted_color_scales$color_scale[
    accepted_color_scales$shortcode == r_tsingle_colorChoice()
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
  
  r_theatmap_neurons <- input$theatmap_neurons %>%
    stringr::str_to_upper() %>%
    split_text_to_vector() %>%
    validate_neurons(neurons_table) %>%
    reactive()
  
  r_theatmap_scaleFill <- reactive(accepted_color_scales$fill_scale[
    accepted_color_scales$shortcode==input$theatmap_useColorScale
  ])
  r_theatmap_scaleColor <- reactive(accepted_color_scales$color_scale[
    accepted_color_scales$shortcode==input$theatmap_useColorScale
  ])
  
  
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
  ##~ proportions ----
  output$single_gene_proportions <- plotly::renderPlotly({
    
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
    
    plotly::ggplotly(gg)
  })
  
  ##~ tpm ----
  output$single_gene_tpm <- plotly::renderPlotly({
    
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
    
    
    plotly::ggplotly(gg)
  })
  
  
  
  
  
  
  
  # Heatmap ----
  output$heatmap <- plotly::renderPlotly({
    
    
    plot_data <- tx_long %>%
      filter(gene_id %in% !!r_theatmap_genes_id(),
             neuron_id %in% !!r_theatmap_neurons()) %>%
      collect()
    
    gg <- plot_data %>%
      mutate(gene_id = i2s(gene_id, gids)) %>%
      group_by(gene_id, transcript_id, neuron_id) %>%
      summarize(`log(TPM)` = log1p(mean(TPM)),
                TPM = mean(TPM),
                .groups = "drop") %>%
      rename(Transcript = transcript_id,
             Neuron = neuron_id,
             Gene = gene_id) %>%
      arrange(Gene, Transcript, Neuron) %>%
      mutate(Transcript = forcats::fct_inorder(Transcript)) %>%
      ggplot() +
      theme_minimal() +
      {if(input$theatmap_log_scale) {
        geom_tile(aes(x = Neuron, y =  Transcript, fill = `log(TPM)`))} else {
          geom_tile(aes(x = Neuron, y =  Transcript, fill = TPM))}} +
      geom_text(aes(x = 'Gene',
                    y = Transcript,
                    color = Gene,
                    label = Gene)) +
      scale_x_discrete(limits = c(unique(plot_data$neuron_id), '','Gene',' ')) +
      # viridis::scale_fill_viridis(direction = 1, option = "magma") +
      scale_fill_gradientn(colors = MetBrewer::met.brewer("OKeeffe2",
                                                          direction = -1)) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5)) +
      {if(input$theatmap_log_scale) {
        labs(title = "Transcript expression level (log TPM — mean per neuron)", x=NULL, y=NULL)
      } else {
        labs(title = "Transcript expression level (mean TPM per neuron)", x=NULL, y=NULL)}}
    
    
    plotly::ggplotly(gg)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)




















