


# Initializations ----

#~ Packages ----
library(shiny)
library(dplyr)
library(ggplot2)


library(wbData)
gids <- wb_load_gene_ids(281)

#~ Data ----

t_exp_db <- DBI::dbConnect(duckdb::duckdb(),
                           "data/t_exp.duckdb.db",
                           read_only = TRUE)

tx_long <- tbl(t_exp_db, "t_exp")


measured_neurons <- tx_long |>
  select(neuron_id) |>
  distinct() |>
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


#~ Functions ----
source("R/utils.R")



# Define UI ----
ui <- fluidPage(
  
  # Application title
  titlePanel("Transcript-level quantification"),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "gene_id_or_name",
                label = "Gene:",
                value = "ric-4"),
      fluidRow(
        column(width = 10,
               textInput(inputId = "neurons",
                         label = "Neurons:",
                         value = "ALL")),
        column(width = 1,
               actionButton("explainNeurons",
                            label = "",
                            icon = icon("question")))
        ),
   checkboxInput("plotIndividualSamples", "Plot individual samples"),
   checkboxInput("log_scale", "Use log scale for TPM"),
   selectInput("useColorScale", "Color scale", accepted_color_scales$shortcode)
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Single gene",
               conditionalPanel('input.plotIndividualSamples == 0',
                                plotly::plotlyOutput("tx_neuron_proportions")),
               conditionalPanel('input.plotIndividualSamples == 0',
                                plotly::plotlyOutput("tx_neuron_tpm")),
               conditionalPanel('input.plotIndividualSamples == 1',
                                plotly::plotlyOutput("tx_samples_proportions")),
               conditionalPanel('input.plotIndividualSamples == 1',
                                plotly::plotlyOutput("tx_samples_tpm")),
               width = 3
      ),
      
      tabPanel("Heatmap of transcript expression",
               plotly::plotlyOutput("heatmap")
               
      )
    )
  )
)
)


server <- function(input, output, session) {
  
  
  # Define reactives ----
  my_gene_id <- reactive({input$gene_id_or_name %>%
      split_text_to_vector() %>%
      convert_to_wb_id(gids)})
  my_gene_id_1 <- reactive(my_gene_id()[[1]])
  my_gene_name_1 <- reactive(i2s(my_gene_id_1(), gids))
  my_neurons <- reactive(input$neurons %>%
                           stringr::str_to_upper() %>%
                           split_text_to_vector() %>%
                           validate_neurons(neurons_table))
  
  
  scaleFill <- reactive(accepted_color_scales$fill_scale[accepted_color_scales$shortcode==input$useColorScale])
  scaleColor <- reactive(accepted_color_scales$color_scale[accepted_color_scales$shortcode==input$useColorScale])
  
  
  observeEvent(input$explainNeurons, {
    showModal(modalDialog(
      title = "Neurons",
      'Use ALL for all neurons, individual neuron names (e.g. "AWA", "ASEL", or "OLL"), or keywords such as "ACh", "motor", "sensory", ... The combination of all neurons corresponding to these keywords will be displayed.',
      easyClose = TRUE
    ))
  })
  
  # Single gene plots ----
  output$tx_neuron_proportions <- plotly::renderPlotly({
    
    ind_dta <- tx_long %>%
      filter(gene_id == !!my_gene_id_1(),
             neuron_id %in% !!my_neurons()) |>
      collect()
    
    ind_graph <- ind_dta %>%
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
      ggplot(aes(x=Neuron,y=`Transcript Usage`, fill = Transcript)) +
      geom_col(position = position_stack()) +
      # geom_errorbar(aes(ymin=prop,ymax=prop+sd), position = position_stack()) +
      theme_classic() +
      coord_flip() +
      labs(x=NULL,y="Mean proportion",title=paste0(my_gene_name_1(), "/", my_gene_id_1())) +
      scale_y_continuous(labels = pct) +
      scaleFill()
    
    plotly::ggplotly(ind_graph)
    
  })
  
  
  output$tx_neuron_tpm <- plotly::renderPlotly({
    
    
    ind_dta <- tx_long %>%
      filter(gene_id == !!my_gene_id_1(),
             neuron_id %in% !!my_neurons()) |>
      collect()
    
    ind_graph <- ind_dta %>%
      group_by(transcript_id, neuron_id) %>%
      summarize(`Mean TPM` = mean(TPM),
                sd = sd(TPM),
                .groups = "drop") %>%
      rename(`Neuron` = neuron_id,
             `Transcript` = transcript_id) %>%
      ggplot(aes(x=`Neuron`, y=`Mean TPM`, color = `Transcript`)) +
      geom_point(position=position_dodge(.2), pch=15,size=2) +
      geom_errorbar(aes(ymin=`Mean TPM`-sd,ymax=`Mean TPM`+sd),width=0,position=position_dodge(.2)) +
      theme_classic() +
      scaleColor() +
      coord_flip() +
      {if(input$log_scale) scale_y_log10()} +
      labs(x=NULL,y="Mean transcript TPM (+/-sd)",title=paste0(my_gene_name_1(), "/", my_gene_id_1()))
    
    plotly::ggplotly(ind_graph)
    
  })
  
  
  
  output$tx_samples_proportions <- plotly::renderPlotly({
    
    ind_dta <- tx_long %>%
      filter(gene_id == !!my_gene_id_1(),
             neuron_id %in% !!my_neurons()) |>
      collect()
    
    ind_graph <- ind_dta %>%
      group_by(sample_id, gene_id) %>%
      mutate(sample_prop = round(100*TPM/sum(TPM), 1)) %>%
      rename(Neuron = neuron_id,
             Transcript = transcript_id,
             `Transcript Proportion` = sample_prop) %>%
      ggplot(aes(x=sample_id, y=`Transcript Proportion`, fill = Transcript)) +
      geom_col(position = position_stack()) +
      # geom_errorbar(aes(ymin=prop,ymax=prop+sd), position = position_stack()) +
      theme_classic() +
      coord_flip() +
      labs(x=NULL,y="Proportion usage",title=paste0(my_gene_name_1(), "/", my_gene_id_1())) +
      scale_y_continuous(labels = pct) +
      scaleFill()
    
    plotly::ggplotly(ind_graph)
    
  })
  
  output$tx_samples_tpm <- plotly::renderPlotly({
    
    ind_dta <- tx_long %>%
      filter(gene_id == !!my_gene_id_1(),
             neuron_id %in% !!my_neurons()) %>%
      rename(Neuron = neuron_id,
             Transcript = transcript_id,
             `Transcript Usage` = TPM) |>
      collect()
    
    ind_graph <- ind_dta %>%
      ggplot(aes(x=sample_id, y=`Transcript Usage`, fill = Transcript)) +
      theme_classic() +
      geom_point(position=position_dodge(.2), pch=15,size=2) +
      coord_flip() +
      {if(input$log_scale) scale_y_log10()} +
      labs(x=NULL,y="Transcript TPM",title=paste0(my_gene_name_1(), "/", my_gene_id_1())) +
      scaleFill()
    
    plotly::ggplotly(ind_graph)
    
  })
  
  
  
  # Heatmap ----
  output$heatmap <- plotly::renderPlotly({
    
    plot_data <- tx_long %>%
      filter(gene_id %in% !!my_gene_id(),
             neuron_id %in% !!my_neurons()) |>
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
      mutate(Transcript = fct_inorder(Transcript)) %>%
      ggplot() +
      theme_minimal() +
      {if(input$log_scale) {
        geom_tile(aes(x = Neuron, y =  Transcript, fill = `log(TPM)`))
      } else {
        geom_tile(aes(x = Neuron, y =  Transcript, fill = TPM))}} +
      geom_text(aes(x = 'Gene', y = Transcript, color = Gene, label = Gene)) +
      scale_x_discrete(limits = c(unique(plot_data$Neuron), '','Gene',' ')) +
      # viridis::scale_fill_viridis(direction = 1, option = "magma") +
      scale_fill_gradientn(colors = MetBrewer::met.brewer("OKeeffe2",direction = -1)) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5)) +
      {if(input$log_scale) {
        labs(title = "Transcript expression level (log TPM — mean per neuron)", x=NULL, y=NULL)
      } else {
        labs(title = "Transcript expression level (mean TPM per neuron)", x=NULL, y=NULL)}}
    
    
    plotly::ggplotly(gg)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)




















