


# Initializations ----

#~ Packages ----
library(shiny)
library(plotly)
library(tidyverse)

library(wbData)
gids <- wb_load_gene_ids(281)

#~ Data ----
tx_long <- read_tsv("data/t_exp.tsv",
                    col_types = cols(transcript_id = col_character(),
                                     gene_id = col_character(),
                                     sample_id = col_character(),
                                     neuron_id = col_character(),
                                     TPM = col_double()))
stop_for_problems(tx_long)

neurons_table <- read_csv("data/neuron_properties.csv",
                          col_types = "cccc") %>%
  filter(Neuron_type %in% as.character(unique(tx_long$neuron_id))) %>%
  mutate(across(Modality:Neurotransmitter, str_to_lower))
stop_for_problems(neurons_table)


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
                                plotlyOutput("tx_neuron_proportions")),
               conditionalPanel('input.plotIndividualSamples == 0',
                                plotlyOutput("tx_neuron_tpm")),
               conditionalPanel('input.plotIndividualSamples == 1',
                                plotlyOutput("tx_samples_proportions")),
               conditionalPanel('input.plotIndividualSamples == 1',
                                plotlyOutput("tx_samples_tpm")),
               width = 3
      ),
      
      tabPanel("Heatmap of transcript expression",
               # plotOutput("heatmapCplx"),
               # plotOutput("heatmapPhmp", width = "90%"),
               plotlyOutput("heatmap")
               
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
                           str_to_upper() %>%
                           split_text_to_vector() %>%
                           validate_neurons(neurons_table))
  
  
  scaleFill <- reactive(accepted_color_scales$fill_scale[accepted_color_scales$shortcode==input$useColorScale])
  scaleColor <- reactive(accepted_color_scales$color_scale[accepted_color_scales$shortcode==input$useColorScale])
  
  
  # shinyBS::addPopover(session,
  #                     id = "neurons",
  #                     title = "Neurons",
  #                     content = '<p>Use ALL for all neurons, individual neuron names (e.g. "AWA", "ASEL", or "OLL"), or keywords such as "ACh", "motor", "sensory", ...</p><p>You can combine keywords and neuron names as needed.</p>',
  #                     placement = "right",
  #                     trigger = "click")
  
  observeEvent(input$explainNeurons, {
    showModal(modalDialog(
      title = "Neurons",
      'Use ALL for all neurons, individual neuron names (e.g. "AWA", "ASEL", or "OLL"), or keywords such as "ACh", "motor", "sensory", ...',
      easyClose = TRUE
    ))
  })
  
  # Single gene plots ----
  output$tx_neuron_proportions <- renderPlotly({
    
    ind_graph <- tx_long %>%
      filter(gene_id == my_gene_id_1(),
             neuron_id %in% my_neurons()) %>%
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
    
    ggplotly(ind_graph)
    
  })
  
  
  output$tx_neuron_tpm <- renderPlotly({
    
    
    ind_graph <- tx_long %>%
      filter(gene_id == my_gene_id_1(),
             neuron_id %in% my_neurons()) %>%
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
    
    ggplotly(ind_graph)
    
  })
  
  
  
  output$tx_samples_proportions <- renderPlotly({
    
    ind_graph <- tx_long %>%
      filter(gene_id == my_gene_id_1(),
             neuron_id %in% my_neurons()) %>%
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
    
    ggplotly(ind_graph)
    
  })
  
  output$tx_samples_tpm <- renderPlotly({
    
    
    ind_graph <- tx_long %>%
      filter(gene_id == my_gene_id_1(),
             neuron_id %in% my_neurons()) %>%
      rename(Neuron = neuron_id,
             Transcript = transcript_id,
             `Transcript Usage` = TPM) %>%
      ggplot(aes(x=sample_id, y=`Transcript Usage`, fill = Transcript)) +
      theme_classic() +
      geom_point(position=position_dodge(.2), pch=15,size=2) +
      coord_flip() +
      {if(input$log_scale) scale_y_log10()} +
      labs(x=NULL,y="Transcript TPM",title=paste0(my_gene_name_1(), "/", my_gene_id_1())) +
      scaleFill()
    
    ggplotly(ind_graph)
    
  })
  
  
  
  # Heatmaps ----
  # output$heatmapCplx <- renderPlot({
  #     
  #     message("Plot Cplx!")
  #     tx_long %>%
  #         filter(gene_id %in% my_gene_id(),
  #                neuron_id %in% my_neurons()) %>%
  #         mutate(gene_id = i2s(gene_id, gids)) %>%
  #         group_by(gene_id, transcript_id, neuron_id) %>%
  #         summarize(mean_tpm = mean(TPM)) %>%
  #         rename(TPM = mean_tpm,
  #                Transcript = transcript_id,
  #                Neuron = neuron_id,
  #                Gene = gene_id) %>%
  #         group_by(Gene) %>%
  #         tidyHeatmap::heatmap(Transcript, Neuron, TPM,
  #                              transform = log1p, .scale = "none",
  #                              cluster_columns=FALSE,
  #                              cluster_rows = FALSE,
  #                              column_names_gp = grid::gpar(fontsize = 9))
  #     
  #     
  # })
  
  # output$heatmapPhmp <- renderPlot({
  #   
  #   tx_long %>%
  #     filter(gene_id %in% my_gene_id(),
  #            neuron_id %in% my_neurons()) %>%
  #     mutate(gene_id = i2s(gene_id, gids)) %>%
  #     group_by(gene_id, transcript_id, neuron_id) %>%
  #     summarize(mean_tpm = log1p(mean(TPM)),
  #               .groups = "drop") %>%
  #     rename(`log(TPM)` = mean_tpm,
  #            Transcript = transcript_id,
  #            Neuron = neuron_id,
  #            Gene = gene_id) %>%
  #     mutate(Transcript = factor(Transcript, levels = sort(as.character(unique(Transcript))))) %>%
  #     arrange(Gene, Transcript, Neuron) %>%
  #     tidyheatmap::tidy_heatmap(rows = Transcript,
  #                               columns =  Neuron,
  #                               values =  `log(TPM)`,
  #                               annotation_row = Gene,
  #                               gaps_row = Gene,
  #                               fontsize = 9,
  #                               main = "Transcript expression level (log TPM — mean per neuron)")
  # })
  
  
  output$heatmap <- renderPlotly({
    
    plot_data <- tx_long %>%
      filter(gene_id %in% my_gene_id(),
             neuron_id %in% my_neurons()) %>%
      mutate(gene_id = i2s(gene_id, gids)) %>%
      group_by(gene_id, transcript_id, neuron_id) %>%
      summarize(`log(TPM)` = log1p(mean(TPM)),
                TPM = mean(TPM),
                .groups = "drop") %>%
      rename(Transcript = transcript_id,
             Neuron = neuron_id,
             Gene = gene_id) %>%
      mutate(Transcript = factor(Transcript, levels = sort(as.character(unique(Transcript))))) %>%
      arrange(Gene, Transcript, Neuron)
    
    gg <- ggplot(plot_data) +
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
    
    
    ggplotly(gg)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)




















