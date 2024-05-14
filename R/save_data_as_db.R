# library(readr)
# library(dplyr)
# library(wbData)
# gcoords <- wb_load_gene_coords(289)
# 
# tx_long <- read_tsv("data/t_exp.tsv",
#                     col_types = cols(transcript_id = col_character(),
#                                      gene_id = col_character(),
#                                      sample_id = col_character(),
#                                      neuron_id = col_character(),
#                                      TPM = col_double()))
# stop_for_problems(tx_long)
# 
# # Add URL to Wormbase browser
# tx_long2 <- gcoords %>%
#   mutate(padding = round(0.15*(end-start)),
#          position = paste0(chr,":",
#                            pmax(0,start-padding),
#                            "..",
#                            end+padding),
#          browser_url = paste0("https://wormbase.org/tools/genome/jbrowse-simple/?loc=",position,"&data=data%2Fc_elegans_PRJNA13758&tracks=Curated_Genes%2CRNASeq%2CRNASeq%20Splice%20Junctions%20(common)")) %>%
#   select(gene_id, browser_url) %>%
#   right_join(tx_long,
#              by = "gene_id")
# 
# 
# 
# #
# # # save as SQLite
# con <- DBI::dbConnect(RSQLite::SQLite(), "data/t_exp.sqlite.db")
# DBI::dbWriteTable(con, name = "t_exp",value = tx_long2)
# DBI::dbDisconnect(con)
# 
# unique(tx_long$neuron_id) |>
#   qs::qsave("data/measured_neurons.qs")


# 
# # save as DuckDB
# con <- DBI::dbConnect(duckdb::duckdb(), "data/t_exp.duckdb.db")
# DBI::dbWriteTable(con, name = "t_exp",value = tx_long)
