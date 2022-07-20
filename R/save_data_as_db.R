# library(readr)
# 
# tx_long <- read_tsv("data/t_exp.tsv",
#                     col_types = cols(transcript_id = col_character(),
#                                      gene_id = col_character(),
#                                      sample_id = col_character(),
#                                      neuron_id = col_character(),
#                                      TPM = col_double()))
# stop_for_problems(tx_long)
# 
# # save as SQLite
# con <- DBI::dbConnect(RSQLite::SQLite(), "data/t_exp.sqlite.db")
# DBI::dbWriteTable(con, name = "t_exp",value = tx_long)
# DBI::dbDisconnect(con)
# 
# # save as DuckDB
# con <- DBI::dbConnect(duckdb::duckdb(), "data/t_exp.duckdb.db")
# DBI::dbWriteTable(con, name = "t_exp",value = tx_long)
