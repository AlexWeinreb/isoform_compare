
# add the % sign on the axis (for values already between 1-100)
pct <- function(...) scales::number(..., suffix = "%")


# copy relevant functions from wbData
wb_load_gene_ids <- function(WS){
  
  # validate input
  dir_cache <- "."
  
  cached_file <- paste0("data/c_elegans.PRJNA13758.WS",WS,".geneIDs.txt.gz")
  
  
  
  # read data, note format changed between versions
  if(WS >= 264){
    geneIDs <- readr::read_csv(cached_file,
                               col_names = c("X", "gene_id","symbol",
                                             "sequence","status","biotype"),
                               col_types = "cccccc")
  } else if(WS >= 236){
    geneIDs <- readr::read_csv(cached_file,
                               col_names = c("X", "gene_id","symbol",
                                             "sequence","status"),
                               col_types = "ccccc")
  } else{
    geneIDs <- readr::read_csv(cached_file,
                               col_names = c("gene_id","symbol","sequence"),
                               col_types = "ccc")
  }
  
  geneIDs$name <- ifelse(!is.na(geneIDs$symbol), geneIDs$symbol, geneIDs$sequence)
  geneIDs
}

i2s <- function(gene_id, geneIDs, warn_missing = FALSE){
  res <- geneIDs$name[match(gene_id, geneIDs$gene_id, incomparables = NA)]
  if(warn_missing && any(is.na(res))){
    warning("i2s: ",sum(is.na(res))," gene IDs could not be converted. NA are returned.")
  }
  res
}
wb_seq2id <- function(seq_id, geneIDs, warn_missing = FALSE){
  res <- geneIDs$gene_id[match(seq_id, geneIDs$sequence, incomparables = NA)]
  if(warn_missing && any(is.na(res))){
    warning("wb_seq2symbol: ",sum(is.na(res))," sequence IDs could not be converted. NA are returned.")
  }
  res
}
s2i <- function(symbol, geneIDs, warn_missing = FALSE){
  res <- geneIDs$gene_id[match(symbol, geneIDs$name, incomparables = NA)]
  if(warn_missing && any(is.na(res))){
    warning("s2i: ",sum(is.na(res))," gene symbols could not be converted. NA are returned.")
  }
  res
}



# Convert vector of any gene format into WBGene ID
convert_to_wb_id <- function(gene_id_or_name, gids){
  is_wb_id <- stringr::str_detect(gene_id_or_name, "^WBGene[\\d]{8}$")
  is_seq_name <- stringr::str_detect(gene_id_or_name, "^[A-Z0-9.]+$")
  is_symbol <- stringr::str_detect(gene_id_or_name, "^[a-z]{3,4}-[\\d]{1,4}")
  
  # Check
  nb_class <- (is_wb_id + is_seq_name + is_symbol)
  multi_class <- gene_id_or_name[nb_class > 1]
  unclassed <- gene_id_or_name[nb_class == 0]
  
  if(length(multi_class) != 0L){
    warning("Genes were misclassified during name conversion: ", multi_class)
  }
  if(length(unclassed) != 0L){
    warning("Genes were not classified during name conversion: ", unclassed)
  }
  
  # create out
  out_gene_id <- character(length(gene_id_or_name))
  out_gene_id[is_wb_id] <- gene_id_or_name[is_wb_id]
  out_gene_id[is_seq_name] <- wb_seq2id(gene_id_or_name[is_seq_name], gids, warn_missing = TRUE)
  out_gene_id[is_symbol] <- s2i(gene_id_or_name[is_symbol], gids, warn_missing = TRUE)
  
  out_gene_id
}

# split gene or neuron list into vector
split_text_to_vector <- function(text){
  stringr::str_split(text,
            pattern = "[[:space:],;]")[[1]] %>%
    stringi::stri_remove_empty()
}

# ensure we have a list of valid neurons
validate_neurons <- function(neurs, neurons_table){
  
  # replace synonyms
  neurs <- recode(neurs,
                  SENS="SENSORY",
                  INTER="INTERNEURON",
                  INTERNEURONS="INTERNEURON",
                  MOTORNEURON="MOTOR",
                  MOTORNEURONS="MOTOR",
                  MOTONEURON="MOTOR",
                  MOTONEURONS="MOTOR",
                  PHARYNX="PHARYNGEAL",
                  PHA="PHARYNGEAL",
                  CILIA="CILIATED",
                  CHOLINERGIC="ACH",
                  ACHERGIC="ACH",
                  ACETYLCHOLINERGIC="ACH",
                  ACETYLCHOLINE="ACH",
                  CHO="ACH",
                  SEROTONERGIC="SEROTONIN",
                  SEROTONINERGIC="SEROTONIN",
                  `5-HT`="SEROTONIN",
                  `5HT`="SEROTONIN",
                  GLUTAMATERGIC="GLUTAMATE",
                  GLU="GLUTAMATE",
                  GLUERGIC="GLUTAMATE",
                  GLUTA="GLUTAMATE",
                  GABAERGIC="GABA",
                  OCTO="OCTOPAMINE",
                  OCTOPAMINERGIC="OCTOPAMINE")
  
  
  if(length(neurs) == 1L && stringr::str_to_upper(neurs) == "ALL"){
    return(neurons_table$Neuron_type)
  }
  
  
  # Subcategories
  if(any(neurs == "SENSORY")){
    neurs <- setdiff(neurs, "SENSORY") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "sensory")])
  }
  if(any(neurs == "INTERNEURON")){
    neurs <- setdiff(neurs, "INTERNEURON") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "interneuron")])
  }
  if(any(neurs == "MOTOR")){
    neurs <- setdiff(neurs, "MOTOR") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "motor")])
  }
  if(any(neurs == "PHARYNGEAL")){
    neurs <- setdiff(neurs, "PHARYNGEAL") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "pharyngeal")])
  }
  if(any(neurs == "CILIATED")){
    neurs <- setdiff(neurs, "CILIATED") %>%
      c(neurons_table$Neuron_type[neurons_table$Ciliated == "yes"])
  }
  if(any(neurs == "ACH")){
    neurs <- setdiff(neurs, "ACH") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "ach")])
  }
  if(any(neurs == "SEROTONIN")){
    neurs <- setdiff(neurs, "SEROTONIN") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "serotonin")])
  }
  if(any(neurs == "GLUTAMATE")){
    neurs <- setdiff(neurs, "GLUTAMATE") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "glutamate")])
  }
  if(any(neurs == "GABA")){
    neurs <- setdiff(neurs, "GABA") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "gaba")])
  }
  if(any(neurs == "OCTOPAMINE")){
    neurs <- setdiff(neurs, "OCTOPAMINE") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "octopamine")])
  }
  
  
  
  
  
  
  # Remove unknown neurons
  unknown_neurs <- setdiff(neurs, neurons_table$Neuron_type)
  
  if(length(unknown_neurs) != 0L){
    warning("Neuron name not recognized: ", unknown_neurs)
  }
  
  setdiff(neurs, unknown_neurs)
}
