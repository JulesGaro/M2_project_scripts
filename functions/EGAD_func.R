# Function Use in the GBA analysis

load_go_annotation <- function(file_path = '/space/scratch/jgarreau/data/pro_GO.csv') {
  # Load Gene ontology annotation file
  # 
  # Args:
  #   file_path : annotation file path
  # 
  # Return:
  #   GO_unique_filtered : GO annotation tibble

  GO <- read.delim(file = file_path, sep = ",", stringsAsFactors = TRUE)
  GO_unique <- data.frame(table(GO$GO.ID))
  colnames(GO_unique) <- c('GO', 'count')
  GO_unique_filtered <-  filter(GO_unique, count >=20) %>% as.tibble()
  return(GO_unique_filtered)
}