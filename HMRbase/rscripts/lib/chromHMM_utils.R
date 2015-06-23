emission.load <- function(file.emi) {
  emi_matrix <- read.table(file.emi, header=TRUE, sep='\t',
                   colClasses=c('character',rep('numeric',2), 'character', 'numeric'))
  colnames(emi_matrix)[1] <- "state"
  state <- paste0(rep('E', dim(emi_matrix)[1]), emi_matrix$state)
  emi_matrix$state <- NULL
  emi_matrix <- as.matrix(sapply(emi_matrix, as.numeric))
  rownames(emi_matrix) <- state
  return(emi_matrix)
}

emission.find_core <- function(emi_matrix) {
  sum_emi <- rowSums(emi_matrix)
  return (rownames(emi_matrix)[which.max(sum_emi)])
}