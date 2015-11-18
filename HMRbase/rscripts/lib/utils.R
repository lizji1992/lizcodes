cat_dfs_from_list <- function(data) {
  df <- data[[1]]
  for (i in 2:length(data)) {
    df <- rbind(df, data[[i]])
  }
  df$type <- rep(names(data), lapply(data, nrow))
  return(df)
}

