LOGCOLS_VDHMM <- c("ITR", "F_P", "F_P", "F_ALPHA","F_BETA", "B_ALPHA", "B_BETA", "r",
                  "LLH", "LLH_LAST_ITR", "LLH_DELTA")

LOGCOLS_HMM <- c("ITR", "F_P", "B_P", "F_ALPHA","F_BETA", "B_ALPHA", "B_BETA",
                   "LLH", "LLH_LAST_ITR", "LLH_DELTA")

LOGCOLS_CTHMM <- c("ITR", "F_RATE", "B_RATE", "F_ALPHA","F_BETA", "B_ALPHA", "B_BETA",
                 "LLH", "LLH_LAST_ITR", "LLH_DELTA")

get_params_from_log <- function(file, LOGCOLS) {
  con <- file(file, blocking = FALSE)
  lines <- readLines(con)
  elems <- strsplit(lines, "\\s+")
  elems <- sapply(elems, function(x)  x[x != ""])
  tab <- elems[sapply(elems, function(x) length(x) == length(LOGCOLS))]
  params <- matrix(as.numeric(tail(tab, n=1)[[1]]), nrow=1)
  colnames(params) <- LOGCOLS
  close(con)
  return(tab)
}
