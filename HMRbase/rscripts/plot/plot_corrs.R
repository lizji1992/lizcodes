plot_corrs <- function(dir = "~/panfs/corr", name="Blood_1") {
  f_all <- file.path(dir, paste0(name, "_all.corr"))
  f_near <- file.path(dir, paste0(name, "_near.corr"))
  x_all <- read.table(f_all, header=FALSE, col.names=c("dist", "corr", "num"))
  x_near <- read.table(f_near, header=FALSE, col.names=c("dist", "corr", "num"))
  d_all <- cat_dfs_from_list(list(all=x_all, near=x_near))
  p_all <- ggplot(d_all, aes(x=dist, y=corr, color=type)) + geom_line(alpha=0.5)
  print(p_all)
  
  f_fg <- file.path(dir, paste0(name, "_fg.corr"))
  f_fg_notspan <- file.path(dir, paste0(name, "_fg_notspan.corr"))
  x_fg <- read.table(f_fg, header=FALSE, col.names=c("dist", "corr", "num"))
  x_fg_notspan <- read.table(f_fg_notspan, header=FALSE, col.names=c("dist", "corr", "num"))
  d_fg <- cat_dfs_from_list(list(fg=x_fg, fg_notspan=x_fg_notspan))
  p_fg <- ggplot(d_fg, aes(x=dist, y=corr, color=type)) + geom_line(alpha=0.5)
  print(p_fg)
  
  
  f_bg <- file.path(dir, paste0(name, "_bg.corr"))
  f_bg <- file.path(dir, paste0(name, "_bg_notspan.corr"))
  x_bg <- read.table(f_bg, header=FALSE, col.names=c("dist", "corr", "num"))
  x_bg_notspan <- read.table(f_bg_notspan, header=FALSE, col.names=c("dist", "corr", "num"))
  d_bg <- cat_dfs_from_list(list(bg=x_bg, bg_notspan=x_bg_notspan))
  p_bg <- ggplot(d_bg, aes(x=dist, y=corr, color=type)) + geom_line(alpha=0.5)
  print(p_bg)  
}
  