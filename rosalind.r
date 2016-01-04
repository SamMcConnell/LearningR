count.bases <- function(x){
  return(table(strsplit(x, split="")[[1]]))
}

transcribe <- function(x){
  x <- strsplit(x, split="")[[1]]
  x[x=="T"] <- "U"
  return(paste0(x, collapse=""))
}

revc <- function(x){
  x <- rev(strsplit(x, split="")[[1]])
  x <- c("T", "G", "C", "A")[match(x, c("A", "C", "G", "T"))]
  return(paste0(x, collapse=""))
}