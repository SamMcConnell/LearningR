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

pdom <- function(k,m,n){
  total <- k + m + n
  denom <- total * (total - 1)
  return((k*(k-1) + 2*k*m + 2*k*n + m*(m-1)*0.75 + m*n)/denom)
}

rfib <- function(n,k){
  if (n < 1) return(1)
  a <- 1
  b <- 0
  for (i in 1:(n-1)) {
    temp <- a
    a <- a + k*b
    b <- temp
  }
  return(a)
}

high.gc <- function(x){   #"FastaFile.txt" -> "ID ##.#####"  Produces ID and GC percentage of the sequence with the highest GC content (%)
  gccal <- function(s){   #Characters (only A, C, T, or G's) -> numeric.  Produces % GC content of sequence
    s1 = strsplit(s, split="")[[1]]
    return(100 * length(grep("[GC]", s1)) / length(s1))
  }
  fdf <- function(x){     #"FastaFile.txt" -> data.frame with Sequence IDs as row names.  All Sequences under Sequence column
    rsf = readLines(x)
    rsf = strsplit(gsub("([ACGT])\\n([ACGT])", "\\1\\2", paste0(rsf, collapse="\n")), split="\n")[[1]]
    ids = grepl("^>", rsf)
    return(data.frame(row.names = sub(">", "", rsf[ids]), Sequence = rsf[!ids], stringsAsFactors = F))
  }
  
  rsf = fdf(x)
  rsf[,"Sequence"] <- unlist (lapply (rsf[,"Sequence"], gccal))
  n = apply(rsf, 2, which.max)
  return(cat(row.names(rsf)[n], rsf[n,]))
}