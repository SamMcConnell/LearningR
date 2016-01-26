count.bases <- function(x){
  return(table(strsplit(x, split="")[[1]]))
}

transcribe <- function(x){
  x <- strsplit(x, split="")[[1]]
  x[x=="T"] <- "U"
  return(paste0(x, collapse=""))
}

revc <- function(x){ # Returns the reverse compliment of a string of DNA
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
  if (n < 2) return(1)
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
    rsf = readLines(x)    # rsf = result-so-far
    rsf = strsplit(gsub("([ACGT])\\n([ACGT])", "\\1\\2", paste0(rsf, collapse="\n")), split="\n")[[1]]
    ids = grepl("^>", rsf)
    return(data.frame(row.names = sub(">", "", rsf[ids]), Sequence = rsf[!ids], stringsAsFactors = F))
  }
  
  rsf = fdf(x)
  rsf[,"Sequence"] <- unlist (lapply (rsf[,"Sequence"], gccal))
  n = apply(rsf, 2, which.max)
  return(cat(row.names(rsf)[n], rsf[n,]))
}

translate <- function(x) { # Translates an RNA sequence to an amino acid sequence
  sst <- strsplit(x, split="")[[1]]
  off <- length(sst)%%3
  
  if (off != 0){
    sst <- head(sst, -(off))
  }
  
  sst <- paste0(sst[c(T,F,F)], sst[c(F,T,F)], sst[c(F,F,T)])
  
  codon <- function(x) { # Translates a codon to its appropriate amino acid
    aas <- c("GC[ACGU]", "UG[UC]", "GA[UC]", "GA[AG]", "UU[UC]", "GG[ACGU]", "CA[UC]", "AU[ACU]", "AA[AG]", "UU[AG]|CU[ACGU]", "AUG", "AA[UC]", "CC[ACGU]", "CA[AG]", "CG[ACGU]|AG[AG]", "UC[ACGU]|AG[UC]", "AC[ACGU]", "GU[ACGU]", "UGG", "UA[UC]", "UA[AG]|UGA")
    aas1 <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "")
    temp = as.logical(sapply(aas, regexpr, x) + 1)
    return(aas1[temp])
  }
  
  return(paste0(sapply(sst, codon), collapse=""))
}

hamm <- function(s, t) { # Returns the number of substitutions between two sequences of the same length
  s <- strsplit(s, split="")[[1]]
  t <- strsplit(t, split="")[[1]]
  r <- s != t
  return(length(s[r]))
}

subs <- function(x, y){ # Returns the positions the pattern X appears in Y, including overlapping positions (gregexpr does not account for overlapping patterns in Y).
  xlen <- nchar(x)
  ylen <- nchar(y)
  count <- numeric()
  count2 <- 0
  
  for (xlen in 1:ylen) {
  if (regexpr(x, y)[1] == -1) break
  count2 <- count2 + regexpr(x, y)[1]
  count <- append(count, count2)
  y <- substr(y, (regexpr(x, y)[1] + 1), nchar(y))
  ylen <- nchar(y)
  }
  
  return(count)
}

iev <- function(a, b, c, d, e, f){
  return(2*(a+b+c) + 1.5*d + e)
}

fibd <- function(n, m){ # Returns the population after n rounds of growth
  # Each round of growth takes one month.  Each organism dies after m months
  # An organism can only reproduce if it is older than 0 months
  # The number of offspring introduced in a month is equal to the total of all reproductively active organisms
  
  if (n == 1) return(1)
  age.dist = as.bigz(vector(mode = "integer", length = m)) # Age distribution over the ages of 0, 1, 2, ... m-2, m-1 months (m not included, as they would be dead)
  age.dist[1] = 1L    # All populations are assumed to start with 1 reproductively immature organism
  
  for (month in 1:(n-1)){
    offspring = sum(tail(age.dist, -1))
    age.dist = append(offspring, head(age.dist, -1))
  }
  
  return(sum(age.dist))
}

codon.count <- function(x) {
  x <- append(strsplit(x, split = "")[[1]], "stop")
  codons <- c(4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 1, 2, 4, 2, 6, 6, 4, 4, 1, 2, 3)
  amino.acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "stop")
  x <- codons[match(x, amino.acids)]
  x <- unlist(lapply(split(x, ceiling(seq_along(x)/14)), prod)) %% 1000000
  if(length(x)==1) return(unname(x))
  for (i in 2:length(x)) {
    x[i] <- (x[i]*x[i-1]) %% 1000000
    x[i-1] <- 0
  }
  return(unname(x[length(x)]))
}

