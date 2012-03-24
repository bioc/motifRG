library(Biostrings)
entropy <- function(counts)
  {
    freq <- counts/sum(counts)
    sum(-freq * log2(freq),na.rm=T)
  }

relativeEntropy <- function(counts, null)
  {
    freq <- counts/sum(counts)
    sum(freq * log2(freq/null),na.rm=T)    
  }

reduceReversecomplement<- function(nmers)
{
  nmers <- DNAStringSet(nmers)
  nmers.rev <- reverseComplement(nmers)
  to.remove <- nmers.rev %in% nmers & as.character(nmers.rev) < as.character(nmers)
  return(nmers[!to.remove])
}

addReversecomplement<- function(nmers)
  {
    nmers <- DNAStringSet(nmers)
    nmers.rev <- reverseComplement(nmers)
    return(union(nmers, nmers.rev))
  }


combineReversecomplement<- function(data, names=NULL)
{
  if(is.matrix(data)){
    pattern=colnames(data)
  }
  else{
    pattern=names(data)
  }    
  rev <- as.character(reverseComplement(DNAStringSet(pattern)))
  select <- which(rev  %in% pattern & pattern != rev)
  
  if(is.matrix(data)){
    for(i in select){
      data[,i] <- data[,i] + data[,rev[i]]
    }
    if (is.null(names)){
      data <- data[, pattern <= rev,drop=F]
    }
    else{
      data <- data[,pattern %in% names,drop=F]
    }
  }
  else{
    data[select] <- data[select] + data[rev[select]]
    if (is.null(names)){
      data <- data[pattern <= rev,drop=F]
    }
    else{
      data <- data[pattern %in% names,drop=F]
    }
  }
  return(data)
}

pasteSeq <- function(seqs, linker=10)
{
  len = nchar(seqs)
  delim = paste( rep("+", linker), collapse="")
  new.seq <- DNAString(paste(as.character(seqs), collapse=delim))     
  start = c(1, cumsum(len+linker)+1)[1:length(len)]
  end = start + len - 1
  v <- Views(new.seq, start=start, end=end)
  return(v)
}

trimPattern <- function(x) #remove the leading or trailing "N"
  {
    return(gsub("^N*(.*[^N]+)N*$", "\\1", x))
  }

countNmer <- function(seq, n, both.strand=T, collapse=F)
  {
    simplify.as <- "matrix"
    if(collapse){
      simplify.as <- "collapse"
    }
    tmp <- oligonucleotideFrequency(seq, n, simplify.as=simplify.as)
    if(both.strand){
      tmp <- combineReversecomplement(tmp)      
    }
    tmp
  }


countSet <- function(motifset, seqs, both.strand=T, collapse=F, extended=F,...)
  {
    orgset <- motifset
    motifset <- DNAStringSet(motifset)
    seqs <- DNAStringSet(seqs)
    if(collapse){
      collapse=1
    }
    if(both.strand){
      motifset <- addReversecomplement(motifset)
    }
    if(!extended){
      dict <- PDict(motifset)      
      counts <- vcountPDict(dict, seqs, collapse=collapse,...)
    }
    else{
      counts <- vcountPDict(motifset, seqs, collapse=collapse, fixed=F,...)
    }
    if(!collapse){
      if(is.matrix(counts)){
        counts <- t(counts)
      }
      else{
        counts <- as.matrix(counts)
      }
      colnames(counts) <- as.character(motifset)
    }
    else{      
      names(counts) <- as.character(motifset)
    }
    if(both.strand){
      counts <- combineReversecomplement(counts,orgset)
    }
    return(counts)
  }



##motifset can be different length, allow degenerate letters
countSetGeneral <- function(motifset, seqs, collapse=F, ...)
  {
    extended <- 1:length(motifset) %in% grep("[^ACGT]", as.character(motifset))
    counts <- NULL
    if(sum(!extended)>0){
      counts <- countSet(motifset[!extended], seqs,collapse=collapse,...)
    }
    if(sum(extended)>0){
      if(!collapse){
        counts <- cbind(counts, countSet(motifset[extended], seqs,  extended=T,collapse=collapse, ...))
      }
      else{
        counts <- c(counts, countSet(motifset[extended], seqs, extended=T,collapse=collapse,...))
      }
    }
    counts
  }


enumerateNmer <- function(width, alphabet = c("A", "C", "G", "T"))
  {
    pattern <- c()
    for(i in 1:width){
      pattern <- unlist(lapply(alphabet, function(x){paste(pattern, x, sep="")}))
    }
    return(pattern)
  }


appendNucleotides <- function(patterns, left=1, right=0)
  {
    if(left > 0){
      left.patterns <- enumerateNmer(left)
      patterns <- paste(rep(left.patterns, length(patterns)),
                        rep(patterns, rep(length(left.patterns), length(patterns))), sep="")      
    }
    if(right > 0){
      right.patterns <- enumerateNmer(right)
      patterns <- paste(rep(patterns, rep(length(right.patterns), length(patterns))),                        
                        rep(right.patterns, length(patterns)),
                        sep="")      
    }
    return(patterns)
  }

##detect pattern1 is a subpattern of pattern2 (with the same length)
isContained <- function(pattern1, pattern2)
{
  v1 <- unlist(strsplit(pattern1, ""))
  v2 <- unlist(strsplit(pattern2,""))
  tmp <- lapply(IUPAC_CODE_MAP, function(x){unlist(strsplit(x, ""))})
  for(i in 1:length(v1)){
    if( any(!tmp[[v1[i]]] %in% tmp[[v2[i]]])){
      return(FALSE)
    }
  }
  return(TRUE)
}

expandPattern<- function(patterns)
  {
    final.pattern <- c()
    for(pattern in patterns){
      all.pattern <- pattern
      for(i in 1:nchar(pattern)){
        nt <- substr(pattern,i,i)
        if (nt %in% names(IUPAC_CODE_MAP)){
          dest <- unlist(strsplit(IUPAC_CODE_MAP[[nt]], ""))
          all.pattern <- unlist(lapply(dest, function(x){substr(all.pattern, i,i) <- x
                                                         all.pattern}))
        }
      }
      final.pattern <- c(final.pattern,all.pattern)
    }    
    return(final.pattern)
  }


degeneratePattern <- function(patterns, hamming.distance=1, alphabet=DNA_BASES)
  {
    to.test <- c(patterns)
    for(i in 1:hamming.distance){
      temp <- c()
      for(j in 1:nchar(patterns[1])){
        for (nt in alphabet){
          temp <- c(temp,unlist(lapply(to.test, function(x){substr(to.test, j, j) <- nt
                                                            to.test})))
        }
      }
      to.test <- unique(temp)
    }
    return(to.test)
  }



findPatternView <- function(patterns, seqs.view,  both.strand=T, flank=0, rm.dup=T, ...)
{
  seq <- subject(seqs.view)
  if(length(patterns) > 1){
    dict <- PDict(patterns)
    match <- matchPDict(dict, seq,...)
    match.ranges <- IRanges(unlist(startIndex(match)) - flank,
                            unlist(endIndex(match))+ flank)
    match.view <- Views(seq, match.ranges)
    match.df <- data.frame(match.strand=rep("+", length(match.view)), pattern = as.character(match.view), stringsAsFactors =F)
  }
  else{
    match <- matchPattern(patterns, seq,fixed=c(F,T),...);
    match.ranges <- as(match, "IRanges")
    start(match.ranges) <- start(match.ranges) - flank
    end(match.ranges) <- end(match.ranges) + flank
    match.view <- Views(seq, match.ranges)    
    match.df <- data.frame(match.strand=rep("+", length(match.view)), pattern=as.character(match.view),  stringsAsFactors =F)
  }
  if(both.strand){        
    if(length(patterns) > 1){
      dict <- PDict(reverseComplement(DNAStringSet(patterns)))
      match <- matchPDict(dict, seq)
      tmp.ranges <- IRanges(unlist(startIndex(match)) - flank,
                            unlist(endIndex(match)) + flank)
      tmp <- Views(seq, tmp.ranges)
      match.df.rev <- data.frame(match.strand=rep("-", length(tmp)),
                                 pattern = as.character(reverseComplement(DNAStringSet(tmp))),
                                 stringsAsFactors =F)
      match.df <- rbind(match.df, match.df.rev)
      match.ranges <- append(match.ranges, tmp.ranges)
    }
    else{
      patterns <- reverseComplement(DNAString(patterns))
      match <- matchPattern(patterns, seq,fixed=c(F,T),...)
      tmp.ranges <- as(match, "IRanges")
      start(tmp.ranges) <- start(tmp.ranges) - flank
      end(tmp.ranges) <- end(tmp.ranges) + flank
      tmp <- Views(seq, tmp.ranges)    
      match.df <- rbind(match.df,  data.frame(match.strand=rep("-", length(tmp)),
                                              pattern=as.character(reverseComplement(DNAStringSet(tmp))),
                                              stringsAsFactors =F))
      match.ranges <- append(match.ranges, tmp.ranges)          
    }
  }
  ord <- order(start(match.ranges), match.df$match.strand)
  match.ranges <- match.ranges[ord]
  match.df <- match.df[ord,]
  if(rm.dup){
    dup <- duplicated(start(match.ranges))
    match.ranges <- match.ranges[!dup]
    match.df <- match.df[!dup,]
  }
  l <- findOverlaps(match.ranges, IRanges(start(seqs.view),end(seqs.view)), select="first")
  match.df$seq.id <- l
  match.df$pos <- start(match.ranges)-start(seqs.view)[l]
  return(RangedData(match.ranges, match.df))
}

isPalindrome <- function(x)
  {
    a <- DNAStringSet(unique(x))
    b <- reverseComplement(a)
    x[x %in% as.character(a) == as.character(b)]
  }