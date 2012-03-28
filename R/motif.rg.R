library(Biostrings)
library (seqLogo)

setClass("Motif", representation(score="numeric", sign="logical", count="numeric", match="data.frame", pattern="character", consensus="character"))

my.pvec <- function (v, FUN, simplify=F,  mc.cores,..., mc.set.seed = TRUE, mc.silent = FALSE, mc.cleanup = TRUE) 
{
    if (!is.vector(v)) 
        stop("'v' must be a vector")
    env <- parent.frame()
    cores <- as.integer(mc.cores)
    if (cores < 1L) 
        stop("'mc.cores' must be >= 1")
    if (cores == 1L) 
        return(FUN(v, ...))
    if (mc.set.seed) 
        mc.reset.stream()
    n <- length(v)
    l <- if (n <= cores) 
        as.list(v)
    else {
        il <- as.integer(n/cores)
        xc <- n - il * cores
        sl <- rep(il, cores)
        if (xc) 
            sl[1:xc] <- il + 1L
        si <- cumsum(c(1L, sl))
        se <- si + c(sl, 0L) - 1L
        lapply(seq_len(cores), function(ix) v[si[ix]:se[ix]])
    }
    jobs <- NULL
    ## cleanup <- function() {
    ##     if (length(jobs) && mc.cleanup) {
    ##         mccollect(children(jobs), FALSE)
    ##         mckill(children(jobs), if (is.integer(mc.cleanup)) 
    ##             mc.cleanup
    ##         else 15L)
    ##         mccollect(children(jobs))
    ##     }
    ##     if (length(jobs)) {
    ##         mccollect(children(jobs), FALSE)
    ##     }
    ## }
    ## on.exit(cleanup())
    FUN <- match.fun(FUN)
    jobs <- lapply(seq_len(min(n, cores)), function(i) mcparallel(FUN(l[[i]], 
        ...), name = i, mc.set.seed = mc.set.seed, silent = mc.silent))
    res <- mccollect(jobs)
    names(res) <- NULL
    if(simplify){
      res <- do.call(c, res)
      if (length(res) != n) 
        warning("some results may be missing, folded or caused an error")
    }
    res
}

polyn <- function(letter, ntimes)
  {
    paste(rep(letter, ntimes),collapse="")
  }

summaryMotif <- function(motifs, category){ 
  if(is.null(motifs)) {return(NULL)}
  cons <- c()
  scores <- c()
  signs <- c()
  hits.counts1 <- c()
  hits.counts2 <- c()
  seq.counts1 <- c()
  seq.counts2 <- c()  
  fg.set <- category == 1
  bg.set <- !fg.set
  fg.size <- sum(fg.set)
  bg.size <- sum(bg.set)
                                        #summarize motif 
  for(i in 1:length(motifs)){
    scores[i] <- motifs[[i]]@score
    signs[i] <- motifs[[i]]@sign    
    pwm <-  getPWM( (motifs[[i]]@match)$pattern, pseudo=0)
    cons[i] <- consensus(pwm$prob)
    count <- motifs[[i]]@count
    hits.counts1[i] <- sum(count[fg.set])
    hits.counts2[i] <- sum(count[bg.set])
    seq.counts1[i] <-  sum(count[fg.set] > 0)
    seq.counts2[i] <-  sum(count[bg.set] > 0)
  }
  ratio <- (hits.counts1/hits.counts2)/(fg.size /bg.size)
  frac1 <- seq.counts1/fg.size
  frac2 <- seq.counts2/bg.size
  summary <- data.frame(scores=scores,
                        signs= signs,
                        fg.hits=hits.counts1, bg.hits= hits.counts2,
                        fg.seq=seq.counts1, bg.seq = seq.counts2,
                        ratio=ratio, fg.frac=frac1, bg.frac=frac2)       
  row.names(summary) <- cons
  summary
}


patternInfo <- function(x){
  v <- unlist(strsplit(x, ""))
  sum(sapply(v, function(y){2-log2(length(IUPAC[[y]]))}))             
}

IUPAC <- lapply(IUPAC_CODE_MAP, function(x){unlist(strsplit(x, ""))})
IUPAC.mat <- sapply(names(IUPAC), function(x){
  sapply(names(IUPAC), function(y){
    all(IUPAC[[x]] %in% IUPAC[[y]])
  })})


findMotif <-function(all.seq,                                #all sequences for motif search
                     category,                               #binary vector: sequences are in foreground or background       
                     start.width=6,                          #the starting width for motif
                     min.cutoff=5,                           #Z-value threshold for motifs
                     min.ratio = 1.3,                        #Minimal fold change required for motifs.
                     min.frac = 0.01,                                              
                     both.strand=T,
                     flank=2,
                     max.motif = 5,                                        
                     mask=T,                                 #whether we should mask the previously discovered mtifs
                     other.data=NULL,                        #additional terms (formula and data) for regression
                     start.nmer= NULL,                       #matrix for precomputed nmer counts
                     enriched.only=F,                        #compute only the enriched motifs
                     n.bootstrap = 5,                        #bootstrapping times
                     bootstrap.pvalue=0.1,                   #bootstrapping pvalue
                     is.parallel = T,
                     mc.cores = 4,
                     min.info = 10,
                     max.width=15)  
  {
    if(is.parallel){
      if(Sys.info()[["sysname"]] == "Windows"){
        is.parallel=F
      }
      else{
        if(!require("parallel")){
          is.parallel=F
        }
      }
    }
    fg.set <- category==1
    bg.set <- category==0
    fg.num <- sum(fg.set)
    bg.num <- sum(bg.set)
    fg.size <- sum(width(all.seq[fg.set]))
    bg.size <- sum(width(all.seq[bg.set]))
    #partitions <- lapply(1:npartition, function(x){sample(1:length(all.seq), as.integer(length(all.seq)* partition.frac))})    
    
    filterByRatio <- function(fg.total, bg.total, nmer=names(fg.total))
      {
        frac1 <- fg.total /fg.num
        frac2 <- fg.total /bg.num      
        ratio <- (fg.total+1)/fg.size / ((bg.total+1)/bg.size)
        filter <- abs(log(ratio)) > log(min.ratio)  & pmax(frac1, frac2) >  min.frac
        as.character(nmer) %in% names(fg.total)[filter]
      }
    
                                        #enumerate all nmers
    if (is.null(start.nmer)){
      fg.total <- countNmer( all.seq[fg.set], start.width, collapse=T, both.strand=both.strand)
      bg.total <- countNmer( all.seq[bg.set], start.width, collapse=T, both.strand=both.strand)
    }
    else{
        fg.total <- countSetGeneral(start.nmer, all.seq[fg.set],collapse=T)
        bg.total <- countSetGeneral(start.nmer, all.seq[bg.set],collapse=T)
    }
    start.nmer <- names(fg.total)[filterByRatio(fg.total, bg.total)]
    if(is.null(start.nmer) | length(start.nmer)==0){
      return(NULL)
    }

    if(is.parallel){
      tmp1 <- my.pvec(as.character(start.nmer),mc.cores=mc.cores, simplify=F,function(m){
        counts <- countSet(m, all.seq, both.strand=both.strand)
        zvalue <- getScore(category, counts, colnames(counts) , other.data, discretize=T)
        list(counts, zvalue)
      })
      counts <- do.call("cbind", lapply(tmp1, function(x){x[[1]]}))
      zvalue <- do.call("c", lapply(tmp1, function(x){x[[2]]}))
    }
    else{
      counts <- countSet(start.nmer, all.seq, both.strand=both.strand)
      zvalue <- getScore(category, counts, colnames(counts) , other.data, discretize=T)
    }
    
    pattern.counts <- counts
    candidates <- intersect(names(zvalue)[abs(zvalue) > min.cutoff], start.nmer)
    if(both.strand){
      candidates <- addReversecomplement(candidates)
    }
    
    org.seq <- all.seq
    org.seq.view <- pasteSeq(org.seq)    
    extendPatternFlank <- function(pattern, flank)
      {
        left.extend <- max(flank - (nchar(pattern) - nchar(gsub("^N+", "", pattern))),0)
        right.extend <- max(flank- (nchar(pattern) - nchar(gsub("N+$", "", pattern))),0)
        left.flank <- right.flank <- NULL
        if(left.extend > 0){
          left.flank <- paste(degeneratePattern(polyn("N", left.extend), alphabet=DNA_ALPHABET[1:14]), pattern,
                              polyn("N", right.extend), sep="")         
        }
        if(right.extend > 0){
          right.flank <- paste(polyn("N",left.extend), pattern,
                               degeneratePattern(polyn("N", right.extend), alphabet=DNA_ALPHABET[1:14]), sep="")                                
        }
        return(c(left.flank, right.flank))
      }
                                        #Refine the selected motifs
    refinePattern<- function(init.zvalue,
                             trim.mask=flank)
      {                 
        good.motifs <- list()
        nmotif <- 0        
                                        #last regression scoring results
        result <- data.frame(zvalue = init.zvalue, sign = init.zvalue > 0, score=abs(init.zvalue))
        row.names(result) <- names(init.zvalue)
        all.seq <- pasteSeq(all.seq)        
        repeat{
                                        #find the highest scoring pattern          
          patterns <- row.names(result)       
                                        #Skip the pattern contained within a selected motif, or if the pattern score is too low.
          patterns <- patterns[patterns %in% candidates]          
          if(length(patterns)==0){
            break
          }
          pattern <- ""
          for(i in 1:length(patterns)){
            if (result[patterns[i],"sign"] | !enriched.only){
              bs.scores <- getScore.bootStrap(category, pattern.counts[,patterns[i]] , other.data=other.data, discretize=T,n=n.bootstrap,is.parallel=is.parallel, mc.cores=mc.cores)
              pval <- t.test(bs.scores)$p.value
              if(is.na(pval)){next}
              cat(patterns[i], pval, "\n")
              if(pval < bootstrap.pvalue){
                pattern <- patterns[i]
                score <- result[patterns[i], "score"]
                sign <- result[pattern, "sign"]
                break
              }
            }
          }
          if(pattern == ""){break}
          if(i > 1){
            cat("Skip pattern ", patterns[1:(i-1)],"\n");            
          }
          if(length(patterns)==1){
            patterns <- NULL
          }
          else{
            patterns <- patterns[(i+1):length(patterns)]
          }
          seed <- pattern          
          motif.counts <- pattern.counts[,pattern]                    
          nmotif <- nmotif + 1
          pattern.fg.total <- sum(motif.counts[fg.set])
          pattern.bg.total <- sum(motif.counts[bg.set])
          cat(" Refine ", pattern, score, ":", bs.scores, sign, pattern.fg.total, pattern.bg.total,
              sum(motif.counts[fg.set]> 0), sum(motif.counts[bg.set]> 0), "\n")
                    
          filterCandidates <- function(to.test, fg.total, bg.total)
            {
              if (sign){
                filter1 <- fg.total >pattern.fg.total | bg.total < pattern.bg.total              
              }
              else{
                filter1 <- fg.total <pattern.fg.total | bg.total > pattern.bg.total
              }              
              to.test.info <- sapply(as.character(to.test), patternInfo)
              filter2 <- to.test.info> min.info
              filter3 <- filterByRatio(fg.total, bg.total)
              to.test <- to.test[filter1 & filter2 & filter3]
              if(both.strand){
                to.test <- reduceReversecomplement(to.test)
              }
            }
          if(both.strand){
            tested <- addReversecomplement(pattern)
          }
          start.flag <- T
          repeat{
            for(flag in c("extend", "mutate")){
              nochange <- 1
              repeat{
                if(flag == "extend"){ ##Try extend the pattern first.
                  if (nchar(pattern) < max.width){
                    to.test <- extendPatternFlank(pattern, flank)
                  }
                  else{
                    break
                  }
                }              
                else{ ##search the neighborhood with hamming distance de.allowed                  
                  to.test <- unique(degeneratePattern(pattern, alphabet= DNA_ALPHABET[1:15]))                  
                }
                if(is.null(to.test)| length(to.test) ==0){break}
                to.test <- setdiff(DNAStringSet(to.test), tested)
                if(is.null(to.test)| length(to.test) ==0){break}                
                if(nchar(seed) <= width(to.test)[1]){ ###Make sure that the seed is included in the pattern. 
                  to.test <- to.test[vcountPattern(seed, to.test, fixed=F)>0]
                }
                if(is.null(to.test)| length(to.test) ==0){break}
                fg.total <- countSetGeneral(to.test, all.seq[fg.set],collapse=T)
                bg.total <- countSetGeneral(to.test, all.seq[bg.set],collapse=T)
                to.test <- tryCatch(filterCandidates(to.test, fg.total, bg.total), error=function(e){
                  print("Filter error")
                  print(e)
                  NULL
                })
                if(is.null(to.test)| length(to.test) ==0){break}
                if(is.parallel){
                  tmp1 <- my.pvec(as.character(to.test), mc.cores=mc.cores,function(m){
                    counts <- countSetGeneral(m, all.seq, both.strand=both.strand)
                    zvalue <- getScore(category, counts, colnames(counts) , other.data, discretize=T)
                    select <- (zvalue > 0) == sign & abs(zvalue) > score
                    tmp.result <- list(counts[,names(select)[select],drop=F], zvalue=zvalue[select])
                    rm(counts)
                    rm(zvalue)
                    gc()
                    tmp.result
                  }, simplify=F)
                  test.counts <- do.call("cbind", lapply(tmp1, function(x){x[[1]]}))
                  zvalue <- do.call("c", lapply(tmp1, function(x){x[[2]]}))
                  ord <- order(abs(zvalue), decreasing=T)
                  zvalue <- zvalue[ord]
                }
                else{
                  test.counts <-  countSetGeneral(to.test,all.seq)
                  zvalue <- getScore(category, test.counts, terms=as.character(to.test), other.data= other.data, discretize=T)
                  select <- (zvalue > 0) == sign & abs(zvalue) > score
                  test.counts <- test.counts[, names(select)[select], drop=F]
                  zvalue <- zvalue[select]
                }
                result <- data.frame(zvalue = zvalue, sign = zvalue > 0, score=abs(zvalue))
                if(is.null(result) || nrow(result)==0) {break}
                for(p in row.names(result)){
                  new.bs.scores <- getScore.bootStrap(category, test.counts[,p] , other.data=other.data, discretize=T,n=n.bootstrap, mc.cores=mc.cores)
                  if(t.test(new.bs.scores, bs.scores)$p.value < bootstrap.pvalue){
                    pattern= p
                    bs.scores <- new.bs.scores
                    score <- result[p, "score"]
                    motif.counts <- test.counts[,pattern]
                    pattern.fg.total <- sum(motif.counts[fg.set])
                    pattern.bg.total <- sum(motif.counts[bg.set])
                    tmp <- trimPattern(pattern)
                    if (countPattern(seed, tmp) > 0){
                      pattern <- tmp
                    }
                    else{
                      seed <- pattern
                    }                  
                    cat(pattern, score, ":", bs.scores,  pattern.fg.total, pattern.bg.total, "\n")
                    nochange <- 0
                    break
                  }
                }                
                if(both.strand){
                  tested <- c(tested, addReversecomplement(to.test))
                }
              }
              if(nochange & !start.flag){break}
              start.flag=F
            }
            if(nochange){break}
          }
          pattern <- trimPattern(pattern)
          pattern.match <- findPatternView(pattern, all.seq,  both.strand=both.strand, flank=flank, rm.dup=F)          
          #collect motifs statistics
          tmp <- do.call("data.frame", list(sapply(colnames(pattern.match), function(x)pattern.match[[x]], simplify=F), stringsAsFactors=F))
          motif <- new("Motif", score=score, sign=sign, count=motif.counts, match=tmp,pattern=pattern)
          good.motifs[[pattern]] <- motif
          cat("New motif: ", pattern, "\n")
          if(is.null(patterns) || length(patterns) == 0 || nmotif >= max.motif){
            break
          }          
          if(mask){#####Mask existing motif matches #########
            tmp <- subject(all.seq)
            r <- reduce(ranges(pattern.match)[[1]])
            cat("match range ", length(r), "\n")
            start <- start(r)
            end <- end(r)
            if(!is.null(trim.mask) && 2* trim.mask < start.width + 2*flank){
              start <- start + trim.mask
              end <- end - trim.mask
            }
	    print(Views(tmp, start[1:3], end[1:3]))
            masks(tmp) <- Mask(length(tmp), start, end)
            tmp <- injectHardMask(tmp, "+")            
                                        #re-count the patterns
            all.seq <- Views(tmp, start(all.seq), end(all.seq))
            pattern.counts <- countSet(patterns, all.seq, both.strand=both.strand)
          }
          else{
            Pattern.counts <- countSet(patterns, all.seq, both.strand=both.strand)
          }                  
          print("Rescore")
          zvalue <- getScore(category, pattern.counts, patterns, other.data=other.data, discretize=T)
          result <- data.frame(zvalue = zvalue, sign = zvalue > 0, score=abs(zvalue))
          row.names(result) <- names(zvalue)          
          print("Finished Rescore")
        }
        return(good.motifs)
      }    
    motifs <- refinePattern(zvalue)    
    ##Rescore the motifs on unmasked sequence####
    if(mask){
      mask.motifs <- motifs
      motifs <- sapply(mask.motifs, function(x){
        y <- x
        pattern <- y@pattern                                        #find Match in the original sequence
        pattern.match <- findPatternView(pattern, org.seq.view,  both.strand=both.strand, flank=flank, rm.dup=F)
        y@match <- do.call("data.frame", list(sapply(colnames(pattern.match), function(x)pattern.match[[x]], simplify=F), stringsAsFactors=F))
        tmp <- rep(0, length(all.seq))
        l <- table(pattern.match$seq.id)
        tmp[as.numeric(names(l))] <- l
        y@count <- tmp
        tmp <- matrix(tmp,ncol=1, dimnames=list(NULL, pattern))
        zvalue <- getScore(category, tmp, pattern, other.data= other.data, discretize=T)
        y@score <- abs(zvalue)
        y
      })
      return(list(motifs=motifs, mask.motifs = mask.motifs, category=category,other.data=other.data))
    }
    return(list(motifs=motifs, category=category,other.data=other.data))
  }


consensus <- function(pwm)
  {
    pattern <- c()
    nt <- row.names(pwm)
    for(i in 1:ncol(pwm)){
      set <- c()
      ord <- order(pwm[,i],decreasing=T)
      if(pwm[ord[1],i] > 0.6){
        pattern[i] <- nt[ord[1]]
      }
      else{
        if(sum(pwm[ord[1:2], i]) >0.8){
          set <- nt[ord[1:2]]
        }
        else if(sum(pwm[ord[1:3], i]) > 0.95){
          set <- nt[ord[1:3]]
        }
        else{
          set <- nt
        }
        pattern[i] <- names(IUPAC_CODE_MAP)[sapply(IUPAC_CODE_MAP, function(x){nchar(x)==length(set) && all(set %in% unlist(strsplit(x, "")))})]
      }
    }
    pattern <- paste(pattern, collapse="")
    #pattern <- gsub('N+$', '',pattern)
    #pattern <- gsub('^N+', '',pattern)
    return(pattern)
  }

plotPWM <- function(pwm, alpha, has.box=T,...)
  {
    alphabet = row.names(pwm)
    width <- ncol(pwm)
    if (!has.box){
      xaxt <-"n"
    }
    else{
      xaxt <- "s"
    }
    if(has.box){
      plot(1:width, rep(0, width), ylab="", xlab="", ylim=c(0,length(alphabet)+1),xaxt=xaxt, yaxt="n", type="n",...)
      axis(2, 1:length(alphabet), alphabet)
    }
    else{
      opar <- par(mai = c(0, 0, 0, 0))
      on.exit(par(opar))
      plot(1:width, rep(0, width), ylab="", xlab="", ylim=c(0,length(alphabet)+1),xaxt="n", yaxt="n", type="n",...)
    }
    
    v <- rep(0, length(alphabet))
    names(v) <- alphabet
    psize <- 5* sqrt(pwm)
    #cols <- c("red", "blue", "green", "brown")
    cols <- c("green", "blue", "orange", "red")
    for(j in 1:ncol(pwm)){     
      for(i in 1:length(alphabet)){    
        tmp <- col2rgb(cols[i])/255
        if(alpha[j] > 0.1){
          col <- rgb(tmp[1], tmp[2], tmp[3], alpha=alpha[j])
          points(j, i, cex= psize[i,j],pch=alphabet[i],col=col)
        }
      }
    }
  }
    

plotPair <- function(i,j, joint.prob,  entropy=F, logodds=F, bg.ld=NULL)
  {
    alpha <- joint.prob
    if(logodds){
      bg.prob <- t(t(rowSums(joint.prob))) %*% colSums(joint.prob)
      ld <- log2(joint.prob/ bg.prob)
      if(!is.null(bg.ld)){
        ld <- ld - bg.ld
      }
      alpha <- pmax(joint.prob, bg.prob)            
    }
    if(entropy){
      alpha <-alpha *  ((4 - entropy(as.table(joint.prob)))/4)^(1/3)
    }
    for(k in 1:nrow(joint.prob))
      for(l in 1:ncol(joint.prob))
        if(alpha[k,l] > 0.05){
          if(!logodds){
            col <- rgb(0,0,0, alpha=alpha[k,l])
          }
          else{
            r<-b <- 0
            if (ld[k,l] > 0){
              b <- min(ld[k,l]*2,1)
            }
            if(ld[k,l] < 0){
              r <- min(-ld[k,l]*2,1)
            }              
            col <- rgb(r,0,b, alpha=alpha[k,l])            
          }
          lines(c(i,j), c(k, l), col=col,lwd=6)
        }
  }

plotMotif <-function(match, logodds=F, entropy=F, bg.ld=NULL, alphabet=c("A", "C", "G", "T"), has.box=T,...)
{
  match<- DNAStringSet(match)
  pwm = consensusMatrix(match)[alphabet,]
  pwm <- pwm/colSums(pwm)
  width <- ncol(pwm)
  if(entropy){
    alpha <- apply(pwm, 2, function(x){sqrt((2 - entropy(x))/2)})
  }
  else{
    alpha <- rep(1, ncol(pwm))
  }  
  plotPWM(pwm, alpha=alpha, has.box=has.box,...)
  for(i in 1:(width-1)){
    di <- suppressWarnings(nucleotideFrequencyAt(match, at=c(i, i+1), as.prob=T))
    plotPair(i, i+1, di,  entropy=entropy, logodds=logodds, bg.ld=bg.ld)
  }
}

motifScore <- function(motifs, seqs, category, other.data=NULL)
  {
    counts <-  rowSums(countSetGeneral(motifs,seqs))
    counts <- as.matrix(counts)
    colnames(counts) <- "motif"
    getScore(category, counts, colnames(counts) , other.data, discretize=T)    
  }

getScore <- function(response, data, terms=colnames(data), other.data=NULL, discretize=T, noutlier=6, family=binomial("logit"),sorted=T)
  {
    if(is.vector(data)){
      data <- matrix(data, ncol=1)
      terms <- 1
    }
    if(length(terms)==0){
      return(NULL);
    }
    if(is.null(other.data)){
      mat <- matrix(0, nrow=nrow(data), ncol=2)
    }
    else{
      l <- ifelse(is.vector(other.data), 1, ncol(other.data))
      mat <- matrix(0, nrow=nrow(data), ncol=2 + l)
      mat[,3:ncol(mat)] <- other.data
    }
    mat[,1] <- 1
    if(discretize){
      df <- cbind(response, as.data.frame(mat))
    }
    zvalues <- do.call("c", lapply(terms, function(v){
    #### create design matrix
      mat[,2] <- data[,v]     
      #If most of the data is the same value, skip regression
      tmp <- quantile(data[,v], c(noutlier/nrow(data), 1 - noutlier/nrow(data)))
      if(tmp[2] - tmp[1] <= 0){
        return(0)
      }            
      #### logistic regression
      if(discretize){
        df[,3] <- data[,v]
        tmp <- as.data.frame(table(df))
        tab.response <- as.numeric(levels(tmp[,1]))[tmp[,1]]
        tab.mat <- matrix(0, nrow=nrow(tmp), ncol=ncol(mat))                         
        for(i in 2:(ncol(tmp)-1)){
          #tab.mat[,i-1] <- as.numeric(levels(tmp[,i]))[tmp[,i]]
          tab.mat[,i-1] <- as.numeric(tmp[,i])
        }
        colnames(tab.mat) <- colnames(mat)
        weights <- tmp[,ncol(tmp)]        
        fit <- glm.fit(tab.mat, tab.response, family=family, weights=weights, intercept=T)                
      }
      else{
        fit <- glm.fit(mat, response, family=family, intercept=T)
      }
      ss <- summary.glm(fit)
      z <- ss$coeff[2,3]
      if(abs(z) > 1000 && discretize){
         fit <- glm.fit(mat, response, family=family, intercept=T)
         ss <- summary.glm(fit)
         z <- ss$coeff[2,3]
        if(abs(z) > 10000){
          #Glm error ignore the values
          z <- 0
        }
      }
      z
    }))
    rm(mat)
    gc()
    names(zvalues) <- terms
    if(sorted){
      ord <- order(abs(zvalues),decreasing=T)
      zvalues[ord]
    }
    else{
      zvalues
    }
  }

motifLatexTable <- function(motifs, main="", prefix="motif", dir=".", height=1, width=3,
                            enriched.only=F, plot.pwm=F,
                            summary.cols=c(1,7,8,9),use.mask=T)
  {
    if(!file.exists(dir)){
      dir.create(dir)
    }
    cat("\\begin{table}[ht]\n")
    cat("\\caption{", main, "}\n")
    cat("\\centering\n")
    #if (plot.pwm){
    #  require("seqLogo")
    #}
    if(use.mask){      
      summ <- summaryMotif(motifs$mask.motifs, motifs$category)
      motifs <- motifs$mask.motifs
    }
    else{
      summ <- summaryMotif(motifs$motifs, motifs$category)
      motifs <- motifs$motifs
    }
    ord <- order(summ$sign, abs(summ$score),decreasing=T)
    if(enriched.only){
      ord <- ord[summ$sign[ord]]
    }
    summ[!summ[,"signs"],"scores"] <- - summ[!summ[,"signs"],"scores"]
    motifs <- motifs[ord]
    summ <- summ[ord, summary.cols]
    summ <- format(summ, digits=2)
    nc <- paste(rep("l", ncol(summ) + 2), collapse="|")
    cat("\\begin{tabular}{|", nc, "|}\n",sep="")
    cat("\\hline\n")
    cat("Consensus", colnames(summ), sep="&")    
    cat("&logo\\\\\n")
    cat("\\hline\n")
    for(i in 1:nrow(summ)){
      cat(row.names(summ)[i], "&")
      cat(paste(summ[i,], collapse="&"))
      logofile <- paste(dir, paste(paste(prefix, i, sep="-"), "pdf", sep="."), sep="/")      
      pdf(logofile, height=height, width=width)
      par(mar=c(1,1,1,1))
      if(plot.pwm){
        seqLogo(getPWM(motifs[[i]]@match$pattern, pseudo=0)$prob, yaxis=F, xaxis=F)
      }
      else{
        plotMotif(motifs[[i]]@match$pattern,has.box=F)
      }
      dev.off()
      
      cat("&\\includegraphics[height=",height, "in, width=", width, "in]{", logofile, "}\\\\\n", sep="")
      cat("\\hline\n")
    }
    cat("\\end{tabular}\n")
    cat("\\end{table}\n")    
  }


motifHtmlTable <- function(motifs, dir="html", prefix="motif", enriched.only=F,plot.pwm=F, summary.cols=c(1,7,8,9),use.mask=T)
  {
    if(!file.exists(dir)){
      dir.create(dir)
    }
    header <- system.file("css", "header.html", package="motifRG")
    system(paste("cp", header, dir))
    outfile <- file.path(dir, paste(prefix, "html",sep=""))
    tmp <- paste("cp ",header, outfile)
    print(tmp)
    system(tmp)
    conn<-file(outfile, open="at")
    cat("<body>",file=conn)    
    #if (plot.pwm){
    #  require("seqLogo")
    #}
    if(use.mask){      
      summ <- summaryMotif(motifs$mask.motifs, motifs$category)
      motifs <- motifs$mask.motifs
    }
    else{
      summ <- summaryMotif(motifs$motifs, motifs$category)
      motifs <- motifs$motifs
    }
    ord <- order(summ$sign, abs(summ$score),decreasing=T)
    if(enriched.only){
      ord <- ord[summ$sign[ord]]
    }
    summ[!summ[,"signs"],"scores"] <- - summ[!summ[,"signs"],"scores"]
    motifs <- motifs[ord]
    summ <- summ[ord, summary.cols]
    summ <- format(summ, digits=2)
    nc <- paste(rep("l", ncol(summ) + 2), collapse="|")
    cat("<table border=\"1\" cellpadding=\"3\">\n",file=conn)
    cat("<tr>",file=conn)
    cat("<td> </td>",file=conn)
    for(x in colnames(summ)){
      cat("<th><b ID=\"reallyBold\">", x, "</b> </th>\n",file=conn)
    }
    cat("<th> logo </th>", file=conn)
    cat("</tr>",file=conn)
    for(i in 1:nrow(summ)){
      cat("<tr>",file=conn)
      cat("<td>", row.names(summ)[i], "</td>",file=conn)
      for(x in summ[i,]){
        cat("<td>", x, "</td>",file=conn)
      }
      logofile <- paste(paste(prefix, i, sep="-"), "png", sep=".")
      logofile.fp <- file.path(dir, logofile)
      png(logofile.fp, width=300, height=150)
      par(mar=c(1,1,1,1))
      if(plot.pwm){
        seqLogo(getPWM(motifs[[i]]@match$pattern, pseudo=0)$prob, yaxis=F, xaxis=F)
      }
      else{
        plotMotif(motifs[[i]]@match$pattern,has.box=F)
      }
      dev.off()      
      cat("<td> <table cellpadding=0><img src=\"", logofile, "\"></table></td>\n",file=conn)
      cat("</tr>",file=conn)
    }
    cat("</table>",file=conn)
    cat("</body>",file=conn)
    cat("</html>",file=conn)
    close(conn)
  }




getScore.bootStrap <- function(response, data, other.data=NULL, n = 5, is.parallel=T, mc.cores=mc.cores,...)
{
  tmp.fun <- function(i){
    x <- sample(1:length(response), length(response), replace=T)
    if(!is.null(other.data)){
      select.other.data <- other.data[x,]
    }
    else{
      select.other.data <- NULL
    }
    if(is.vector(data)){
      getScore(response[x], data[x], other.data=select.other.data,sorted=F,...)
    }
    else{
      getScore(response[x], data[x,], other.data=select.other.data,sorted=F,...)
    }
  }
  if(is.parallel){
    if (.Platform$OS.type == "windows") mc.cores <- 1
    result <- mclapply(1:n, tmp.fun,  mc.cores=mc.cores)
    if(!is.vector(data)){
      result <- do.call("cbind", result)
    }
    else{
      result <- do.call("c", result)
    }
  }
  else{
    result <- sapply(1:n, tmp.fun)
  }
}



