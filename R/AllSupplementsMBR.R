################################################################################
#################  Internal functions for RTNduals-methods  ###################
################################################################################
##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##it checks if the input regulatory elements are in the 'TNI' annotation
.checkRegel <- function(tni, regulatoryElements)
{
  #---check regulatoryElements
  tp <- sapply(colnames(tni@annotation), function(i) {
    sum(regulatoryElements%in%tni@annotation[, i])
  })
  colid <- names(tp[which.max(tp)])
  idx <- which(tni@annotation[, colid]%in%regulatoryElements)
  if(length(idx) < length(regulatoryElements)) {
    warning("Not all 'regulatory elements' are available in the 'TNI' annotation!",
            call.=FALSE)
  }
  regulatoryElements <- tni@annotation[idx,]
  idx <- match(rownames(regulatoryElements), tni@transcriptionFactors)
  if(length(idx)==0){
    tp <- paste("NOTE: no 'regulatory element' has been used to call",
                "regulons in the provided 'TNI'!")
    stop(tp, call.=FALSE)
  } else if(length(idx) < length(regulatoryElements)){
    tp <- paste("Not all input 'regulatory elements' have been used to call",
                "regulons in the provided 'TNI'!")
    warning(tp, call.=FALSE)
  }
  regulatoryElements <- tni@transcriptionFactors[idx]
  return (regulatoryElements)
  }

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##creates a matrix of regulons (tnet) for analysis (can be numeric or factors)
.regMatrix <- function(regulons, regel, targets, getNames = TRUE, factors = FALSE)
{
  xmat <- matrix(0, nrow=length(targets), ncol=length(regel))
  rownames(xmat) <- targets
  colnames(xmat) <- regel
  for(i in regel)
  {
    regs <- regulons[[i]]
    if(factors == TRUE) xmat[names(regs), i] <- as.character(regs)
    else {xmat[names(regs), i] <- regs}
  }
  if(getNames) colnames(xmat) <- names(regel)
  if (factors == TRUE) xmat <- as.data.frame(xmat, stringAsFactors=factors)
  return(xmat)
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##gets significant associations inferred for motifs using the 'prob' parameter
.motifsquantile <- function(regcor, prob)
{
  if(length(regcor)>1){
    qtmat <- .getquantiles(regcor)
  } else {
    qtmat <- regcor
    qtmat[] <- 1
  }
  lmat <- array(TRUE,dim=dim(qtmat), dimnames = dimnames(qtmat))
  coord <- which(lmat,arr.ind=TRUE)
  rnames <- rownames(lmat)[coord[, 1]]
  cnames <- colnames(lmat)[coord[, 2]]
  corvalues <- regcor[coord]
  qdat <- qtmat[coord]
  idx <- !is.na(corvalues)
  qdat <- data.frame(Regulon1=rnames[idx],Regulon2=cnames[idx],
                     R=corvalues[idx], Quantile=qdat[idx], 
                     stringsAsFactors=FALSE)
  rownames(qdat) <- paste(qdat$Regulon1, qdat$Regulon2, sep="~")
  #---
  qdat <- qdat[qdat$Quantile>prob,]
  return(qdat)
}

#internal function '.motifsquantile'
.getquantiles <- function(regcor)
{
  qtmat <- as.numeric(abs(regcor))
  n <- sum(!is.na(qtmat))
  if(n>100) n <- 100
  qtmat <- cut(qtmat, breaks=quantile(qtmat, (0:n)/n, na.rm=TRUE), 
               include.lowest=TRUE)
  qtmat <- as.numeric(qtmat)/n
  qtmat <- matrix(qtmat, ncol=ncol(regcor), nrow=nrow(regcor), 
                   dimnames=list(rownames(regcor), colnames(regcor)))
  return(qtmat)
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##it gets the jaccard information between the 'duals'
.jcOverlap <- function(.statlist, regs1, regs2, tnet1, tnet2)
{
  #-----
  regEle1 <- unique(.statlist[, "Regulon1"])
  regEle1 <- regs1[regEle1]
  
  regEle2 <- unique(.statlist[, "Regulon2"])
  regEle2 <- regs2[regEle2]
  
  #----(jaccard)
  jcall <- .jc.overlap(regEle1, regEle2, cbind(tnet1,tnet2), overlap="all")
  if(is.matrix(jcall)){
    tb <- as.matrix(.statlist[, c("Regulon1", "Regulon2")]) 
    jcinf <- apply(tb, 1, function(x) {
      reg1 <- x[1]
      reg2 <- x[2]
      reg1 <- regs1[reg1]
      reg2 <- regs2[reg2]
      jcall[reg2, reg1]
    })
  } else {
    jcinf <- jcall
  }
  .statlist <- cbind(.statlist, Jaccard.coefficient=jcinf)
  return(.statlist)
}

#internal function '.jcOverlap'
.jc.overlap <- function(regEle1, regEle2, tnet, 
                        overlap=c("all","agreement","disagreement")){
  overlap <- match.arg(overlap)
  if(overlap == "all"){
    tnet[tnet != 0] <- 1
    jc <- function(x, xmat){
      c <- x+xmat
      a <- colSums(c == 2)
      b <- colSums(c > 0)
      b[b == 0] <- 1
      a/b
    }
  } else {
    if(overlap == "agreement"){
      ov <- (+1)
    } else {
      ov <- (-1)
    }
    tnet [tnet > 0] <- 1; tnet[tnet < 0] <- -1
    jc <- function(x, xmat){
      c <- x*xmat
      a <- colSums(c == ov)
      c <- abs(x) + abs(xmat)
      b <- colSums(c != 0)
      b[b == 0] <- 1
      a/b
    }
  }
  amap <- apply(tnet[, regEle1, drop=FALSE], 2, jc, 
                xmat=tnet[, regEle2, drop=FALSE])
  return(amap)
}


##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##gets mutual information between 'duals' reg. elements 
##and the p-value when available
.getMI <- function(.statlist, tbmi1, regs1, regs2, size1, size2)
{
  tb <- as.matrix(.statlist[, c("Regulon1", "Regulon2")])
  mutinf <- apply(tb, 1, function(x){
    reg1 <- x[1]
    reg2 <- x[2]
    reg1 <- regs1[reg1]
    reg2 <- regs2[reg2]
    s1 <- size1[reg1]
    s2 <- size2[reg2]
    mi <- abs(tbmi1[reg2, reg1])
    cbind(mi, round(s1), round(s2))
  })
  .statlist <- cbind(.statlist, MI=mutinf[1, ], Size.Regulon1=mutinf[2, ],
                 Size.Regulon2=mutinf[3, ])
  tp <- c("Regulon1", "Size.Regulon1", "Regulon2", "Size.Regulon2", "MI", 
          "R", "Quantile")
  .statlist <- .statlist[, tp]
  .statlist <- .statlist[which(.statlist$MI != 0), ]
  return (.statlist)
}
.getPMI <- function(.statlist, pvmat, regs1, regs2, cutoff)
{
  tb <- as.matrix(.statlist[, c("Regulon1", "Regulon2")])
  pvinf <- apply(tb, 1, function(x){
    reg1 <- x[1]
    reg2 <- x[2]
    reg1 <- regs1[reg1]
    reg2 <- regs2[reg2]
    pvmat[reg2, reg1]
  })
  .statlist <- cbind(.statlist, MI.Adjusted.Pvalue=pvinf)
  .statlist$MI.Adjusted.Pvalue[.statlist$MI.Adjusted.Pvalue<cutoff] <- paste("<", cutoff, sep="")
  tp <- c("Regulon1", "Size.Regulon1", "Regulon2", "Size.Regulon2", "MI", 
          "MI.Adjusted.Pvalue", "R", "Quantile")
  .statlist <- .statlist[, tp]
  return (.statlist)
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##it takes the correlation between regulons
.tni.cor<-function(x, tnet, estimator="pearson",dg=0, asInteger=TRUE, 
                   mapAssignedAssociation=TRUE){
  tfs<-colnames(tnet)
  tar<-rownames(tnet)
  ids<-unique(c(tfs,setdiff(tar,tfs)))
  x=x[ids,]
  x=t(x)
  #--
  pcorm=cor(x[,tfs],x[,tar], method=estimator,use="complete.obs")
  if(asInteger){
    pcorm[pcorm<0]=-1
    pcorm[pcorm>0]=1
  }
  if(length(tfs)>1)diag(pcorm[,tfs])=dg
  #--
  pcorm<-t(pcorm)
  colnames(pcorm)<-tfs
  if(mapAssignedAssociation)pcorm[tnet==0]=0
  pcorm
}

##------------------------------------------------------------------------------
#internal function for 'mbrAssociation'
##This function takes two gene sets (i.e. regulons), a vector 
##containing the size of the gene universe, and compute the number of 
##genes expected to occur in both regulons, the actual observed overlap, 
##and the pvalue from a hypergeometric test.
.mbr.hyper <- function(.statlist, regulons1, regulons2, regs1, regs2,
                       universe, pAdjustMethod, verbose=TRUE)
{
  regpairs <- .statlist[, c("Regulon1", "Regulon2")]
  if(verbose) pb<-txtProgressBar(style=3)
  res <- NULL
  for(i in 1:nrow(regpairs)) {
    if(verbose) setTxtProgressBar(pb, i/nrow(regpairs))
    vecpairs <- as.character(regpairs[i, ])
    ##---
    reg1 <- vecpairs[1]
    reg1 <- regs1[reg1]
    reg2 <- vecpairs[2]
    reg2 <- regs2[reg2]
    ##---
    regulon1 <- names(regulons1[[reg1]])
    regulon2 <- names(regulons2[[reg2]])
    tmp <- .regulon.hyper(regulon1=regulon1, universe=universe, 
                          regulon2=regulon2)
    res <- rbind(res, tmp)
  }
  if(verbose) close(pb)
  results <- cbind(regpairs, res)
  adjPvals <- p.adjust(results[, "Pvalue"], method = pAdjustMethod)
  results <- cbind(results, adjPvals)
  colnames(results)[ncol(results)] <- "Adjusted.Pvalue"
  results <- results[order(results[, "Pvalue"]), , drop=FALSE]
}

#internal function for '.mbr.hyper'
##it makes the hypergeometric test.
.regulon.hyper <- function(regulon1, universe, regulon2) {
  ##number of genes in universe
  N <- length(universe)			
  ##remove genes from gene set that are not in universe			
  regulon1 <- intersect(regulon1, universe) 
  regulon2 <- intersect(regulon2, universe)
  ##size of gene set	
  m <- length(regulon1) 							
  Nm <- N-m	
  ##regulon2 in gene set
  overlap <- intersect(regulon1, regulon2) 	
  ##number of hits between regulons		
  k <- length(overlap) 							
  n <- length(regulon2)	
  HGTresults <- phyper(k-1, m, Nm, n, lower.tail = FALSE)
  ex <- (n/N)*m
  if(m == 0 | n == 0) HGTresults <- 1
  hyp.vec <- c(N, m, n, ex, k, HGTresults)
  names(hyp.vec) <- c("Universe.Size", "R1.Size", "R2.Size", 
                      "Expected.Overlap", "Observed.Overlap", "Pvalue")
  return(hyp.vec)
}
##------------------------------------------------------------------------------
#internal function for 'tni2mbrPreprocess'
##it checks two 'TNIs' objects produced separetely
.checkTNIsProcessing <- function (tni1, tni2, verbose = TRUE)
{
  ## checks gexp consistency
  if (verbose)
    cat("-Checking expression matrix consistency...\n")
  if (!identical(tni1@gexp,tni2@gexp))
    stop("NOTE: TNIs should use the same expression matrix.")
  
  ## checks parameter consistency
  if (verbose)
    cat("-Checking parameter consistency...\n")
  tp1 <- unlist(tni1@para)
  tp2 <- unlist(tni2@para)
  idx <- tp1%in%tp2
  if(!all(idx)){
    tp <- paste("TNIs were not computed using the same parameters!",
                "The following patameters seem to differ between the two:\n",
                paste(names(tp1)[idx], collapse = "\n "))
    warning(tp)
  }
  ## checks whether both TNIs have undergone all methods in RTN from
  ## Permutation to DPI filter
  if (verbose)
    cat("-Checking if all TNI methods are completed...\n")
  if(any(tni1@status[1:4] != "[x]") || any(tni2@status[1:4] != "[x]")){
    ## gives feedback on which methods were not run
    if (verbose){
      cat("TNI1: ")
      print(tni1@status)
      cat("TNI2: ")
      print(tni2@status)
    }
    tp <- paste("NOTE: both TNIs must be evaluated by the RTN pipeline,",
                "up to the DPI filter step!")
    stop(tp, call. = FALSE)
  }
}

##------------------------------------------------------------------------------
#internal function for 'mbrDuals'
##it checks the consistency of supplementaryTable
.checkConsistencySuppTable <- function(object, supplementaryTable, verbose)
  {
  ##---
  evidenceColname <- colnames(supplementaryTable)[3]
  dualsInformation <- mbrGet(object, what="dualsInformation")
  colnms <- colnames(dualsInformation)
  if(evidenceColname%in%colnms)
  {
    cat("-NOTE: evidence table has been already provided, overwriting information...\n")
  }
  ##-----check consistency
  regs <- c(dualsInformation$Regulon1,dualsInformation$Regulon2)
  regs <- unique(regs)
  regsStab <- c(supplementaryTable$Regulon1,supplementaryTable$Regulon2)
  regsStab <- unique(regsStab)
  ##---consistency between supplementaryTable and annotation
  consc <- (sum(regs%in%regsStab)/length(regs))*100
  if(consc<70){
    tp <- paste("Only ",round(consc,1),"% of the regulatory elements ",
                "are represented in the 'supplementaryTable'!\n", sep="")
    warning(tp)
  } else if(consc>=70 && verbose){
    tp <- paste("-",round(consc,1),"% of the regulatory elements ",
                "are represented in the 'supplementaryTable'!\n", sep="")
    cat(tp)
  }
}
##------------------------------------------------------------------------------
#internal function for 'mbrDuals'
##it checks the evidences for 'duals' in 'supplementaryTable'
.updateEvidenceTable <- function (object, supplementaryTable, verbose=TRUE)
{
  evidenceColname <- colnames(supplementaryTable)[3]
  dualsInformation <- mbrGet(object, what="dualsInformation")
  dualsInformation[, evidenceColname] <- NA
  ##---
  if(verbose){
    tp <- paste("-Transferring evidences from 'supplementaryTable'",
                "to inferred duals...\n", sep=" ")
    cat(tp)
  }
  duals <- rownames(dualsInformation)
  ##--- Test A-B ordering, and update if available
  rownames(supplementaryTable) <- paste(supplementaryTable$Regulon1,supplementaryTable$Regulon2, sep = "~")
  suppl <- supplementaryTable[rownames(supplementaryTable)%in%duals,]
  if(nrow(suppl)>0){
    dualsInformation[rownames(suppl),evidenceColname] <- suppl[[evidenceColname]]
  }
  ##--- Test B-A ordering, and update if available
  rownames(supplementaryTable) <- paste(supplementaryTable$Regulon2,supplementaryTable$Regulon1, sep = "~")
  suppl <- supplementaryTable[rownames(supplementaryTable)%in%duals,]
  if(nrow(suppl)>0){
    dualsInformation[rownames(suppl),evidenceColname] <- suppl[[evidenceColname]]
  }
  ##--- update object
  object <- .mbr.set(name="dualsInformation", para=dualsInformation, object=object)
  return (object)
}
##------------------------------------------------------------------------------
#".mbr.set" internal function
##it setts the slots of a MBR object
.mbr.set <- function(name, para, object)
    {
        if(name=="para")
        {
            object@para <- para
        }
        else if(name=="summary")
        {
            object@summary <- para
        }
        else if(name=="status")
        {
            object@status <- para
        }
        else if(name=="TNI1")
        {
            object@TNI1 <- para
        }
        else if(name=="TNI2")
        {
            object@TNI2 <- para
        }
        else if(name=="statusUpdate")
        {
            object@status[para] <- "[x]"
        }
        else if(name=="testedElementsTNI1")
        {
            object@testedElementsTNI1 <- para
        }
        else if(name=="testedElementsTNI2")
        {
            object@testedElementsTNI2 <- para
        }
        else if(name=="dualRegulons")
        {
            object@dualRegulons <- para
        }
        else if(name=="dualsInformation")
        {
            object@results$dualsInformation <- para
        }
        else if(name=="hypergeometricResults")
        {
            object@results$hypergeometricResults <- para
        }
        
        return(object)
    }
