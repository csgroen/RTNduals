
################################################################################
##########################    plot dual regulons    ############################
################################################################################

#' Plot shared target clouds between dual regulons.
#'
#' This function plots the shared target clouds between a regulon pair.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the method \code{\link[RTNduals:mbrAssociation]{mbrAssociation}}.
#' @param names.duals A vector with 'dual regulon' indentifiers from the 
#' 'dualsInformation' table.
#' @param filepath A character string indicating the file path where the plot 
#' should be saved.
#' @param alpha  The alpha transparency, a number in [0,1].
#' @param lncols A vector of length 2 indicating the colors of the negative 
#' and positive target clouds, respectively.
#' @return A plot with the shared target clouds between dual regulons.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrPreprocess
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, rowAnnotation=annot)
#' 
#' ##--- run mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' 
#' ##--- run mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' 
#' ##-- run mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#' 
#' ##--- run mbrAssociation
#' rmbr <- mbrAssociation(rmbr, prob=0.75)
#' 
#' ##--- run mbrDuals
#' rmbr <- mbrDuals(rmbr)
#' 
#' ##--- get inferred duals and plot the shared cloud of targets
#' duals <- mbrGet(rmbr, what="dualRegulons")
#' mbrPlotDuals(rmbr, names.duals=duals[1])
#'
#' @importFrom grDevices adjustcolor dev.off pdf colorRampPalette col2rgb
#' @importFrom graphics abline axis par plot.new plot.window points title 
#' legend
#' @export

##------------------------------------------------------------------------------
mbrPlotDuals <- function(object, names.duals=NULL, filepath=NULL, 
                         alpha=0.80, lncols=c("darkgreen","darkorange3"))
{
  ##----check object class
  mbr.checks(name="object", para=object)
  mbr.checks(name="names.duals", para=names.duals)
  mbr.checks(name="filepath", para=filepath)
  mbr.checks(name="alpha", para=alpha)
  mbr.checks(name="lncols", para=lncols)
  ##---
  rtni <- .merge.tnis(object)
  rtni_para <- tni.get(rtni, what="para")
  estimator <- rtni_para$perm$estimator
  motifstb <- mbrGet(object, what="dualsInformation")
  motifstb <- .namesDual.check(motifstb, names.duals)
  for(i in 1:nrow(motifstb)){
    mtfs <- motifstb[i,]
    mtfs <- as.character(mtfs)
    reg1 <- mtfs[1]
    reg2 <- mtfs[2]
    rval <- as.numeric(mtfs[3])
    pval <- as.numeric(mtfs[4])
    labelMotif <- paste("dual_", reg1,"_vs_" ,reg2, sep="")
    if(!is.null(filepath)){
      filename <- paste(path.expand(filepath), labelMotif, sep="/")
    } else {
      filename <- NULL
    }
    .tni.plot.greement(rtni=rtni, duals=c(reg1, reg2), corVal=rval,
                       pval=pval, filename=filename, lncols=lncols,
                       alpha=alpha, estimator=estimator
                      )
  }
}

##------------------------------------------------------------------------------
##subfunction for 'mbrPlotDuals'
.tni.plot.greement<-function(rtni, duals, corVal, pval, filename=NULL, 
                             lncols=c("blue","red"), bgcols=lncols, 
                             alpha=0.80, estimator='spearman', lwd=0.70,
                             sharedTargets=TRUE, mapAssignedAssociation=TRUE)
{
  ##---
  tfs <- tni.get(rtni, "regulatoryElements")
  idx1 <- match(duals, names(tfs))
  idx2 <- match(duals, tfs)
  idxcheck<-which(is.na(idx1))
  idx1[idxcheck]<-idx2[idxcheck]
  duals<-tfs[idx1]
  ##---
  reftnet <- tni.get(rtni, "refnet")
  gexp <- tni.get(rtni, "gexp")
  tnet<-reftnet[,duals]
  xy<-.tni.cor(gexp,tnet,asInteger=FALSE,estimator=estimator, 
               mapAssignedAssociation=mapAssignedAssociation)
  if(sharedTargets)
  {
    idx<-rowSums(tnet!=0)==2
    tnet<-tnet[idx,]
    xy<-xy[idx,]
  } else {
    idx<-rowSums(xy!=0)>=1
    tnet<-tnet[idx,]
    xy<-xy[idx,]
  }
  ##---
  xlab=paste(names(duals)[1],"targets (R)")
  ylab=paste(names(duals)[2],"targets (R)")
  xlim=c(-1.0,1.0)
  ylim=c(-1.0,1.0)
  bgcols[1]<-colorRampPalette(c(lncols[1],"white"))(30)[15]
  #bgcols[2]<-"white"
  bgcols[2]<-colorRampPalette(c(lncols[2],"white"))(30)[15]
  #bgcols[4]<-"white"
  bgcols<-adjustcolor(bgcols,alpha.f=alpha)
  ##---plot
  if(!is.null(filename)){
    pdf(file=paste(filename,".pdf",sep=""), height=3, width=3)
  }
  par(mgp=c(2.2, 0.5, 0),mar=c(3.5, 3.5, 1, 1) + 0.1)
  plot.new()
  plot.window(ylim=xlim,xlim=ylim)
  axis(2,cex.axis=1,las=1,tcl=-0.15,lwd=2)
  axis(1,cex.axis=1,las=1,tcl=-0.15,lwd=2)
  title(xlab=xlab,ylab=ylab,cex.lab=1)
  
  ##---legend
  if(pval < 2e-16){
    pval <- "< 2e-16"
  } else {
    pval <- paste("= ",signif(pval,2),sep="")
  }
  legs <- paste("Adj.Pval ", pval, sep="")
  
  if(corVal<0){
    ##---negative Dual
    tpp<-xy[(sign(tnet[, 1])==1 & sign(tnet[, 2])==-1),]
    points(tpp, col=lncols[1], pch=21, cex=0.7, bg=bgcols[1], lwd=lwd)
    tpp<-xy[sign(tnet[, 1])==-1 & sign(tnet[, 2])==1,]
    points(tpp,col=lncols[1],pch=21,cex=0.7,bg="white", lwd=lwd)
    legend("topright", legs, bty="n", cex = 0.9)
  } else {
    ##---positive Dual
    tpp<-xy[rowSums(sign(tnet))==2, ]
    points(tpp,col=lncols[2],pch=21,cex=0.7,bg="white", lwd=lwd)
    tpp<-xy[rowSums(sign(tnet))==-2, ]
    points(tpp,col=lncols[2],pch=21,cex=0.7,bg=bgcols[2], lwd=lwd)
    legend("bottomright", legs, bty="n", cex = 0.7)
  }
  
  ##---lines
  abline(h=0, lwd=1.5, lty="12", col="grey70")
  abline(v=0, lwd=1.5, lty="12", col="grey70")
  if(!is.null(filename)){
    dev.off()
    tp <- paste("- file '", filename,".pdf' has been generated!\n", sep="")
    cat(tp)
  }
  ##---report
  colnames(xy)<-paste(names(duals),"(R)",sep="")
  nms<-rownames(xy)
  annot<-rtni@rowAnnotation[nms,]
  report<-cbind(annot,format(round(xy,3)))
  invisible(report)
}

##subfunction for 'mbrPlotDuals'
.merge.tnis <- function (object)
{
  TNI1 <- mbrGet(object, "TNI1")
  TNI2 <- mbrGet(object, "TNI2")
  elreg1 <- tni.get(TNI1, "regulatoryElements")
  elreg2 <- tni.get(TNI2, "regulatoryElements")
  elregs <- c(elreg1, elreg2)
  elregs <- elregs[!duplicated(elregs)]
  rtni_merge <- new("TNI",gexp = tni.get(TNI1, "gexp"), regulatoryElements = elregs)
  rtni_merge@rowAnnotation <- object@TNI1@rowAnnotation
  rtni_merge@para <- tni.get(TNI1, "para")
  #---
  tnet1 <- tni.get(TNI1, "refnet")[, elreg1]
  tnet2 <- tni.get(TNI2, "refnet")[, setdiff(elreg2,elreg1)]
  rtni_merge@results$tn.ref <- cbind(tnet1, tnet2)
  rtni_merge@status [1:3] <- "[x]"
  return (rtni_merge)
}

##subfunction for 'mbrPlotDuals'
.namesDual.check <- function(motifstb, names.duals)
{
  if (!is.null (names.duals))
  {
    ##----checks names.duals
    if(sum(names.duals%in%rownames(motifstb)) == 0) 
      stop("NOTE: 'names.duals' should be in 'dualsInformation!' 
           see 'mbrGet' function. \n", call.=FALSE)
    if(sum(names.duals%in%rownames(motifstb)) != 
       length(names.duals)) 
      stop ("NOTE: Not all names are available for dual regulons! \n", 
            call.=FALSE)
    ##----
    motifstb <- motifstb[names.duals, c("Regulon1","Regulon2", "R", "Hypergeometric.Adjusted.Pvalue")]
  } else {
    motifstb <- motifstb[, c("Regulon1", "Regulon2", "R", "Hypergeometric.Adjusted.Pvalue")]
  }
  motifstb[, 3] <- round(motifstb[, 3], 2)
  return(motifstb)
}
