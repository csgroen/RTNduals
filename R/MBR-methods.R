################################################################################
##########################         MBR-methods      ############################
################################################################################
##------------------------------------------------------------------------------
#########################################
###class builder
#########################################
.MBRmaker <- function(gexp=gexp, regulatoryElements1=regulatoryElements1,
                      regulatoryElements2=regulatoryElements2)
{
  #---creates the 'MBR' object
  object <- newMBR(gexp=gexp,
                   regulatoryElements1=regulatoryElements1,
                   regulatoryElements2=regulatoryElements2
  )
  #---status
  status <- rep('[ ]', 1, 5)
  names(status) <- c('Preprocess', 'Permutation','Bootstrap', 
                     'DPI.filter', 'Association')
  #---parameters
  sum.info.para <- list()
  sum.info.para$TNIs$perm <- NA
  sum.info.para$TNIs$boot <- NA
  sum.info.para$TNIs$dpi <- NA
  sum.info.para$MBR$association <- matrix(NA, 1, 3)
  colnames(sum.info.para$MBR$association) <- c('minRegulonSize','prob',
                                               'estimator')
  rownames(sum.info.para$MBR$association) <- 'Parameter'
  #---summary dualsInformation
  sum.info.summary <- list()
  sum.info.summary$MBR$Duals <- matrix(NA, 1, 2)
  colnames(sum.info.summary$MBR$Duals) <- c('testedDuals','inferredDuals')
  rownames(sum.info.summary$MBR$Duals) <- 'duals'
  
  #---set
  object <- .mbr.set(name="status", para=status, object=object)
  object <- .mbr.set(name="para", para=sum.info.para, object=object)
  object <- .mbr.set(name="summary", para=sum.info.summary, object=object)
  return(object)
}

#----------------------------------------------------------
#it creates the 'MBR' class object
newMBR <- function(gexp, regulatoryElements1, regulatoryElements2)
{
  #---checks
  if(missing(gexp)) 
    stop("NOTE: 'gexp' is missing ", call.=FALSE)
  if(missing(regulatoryElements1)) 
    stop("NOTE: 'regulatoryElements1' is missing", call.=FALSE)
  if(missing(regulatoryElements2)) 
    stop("NOTE: 'regulatoryElements2' is missing", call.=FALSE)
  mbr.checks(name='gexp', gexp)
  mbr.checks(name='regulatoryElements1', regulatoryElements1)
  mbr.checks(name='regulatoryElements2', regulatoryElements2)
  
  #---creating TNIs
  regulonsTNI1 <- new("TNI", gexp=gexp, 
                      transcriptionFactors=regulatoryElements1)
  regulonsTNI2 <- new("TNI", gexp=gexp, 
                      transcriptionFactors=regulatoryElements2)
  
  #---creating MBR-object
  new(Class="MBR", TNI1=regulonsTNI1, TNI2=regulonsTNI2)
}

##------------------------------------------------------------------------------
#' A preprocessing function for objects of class MBR.
#'
#' @param gexp A numerical matrix, typically with mRNA and/or miRNA expression 
#' values.
#' @param regulatoryElements1 A named vector with regulatory elements listed in 
#' 'gexp' rownames.
#' @param regulatoryElements2 A named vector with regulatory elements listed in 
#' 'gexp' rownames.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to 
#' \code{\link[RTN:tni.preprocess]{tni.preprocess}} function.
#' @return A preprocessed 'MBR-class' object.
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
#' regulatoryElements2=tfs2, gexpIDs=annot)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrPreprocess-methods
#' @aliases mbrPreprocess
#' @export

##Regulons pre-processing method
setMethod("mbrPreprocess",
          "matrix",
          function(gexp, regulatoryElements1, regulatoryElements2, 
                   verbose=TRUE,...)
          {
            ##---
            mbr.checks(name="verbose", para=verbose)  
            object <- .MBRmaker(gexp=gexp, 
                                regulatoryElements1=regulatoryElements1,
                                regulatoryElements2=regulatoryElements2)  
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            ##---pre-processing TNIs
            if(verbose) cat("-Preprocessing TNI objects...\n\n")
            TNI1 <- tni.preprocess(TNI1, verbose=verbose,...=...)
            TNI2 <- tni.preprocess(TNI2, verbose=verbose,...=...)
            tfs1 <- tni.get(TNI1, what="tfs")
            tfs2 <- tni.get(TNI2, what="tfs")
            mbr.checks(name="regulatoryElements", para=tfs1)
            mbr.checks(name="regulatoryElements", para=tfs2)
            
            ##---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="Preprocess", 
                               object=object)
            
            return(object)
          }
)
##------------------------------------------------------------------------------
#' Inference of transcriptional networks.
#'
#' This function takes an MBR object and computes two transcriptional networks 
#' inferred 
#' by mutual information (with multiple hypothesis testing corrections).
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed on to the 
#' \code{\link[RTN:tni.permutation]{tni.permutation}} function.
#' @return An \linkS4class{MBR} object with two mutual information matrices, 
#' one in each "TNI" slot.
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
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' 
#' ##--- run mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrPermutation-methods
#' @aliases mbrPermutation
#' @export

## permutation
setMethod("mbrPermutation",
          "MBR",
          function(object, verbose=TRUE, ...)
          {
            ##---checks
            mbr.checks(name="object", para=object)
            mbr.checks(name="verbose", para=verbose)
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            
            ##---permutation TNIs
            if(verbose)
              cat("-Performing permutation analysis for two TNI objects...\n\n")
            TNI1 <- tni.permutation(TNI1, verbose=verbose,...=...)
            TNI2 <- tni.permutation(TNI2, verbose=verbose,...=...)
            #---get
            tni1_summary <- tni.get(TNI1, what="summary")
            tni2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            #---changes
            mbr_para$TNIs$perm <- tni1_summary$para$perm
            mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            
            #---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="Permutation", 
                               object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return(object)
          }
)

#' Inference of consensus transcriptional networks.
#'
#' This function takes an MBR object and computes two consensus transcriptional 
#' networks.
#'
#' @param object A processed objec of class \linkS4class{MBR} evaluated by the 
#' method \code{\link[RTNduals:mbrPermutation]{mbrPermutation}}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the 
#' \code{\link[RTN:tni.bootstrap]{tni.bootstrap}} function.
#' @return An \linkS4class{MBR} object with two consensus mutual information 
#' matrices, one in each "TNI" slot.
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
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' 
#' ##--- run mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' 
#' ##--- run mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrBootstrap-methods
#' @aliases mbrBootstrap
#' @export

##------------------------------------------------------------------------------
## bootstrap method
setMethod("mbrBootstrap",
          "MBR",
          function(object, verbose=TRUE, ...)
          {
            ##---checks
            mbr.checks(name="object", para=object)
            mbr.checks(name="verbose", para=verbose)
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            
            ##---bootstrap TNIs
            if(verbose) cat("-Performing bootstrap analysis for two TNI 
                            objects...\n\n")
            TNI1 <- tni.bootstrap(TNI1, verbose=verbose,...=...)
            TNI2 <- tni.bootstrap(TNI2, verbose=verbose,...=...)
            
            #---get
            tni1_summary <- tni.get(TNI1, what="summary")
            tni2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            #---changes
            mbr_para$TNIs$boot <- tni1_summary$para$boot
            
            mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            
            #---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="Bootstrap", 
                               object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return(object)
          }
)

#' A filter based on the Data Processing Inequality (DPI) algorithm.
#'
#' This function takes an MBR object and computes two transcriptional networks 
#' filtered by the data processing inequality algorithm.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by 
#' the methods
#'  \code{\link[RTNduals:mbrPermutation]{mbrPermutation}} and 
#'  \code{\link[RTNduals:mbrBootstrap]{mbrBootstrap}}.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ... Additional arguments passed to the 
#' \code{\link[RTN:tni.dpi.filter]{tni.dpi.filter}} function.
#' @return An \linkS4class{MBR} object with two DPI-filtered mutual information 
#' matrices, one in each "TNI" slot.
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
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' 
#' ##--- run mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' 
#' ##--- run mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' 
#' ##--- run mbrDpiFilter
#' rmbr <- mbrDpiFilter(rmbr)
#'
#' @import RTN 
#'
#' @import methods
#' @docType methods
#' @rdname mbrDpiFilter-methods
#' @aliases mbrDpiFilter
#' @export

##------------------------------------------------------------------------------
## dpi filter method
setMethod("mbrDpiFilter",
          "MBR",
          function(object, verbose=TRUE, ...)
          {
            ##---checks
            mbr.checks(name="object", para=object)
            mbr.checks(name="verbose", para=verbose)
            ##---get TNIs
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            
            ##---Dpi filter TNIs
            if(verbose) cat("-Applying dpi filter for two TNI objects...\n")
            TNI1 <-tni.dpi.filter(TNI1, verbose=verbose, ...=...)
            TNI2 <-tni.dpi.filter(TNI2, verbose=verbose, ...=...)
            #---get
            tni1_summary <- tni.get(TNI1, what="summary")
            tni2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            #---changes
            mbr_para$TNIs$dpi <- tni1_summary$para$dpi
            
            mbr_summary$TNIs$TNI1 <- tni1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- tni2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            
            #---set
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            object <- .mbr.set(name="statusUpdate", para="DPI.filter", 
                               object=object)
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return(object)
          }
)

#' Motifs analysis and inference of 'dual regulons'.
#'
#' This function takes an MBR object and compares the shared regulon 
#' targets in order to test whether regulon pairs agree on the predicted 
#' downstream effects.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' methods \code{\link[RTNduals:mbrPermutation]{mbrPermutation}}, 
#' \code{\link[RTNduals:mbrBootstrap]{mbrBootstrap}} 
#' and \code{\link[RTNduals:mbrDpiFilter]{mbrDpiFilter}}.
#' @param regulatoryElements1 An optional character vector specifying which 
#' 'TNI1' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param regulatoryElements2 An optional character vector specifying which 
#' 'TNI2' regulatory elements should be evaluated. If 'NULL' all regulatory 
#' elements will be evaluated.
#' @param minRegulonSize A single integer or numeric value specifying the 
#' minimum number of elements in a regulon. Gene sets with fewer than this 
#' number are removed from the analysis.
#' @param prob A numeric value, representing a quantile threshold applyed to 
#' the association metric used to infer 'dual regulons'.
#' @param estimator A character value specifying the estimator used in the 
#' association analysis. One of "spearman" (default), "kendall", or "pearson".
#' @param pAdjustMethod A single character value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' function for details).
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with two data.frames in the slot 
#' 'results' listing the inferred 'dual regulons' and correspoding statistics.
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
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' 
#' ##--- run mbrPermutation
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' 
#' ##--- run mbrBootstrap
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' 
#' ##--- run mbrAssociation
#' rmbr <- mbrAssociation(rmbr, prob=0.75)
#'
#' @import RTN 
#' @importFrom stats p.adjust phyper
#' @importFrom stats cor quantile
#'
#' @import methods
#' @docType methods
#' @rdname mbrAssociation-methods
#' @aliases mbrAssociation
#' @export

##------------------------------------------------------------------------------
##Inference of duals
setMethod("mbrAssociation",
          "MBR",
          function(object, regulatoryElements1=NULL, regulatoryElements2=NULL, 
                   minRegulonSize=30, prob=0.95, estimator='spearman', 
                   pAdjustMethod="BH", verbose=TRUE)
          {
            ##--- gets
            TNI1 <- mbrGet(object, what="TNI1")
            TNI2 <- mbrGet(object, what="TNI2")
            tni_gexp <- tni.get(TNI1, what="gexp")
            tni_para <- tni.get(TNI1, what="para")
            ##--- checks
            mbr.checks(name="minRegulonSize", para=minRegulonSize)
            mbr.checks(name="prob", para=prob)
            mbr.checks(name="estimator", para=estimator)
            mbr.checks(name="pAdjustMethod", para=pAdjustMethod)
            mbr.checks(name="verbose", para=verbose)
            if(is.null(regulatoryElements1))
            {
              if(verbose) 
                cat("-Selecting regulatory elements from TNI1 object...\n")
              regulatoryElements1 <- tni.get(TNI1, "tfs")
            } else
            {
              regulatoryElements1 <- .checkRegel(TNI1, regulatoryElements1)
            }
            if(is.null(regulatoryElements2))
            {
              if(verbose) 
                cat("-Selecting regulatory elements from TNI2 object...\n")
              regulatoryElements2 <- tni.get(TNI2, "tfs")
            } else
            {
              regulatoryElements2 <- .checkRegel(TNI2, regulatoryElements2)
            }
            mbr.checks(name="numberRegElements", para=regulatoryElements1)
            mbr.checks(name="numberRegElements", para=regulatoryElements2)
            
            ##--- get regulons
            what <- "refregulons.and.mode"
            regulons1 <- tni.get(TNI1, what=what)
            regulons2 <- tni.get(TNI2, what=what)
            
            ##--- get regulatory elements
            regulons1 <- regulons1[regulatoryElements1]
            regulons2 <- regulons2[regulatoryElements2]
            
            ##--- get regulons by min size
            size1 <- unlist(lapply(regulons1, length))
            size2 <- unlist(lapply(regulons2, length))
            ##---
            if( sum(size1)==0 | sum(size2)==0){
              stop("NOTE: at least one input regulon should be above the 
                   'minRegulonSize' in both TNIs!",call.=FALSE)
            }
            ##---
            idx <- size1 >= (minRegulonSize)
            regulons1 <- regulons1[idx]
            regulatoryElements1 <- regulatoryElements1[idx]
            size1 <- size1[idx]
            ##---
            idx <- size2 >= (minRegulonSize)
            regulons2 <- regulons2[idx]
            regulatoryElements2 <- regulatoryElements2[idx]
            size2 <- size2[idx]
            
            ##--- group regulons and regulatory elements
            regel <- unique(c(regulatoryElements1, regulatoryElements2))
            targets <- unique(unlist(lapply(c(regulons1, regulons2), names)))
            targets <- unique(c(regel, targets))
            
            ##--- get mi
            if(verbose) 
              cat("-Extrating inferred regulatory associations...\n")
            tbmi1 <- .regMatrix(regulons1, regulatoryElements1, targets, getNames=FALSE)
            tbmi2 <- .regMatrix(regulons2, regulatoryElements2, targets, getNames=FALSE)
            
            ##--- compute correlation between regulators and targets
            if(verbose) 
              cat("-Computing correlation statistics between regulators and targets...\n")
            tnet1 <- .tni.cor(tni_gexp, tbmi1, estimator=estimator, dg=0, 
                              asInteger=FALSE, mapAssignedAssociation=TRUE)
            tnet2 <- .tni.cor(tni_gexp, tbmi2, estimator=estimator, dg=0, 
                              asInteger=FALSE, mapAssignedAssociation=TRUE)
            
            ##--- compute correlation between regulons
            if(verbose) 
              cat("-Computing correlation statistics between regulon pairs...\n")
            regcor <- cor(tnet1, tnet2, method=estimator)
            intregs <- intersect(regulatoryElements1,regulatoryElements2)
            diag(regcor[intregs,intregs]) <- NA
            regcor[intregs,intregs][upper.tri(regcor[intregs,intregs])] <- NA
            rownames(regcor) <- names(regulatoryElements1)
            colnames(regcor) <- names(regulatoryElements2)
            
            ##--- 
            testedDuals <- sum(!is.na(regcor))
            if(testedDuals < 100){
              tp1 <- paste("NOTE: only",testedDuals,"regulon pair(s) is(are) being tested!\n")
              tp2 <- "Ideally, the search space should represent all possible\n"
              tp3 <- "combinations of a given class of regulators! For example,\n"
              tp4 <- "all nuclear receptors annotated for a given species."
              warning(tp1,tp2,tp3,tp4, call.=FALSE)
            }
            ##--- select motifs based on 'prob' quantile
            if(verbose) cat("-Computing quantile statistics...\n")
            statlist <- .motifsquantile(regcor=regcor, prob=prob)
            inferredDuals <- nrow(statlist)
            
            if(inferredDuals > 0){
              ##--- Mutual Information
              if(verbose)cat("-Computing Mutual Information...\n")
              statlist <- .getMI(statlist, tbmi1, regulatoryElements1, 
                                 regulatoryElements2, size1, size2)
              ##--- PadjustValue
              if(!is.null(TNI1@results$adjpv)) {
                pvmat <- TNI1@results$adjpv
                cutoff <- tni_para$perm$pValueCutoff
                statlist <- .getPMI(statlist, pvmat, regulatoryElements1, 
                                    regulatoryElements2, cutoff=cutoff)
              }
              ##--- Jaccard
              if(verbose)cat("-Computing Jaccard similarity...\n")
              statlist <- .jcOverlap(statlist, regulatoryElements1, 
                                     regulatoryElements2, 
                                     tnet1, tnet2)
              ##--- Hypergeometric
              if(verbose)cat("-Running hypergeometric analysis...\n")
              universe <- rownames(tni_gexp)
              hyperresults <- .mbr.hyper(.statlist=statlist, 
                                         regulons1=regulons1, regulons2=regulons2, 
                                         regulatoryElements1, regulatoryElements2,
                                         universe=universe, pAdjustMethod=pAdjustMethod, 
                                         verbose=verbose)
              statlist$Hypergeometric.Pvalue <- hyperresults$Pvalue
              statlist$Hypergeometric.Adjusted.Pvalue <- hyperresults$Adjusted.Pvalue
            } else {
              warning("No 'dual regulon' has been observed for the input parameters.",
                      call.=FALSE)
            }
            ##--- organize statlist
            if("MI.Adjusted.Pvalue"%in%colnames(statlist)){
              statlist <- statlist[,c("Regulon1","Size.Regulon1","Regulon2",
                                  "Size.Regulon2","Jaccard.coefficient", 
                                  "Hypergeometric.Pvalue",
                                  "Hypergeometric.Adjusted.Pvalue", 
                                  "MI","MI.Adjusted.Pvalue","R","Quantile")]
            } else {
              statlist <- statlist[,c("Regulon1","Size.Regulon1","Regulon2",
                                  "Size.Regulon2","Jaccard.coefficient", 
                                  "Hypergeometric.Pvalue",
                                  "Hypergeometric.Adjusted.Pvalue", 
                                  "MI","R","Quantile")]
            }
            
            ##--- para
            mbr_para <- mbrGet(object,what="para")
            sum.info.par <- c(minRegulonSize, prob, estimator)
            mbr_para$MBR$association['Parameter', ] <- sum.info.par
            
            ##---
            mbr_summary <- mbrGet(object, what="summary")
            mbr_summary$MBR$Duals[,'testedDuals'] <- testedDuals
            mbr_summary$MBR$Duals[,'inferredDuals'] <- inferredDuals
            ##--- set
            object <- .mbr.set(name="statusUpdate", 
                               para="Association", object=object)
            object <- .mbr.set(name="para", 
                               para=mbr_para, object=object)
            object <- .mbr.set(name="summary", 
                               para=mbr_summary, object=object)
            object <- .mbr.set(name="testedElementsTNI1", 
                               para=regulatoryElements1, object=object)
            object <- .mbr.set(name="testedElementsTNI2", 
                               para=regulatoryElements2, object=object)
            object <- .mbr.set(name="dualRegulons", 
                               para=rownames(statlist), object=object)
            object <- .mbr.set(name="dualsInformation", 
                               para=statlist, object=object)
            object <- .mbr.set(name="hypergeometricResults", 
                               para=hyperresults, object=object)    
            return(object)
            }
)


#' A summary for results from the MBR methods.
#'
#' This function lists the inferred 'dual regulons' and, if available, 
#' adds external evidences.
#'
#' @param object A processed object of class \linkS4class{MBR} evaluated by the 
#' method \code{\link[RTNduals:mbrAssociation]{mbrAssociation}}.
#' @param supplementary.table An optional 'data.frame' with three columns 
#' representing 
#' (1) regulatory elements of 'TNI1', (2) regulatory elements of 'TNI2', and 
#' (3) external evidences between the regulatory elements.
#' @param evidenceColname A single character value specifying a column in 
#' the 'supplementary.table'.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when verbose=TRUE) or not (when verbose=FALSE).
#' @return An \linkS4class{MBR} object with an updated 'data.frame' in the slot 
#' 'results' listing the input additional evidences.
#' @examples
#' ##--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' ##--- run mbrAssociation
#' rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1 = tfs1, 
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' rmbr <- mbrPermutation(rmbr, nPermutations=10)
#' rmbr <- mbrBootstrap(rmbr, nBootstrap=10)
#' rmbr <- mbrAssociation(rmbr, prob=0.75)
#' rmbr <- mbrDuals(rmbr)
#' 
#' ##--- check results
#' results <- mbrGet(rmbr, what="dualsInformation")
#' 
#' ##--- add supplementary evidences
#' ## here we build a 'toy' example using the 'rnorm' function for 
#' ## demonstration purposes only!
#' supplementaryTable <- results[,c("Regulon1","Regulon2")]
#' supplementaryTable$ToyEvidence <- rnorm(nrow(results))
#' supplementaryTable
#' 
#' ##--- add supplementary evidences with brDuals function
#' rmbr <- mbrDuals(rmbr, supplementary.table=supplementaryTable,
#' evidenceColname = "ToyEvidence")
#' 
#' ##--- check updated results
#' mbrGet(rmbr, what="dualsInformation")
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import methods
#' @docType methods
#' @rdname mbrDuals-methods
#' @aliases mbrDuals
#' @export

##------------------------------------------------------------------------------
##organize duals
setMethod( "mbrDuals",
           "MBR",
           function(object, supplementary.table=NULL, evidenceColname, 
                    verbose=TRUE)
           {
             ##---checks
             mbr.checks(name="object", para=object)
             ##---
             dualsInformation <- mbrGet(object, what="dualsInformation")
             if(is.null(dim(dualsInformation)))
               stop("empty results in the input 'object'!",call.=FALSE)
             if(verbose) cat("-Sorting by the R value...\n")
             idx <- sort(abs(dualsInformation[,"R"]), decreasing=TRUE, 
                         index.return=TRUE)
             dualsInformation <- dualsInformation[idx$ix, ]
             
             #---set
             object <- .mbr.set(name="dualsInformation", 
                                para=dualsInformation, object=object)
             object <- .mbr.set(name="dualRegulons", 
                                para=rownames(dualsInformation), object=object)
             if(!is.null(supplementary.table))
             {
               ##---checks
               if(missing(evidenceColname)) 
                 stop("'evidenceColname' should be a character value present in 
                      colnames of supplementary.table!",call.=FALSE)
               mbr.checks(name="supplementary.table", para=supplementary.table)
               mbr.checks(name="uniqueInput", para=supplementary.table)
               mbr.checks(name="evidenceColname", para=evidenceColname)
               
               ##---consistency
               if(verbose) 
                 cat("-Checking the 'supplementary.table' consistency...\n")
               supplementary.table <- .consisSuppTable(object, 
                                                       supplementary.table, 
                                                       evidenceColname, 
                                                       verbose=verbose)
               ##---find duals
               object <- .checkLoops(object, supplementary.table, 
                                     evidenceColname, verbose=verbose)
             }
             return(object)
           }
)

#' A preprocessing function for objects of class MBR.
#'
#' This function merges two TNI class objects and creates one MBR class object.
#'
#' @param TNI1 A 'TNI' class object.
#' @param TNI2 Another 'TNI' class object
#' @param verbose A single logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' @param regulatoryElements1 A character vector specifying which 
#' 'TNI1' regulatory elements should be evaluated.
#' @param regulatoryElements2 A character vector specifying which 
#' 'TNI2' regulatory elements should be evaluated.
#' @return An \linkS4class{MBR} object.
#' @examples
#' #--- load a dataset for demonstration
#' data("dt4rtn", package = "RTN")
#' gexp <- dt4rtn$gexp
#' annot <- dt4rtn$gexpIDs
#' tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#' tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]
#' 
#' \dontrun{
#' 
#' ##--- compute a TNI for tfs1
#' tni1 <- new("TNI", gexp=gexp, transcriptionFactors=tfs1)
#' tni1 <- tni.preprocess(tni1, gexpIDs=annot)
#' tni1 <-tni.permutation(tni1)
#' tni1 <-tni.bootstrap(tni1)
#' 
#' ##--- compute a TNI for tfs2
#' tni2 <- new("TNI", gexp=gexp, transcriptionFactors=tfs2)
#' tni2 <- tni.preprocess(tni2, gexpIDs=annot)
#' tni2 <-tni.permutation(tni2)
#' tni2 <-tni.bootstrap(tni2)
#' 
#' ##--- run tni2mbrPreprocess
#' rmbr <- tni2mbrPreprocess(tni1, tni2)
#' }
#'
#' @import methods
#' @docType methods
#' @rdname tni2mbrPreprocess-methods
#' @aliases tni2mbrPreprocess
#' @export

##------------------------------------------------------------------------------
##Combine two TNIs produced separately
setMethod("tni2mbrPreprocess",
          "TNI",
          function (TNI1,  TNI2, 
                    regulatoryElements1=NULL, regulatoryElements2=NULL, 
                    verbose=TRUE)
          {
            if(missing(TNI1)) stop("NOTE: 'tni1' is missing ", call.=FALSE)
            if(missing(TNI2)){
              TNI2 <- TNI1 ##stop("NOTE: 'TNI2' is missing ", call.=FALSE)
              mbr.checks(name="regulatoryElements1", para=regulatoryElements1)
              mbr.checks(name="regulatoryElements2", para=regulatoryElements2)
              ##----
              regulatoryElements1 <- .checkRegel(TNI1, regulatoryElements1)
              regulatoryElements2 <- .checkRegel(TNI2, regulatoryElements2)
              ##----remove duplicated Regulatory Elements
              TNI1@transcriptionFactors <- regulatoryElements1
              TNI2@transcriptionFactors <- regulatoryElements2
            }
            
            mbr.checks (name='tni', para=TNI1)
            mbr.checks (name='tni', para=TNI2)
            .combineTNIs (tni1=TNI1, tni2=TNI2, verbose=verbose)
            #---get
            gexp <- tni.get(TNI1, what="gexp")
            regulatoryElements1 <- tni.get(TNI1, what="tfs")
            regulatoryElements2 <- tni.get(TNI2, what="tfs")
            ##---- creates MBR object
            object <- .MBRmaker(gexp=gexp,
                                regulatoryElements1=regulatoryElements1,
                                regulatoryElements2=regulatoryElements2)
            #---TNIs Update
            object <- .mbr.set(name="TNI1", para=TNI1, object=object)
            object <- .mbr.set(name="TNI2", para=TNI2, object=object)
            #---statu update
            TNI_status <- tni.get(TNI1, what="status")
            status <- names(TNI_status[TNI_status=="[x]"])
            object <- .mbr.set(name="statusUpdate", para=status, object=object)
            #---get Updates
            TNI1_summary <- tni.get(TNI1, what="summary")
            TNI2_summary <- tni.get(TNI2, what="summary")
            mbr_summary <- mbrGet(object, what="summary")
            mbr_para <- mbrGet(object, what="para")
            
            ##---permutation
            
            mbr_para$TNIs$perm <- TNI1_summary$para$perm
            mbr_summary$TNIs$TNI1 <- TNI1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- TNI2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            ##---bootstrap
            mbr_para$TNIs$boot <- TNI1_summary$para$boot
            mbr_summary$TNIs$TNI1 <- TNI1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- TNI2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            ##---summary dpi.filter
            mbr_para$TNIs$dpi <- TNI1_summary$para$dpi
            mbr_summary$TNIs$TNI1 <- TNI1_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI1) <- c('RE', 'Targets', 'Edges')
            mbr_summary$TNIs$TNI2 <- TNI2_summary$results$tnet
            colnames(mbr_summary$TNIs$TNI2) <- c('RE', 'Targets', 'Edges')
            #---set
            object <- .mbr.set(name="para", para=mbr_para, object=object)
            object <- .mbr.set(name="summary", para=mbr_summary, object=object)
            return (object)
          }
)

##------------------------------------------------------------------------------
##show summary information on screen
setMethod( "show",
           "MBR",
           function(object)
           {
             cat("an MBR (Motifs Between Regulons) object:\n")
             message("--status:")
             print(object@status, quote=FALSE)
           }
)

#' Get information from individual slots in MBR object.
#' 
#' Get information from individual slots in an MBR object and any available 
#' results from previous analysis.
#' 
#' @param object A preprocessed object of class \linkS4class{MBR}
#' @param what a single character value specifying which information should be 
#' retrieved from the slots. Options: "TNI1", "TNI2", "testedElementsTNI1", 
#' "testedElementsTNI2", "dualRegulons", "results", "para", "summary", 
#' "status", "dualsInformation" and "hyperResults"
#' @return A slot content from a object of class 'MBR' \linkS4class{MBR} object
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
#' regulatoryElements2=tfs2, gexpIDs=annot)
#' 
#' ##--- get the 'TNI1' slot using 'mbrGet'
#' tni1 <- mbrGet(rmbr, what="TNI1")
#' 
#' @import methods
#' @docType methods
#' @rdname mbrGet-methods
#' @aliases mbrGet
#' @export
##------------------------------------------------------------------------------
##get slots from MBR object
setMethod( "mbrGet",
           "MBR", 
           function(object, what="status")
           {
             ##---check input arguments
             mbr.checks(name="object", para=object)
             mbr.checks(name="mbrGet", para=what)
             ##---Association options any change needs update!
             optsAssoci <- c("testedElementsTNI1", "testedElementsTNI2", 
                             "dualRegulons", "dualsInformation", "results", 
                             "hyperResults")
             ##---get query
             if(what=="TNI1")
             {
               query <- object@TNI1
             }
             else if(what=="TNI2")
             {
               query <- object@TNI2
             }
             else if(what=="para")
             {
               query <- object@para
             }
             else if(what=="summary")
             {
               query <- object@summary
             }
             else if(what=="status")
             {
               query <- object@status
             }
             else if(what%in%optsAssoci)
             {
               if(object@status["Association"] != "[x]")
               {
                 warning("NOTE: input 'object' needs 'mbrAssociation' evaluation!",
                         call.=FALSE)
                 query <- NULL
               } else {
                 if(what=="testedElementsTNI1")
                 {
                   query <- object@testedElementsTNI1
                 }
                 else if(what=="testedElementsTNI2")
                 {
                   query <- object@testedElementsTNI2
                 }
                 else if(what=="dualRegulons")
                 {
                   query <- object@dualRegulons
                 }
                 else if(what=="results")
                 {
                   query <- object@results
                 }
                 else if(what=="dualsInformation")
                 {
                   query <- object@results$dualsInformation
                 }
                 else if(what=="hyperResults")
                 {
                   query <- object@results$hypergeometricResults
                 }
               }
             }
             return(query)
           }
)
