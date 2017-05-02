## ----include=FALSE-------------------------------------------------------
library(RTNduals)
data("dt4rtn", package = "RTN")
gexp <- dt4rtn$gexp
annot <- dt4rtn$gexpIDs
tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]

## ----eval=FALSE----------------------------------------------------------
#  ##--- load package and prepare a dataset for demonstration
#  library(RTNduals)
#  data("dt4rtn", package = "RTN")
#  gexp <- dt4rtn$gexp
#  annot <- dt4rtn$gexpIDs
#  tfs1 <- dt4rtn$tfs[c("IRF8","IRF1","PRDM1","AFF3","E2F3")]
#  tfs2 <- dt4rtn$tfs[c("HCLS1","STAT4","STAT1","LMO4","ZNF552")]

## ----eval=TRUE-----------------------------------------------------------
##--- generate a pre-processed BR-class object
rmbr <- mbrPreprocess(gexp=gexp, regulatoryElements1=tfs1, regulatoryElements2=tfs2, gexpIDs=annot, verbose=FALSE)
rmbr

## ----eval=TRUE-----------------------------------------------------------
##--- compute two regulatory networks
##--- (set nPermutations>=1000)
rmbr <- mbrPermutation(rmbr, nPermutations=100, verbose=FALSE)
rmbr

## ----eval=TRUE-----------------------------------------------------------
##--- check the stability of the two regulatory networks
##--- (set nBootstrap>=100)
rmbr <- mbrBootstrap(rmbr, nBootstrap=10, verbose=FALSE)
rmbr

## ----eval=TRUE-----------------------------------------------------------
##---apply DPI algorithm
rmbr <- mbrDpiFilter(rmbr, eps=0.05, verbose=FALSE)
rmbr

## ----eval=TRUE-----------------------------------------------------------
##--- run the main RTNduals methods
rmbr <- mbrAssociation(rmbr, prob=0.75, verbose=FALSE)
rmbr

## ----eval=TRUE-----------------------------------------------------------
##--- check summary
mbrGet(rmbr, what="summary")

## ----eval=TRUE-----------------------------------------------------------
##--- run 'mbrDuals' and get results
rmbr <- mbrDuals(rmbr)
results <- mbrGet(rmbr, what="dualsInformation")

## ----eval=TRUE-----------------------------------------------------------
##--- here we build a 'toy' evidence table using the 'rnorm' function
supplementaryTable <- results[ ,c("Regulon1","Regulon2")]
supplementaryTable$ToyEvidence <- rnorm(nrow(results))

##--- add supplementary evidences with the 'mbrDuals' function
rmbr <- mbrDuals(rmbr, supplementary.table = supplementaryTable, 
                  evidenceColname = "ToyEvidence", verbose = FALSE)

##--- check updated results
mbrGet(rmbr, what="dualsInformation")

## ----eval=FALSE----------------------------------------------------------
#  duals <- mbrGet(rmbr, what="dualRegulons")
#  mbrPlotDuals(rmbr, names.duals = duals[1])

## ------------------------------------------------------------------------
sessionInfo()

