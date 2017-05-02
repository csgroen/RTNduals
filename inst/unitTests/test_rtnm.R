# Unit tests for MBR-class methods
test_mbr <- function()
{
    data("dt4rtn", package = "RTN")
    tfs1 <- dt4rtn$tfs[c("FOXM1", "E2F2")]
    tfs2 <- dt4rtn$tfs[c("PTTG1", "RARA")]
    ##mbrPreprocess
    rmbr <- mbrPreprocess(gexp=dt4rtn$gexp, regulatoryElements1=tfs1, 
                           regulatoryElements2=tfs2, gexpIDs=dt4rtn$gexpIDs)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[1]=="[x]")
    ##mbrPermutation
    rmbr <- mbrPermutation(rmbr, nPermutations=10, estimator="pearson")
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[2]=="[x]")
    ##mbrBootstrap
    rmbr <- mbrBootstrap(rmbr, estimator="pearson", nBootstrap=10)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[3]=="[x]")
    ##mbrDpiFilter
    rmbr <- mbrDpiFilter(rmbr)
    status <- mbrGet(rmbr, what="status")
    checkTrue(status[4]=="[x]")
    ##mbr.combine.TNIs
    tni1 <- mbrGet(rmbr, what="TNI1")
    tni2 <- mbrGet(rmbr, what="TNI2")
    rmbr <- tni2mbrPreprocess(tni1, tni2)
    status <- mbrGet(rmbr, what="status")
    checkTrue(all(status[1:4]=="[x]"))
    ##mbrAssociation
    rmbr <- mbrAssociation(rmbr, prob=0, estimator="pearson")
    status <- mbrGet(rmbr, what="status")
    dualsInformation <- mbrGet(rmbr, what="dualsInformation")
    checkTrue(status[5]=="[x]" && is.data.frame(dualsInformation) )
    ##mbr.motifs
    rmbr <- mbrDuals(rmbr)
    dualsInformation <- mbrGet(rmbr, what="dualsInformation")
    checkTrue(is.data.frame(dualsInformation))
}
