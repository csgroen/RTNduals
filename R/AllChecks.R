
##------------------------------------------------------------------------------
##This function is used for argument checking
mbr.checks <- function(name, para, paraSuppl)
{
  if(name == "gexp") {
    if(!is.matrix(para) || !is.numeric(para)){
      tp <- paste("NOTE: 'gexp' should be a numeric matrix with, genes on rows",
                  "and samples on cols!")
      stop(tp, call.=FALSE)
    }
    if(is.null(rownames(para)) || is.null(colnames(para)) || 
       length(unique(rownames(para))) < length(rownames(para)) || 
       length(unique(colnames(para))) < length(colnames(para))){
      tp <- paste("NOTE: 'gexp' matrix should be named on rows and cols",
                  "(with unique names)")
      stop(tp, call.=FALSE)
    }
  }
  
  ##---
  else if(name == "regulatoryElements1" || name == "regulatoryElements2"){
    if( !is.character(para) || any(is.na(para)) || 
        any(para == "") || is.null(para)  ){
      tp <- paste("NOTE: 'regulatoryElements' should be a character vector,",
                  "without 'NA' or empty names!")
      stop(tp, call.=FALSE)
    }
    if( any(duplicated(para)) )
      stop("NOTE: 'regulatoryElements' should have unique identifiers!", 
           call.=FALSE)
  }
  
  ##---
  else if(name=="object"){
    if(class(para) != 'MBR')
      stop("NOTE: 'object' should be a 'MBR' class object", call.=FALSE)
  }
  
  ##---
  else if(name=="verbose"){
    if(!is.singleLogical(para))
      stop("NOTE: 'verbose' should be a logical value!", call.=FALSE)
  }
  
  ##---
  else if(name == "minRegulonSize"){
    if(!is.singleNumber(para) || !para>0)
      stop("NOTE: 'minRegulonSize' should be numeric value > 0", call.=FALSE)
  }
  
  ##---
  else if(name == "prob"){
    if(!is.singleNumber(para) || (!para>=0) && (!para<=1))
      stop("NOTE: 'prob' should be a numeric value >= 0 and <= 1!", call.=FALSE)
  }
  
  ##---
  else if(name == "alpha"){
    if(!is.singleNumber(para) || (!para>=0) && (!para<=1))
      stop("NOTE: 'alpha' should be a numeric value >= 0 and <= 1!", call.=FALSE)
  }
  
  ##---
  else if(name == "estimator"){
    if(!is.singleString(para) || !para %in% c("spearman", "kendall", "pearson"))
      stop("NOTE: 'estimator' should be one of 'spearman', 'kendall', 'pearson'!", 
           call.=FALSE)
  }
  
  ##---
  else if(name == "filepath") {
    if(!is.null(para)){
      if(!is.singleString(para) || !dir.exists(para))
        stop("NOTE: 'filepath' should be a valid single path name!", call.=FALSE)
    }
  }
  
  ##---
  else if(name == "names.duals") {
    if(!all.characterValues(para))
      stop("NOTE: 'names.duals' should of character vector!", 
           call.=FALSE)
  }

  ##---
  else if(name == "lncols") {
    if(!is.color(para) || length(para)!=2)
      stop("NOTE: 'lncols' should of a vector (length = 2) with valid colors!", 
           call.=FALSE)
  }
  
  ##---
  else if(name == "regulatoryElements"){
    if(!all.characterValues(para) || any(duplicated(para)) ){
      stop("NOTE: 'regulatoryElements' should be unique character values !", call. = FALSE)
    }
  }
  
  ##---
  else if(name=="numberRegElements"){
    if(length(para)==0)
      stop("NOTE: at least 1 regulatory element should be listed in both 
           'regulatoryElement1' and 'regulatoryElement2'!", call.=FALSE)
  }
  
  ##---
  else if(name=="pAdjustMethod"){
    tp <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    if(!is.character(para) || length(para)!=1 || !(para %in% tp))
      stop("NOTE: 'pAdjustMethod' should be any one of: ",
           paste(tp, collapse = ", "),call.=FALSE)
  }
  
  ##---
  else if(name == "supplementaryTable"){
    #--- general checks
    if(!is.data.frame(para) || !ncol(para)==3 || !nrow(para)>=1 || 
       is.null(dim(para))){
      tp <- paste("NOTE: 'supplementaryTable' should be a 'data.frame' object",
                  "with 3 columns (Regulon1, Regulon2, Evidence)!")
      stop(tp, call.=FALSE)
    }
    if(is.null(colnames(para))){
      stop("NOTE: columns in 'supplementaryTable' should be named!", call.=FALSE)
    }
    bl <- paraSuppl %in% colnames(para)
    if(!is.character(paraSuppl) || is.null(paraSuppl) || 
       !is.singleString(paraSuppl) || !bl){
      tp <- paste("NOTE: 'evidenceColname' should be a character value listed in",
                  "the 'supplementaryTable' colnames!")
      stop(tp, call.=FALSE)
    }
    if(!is.numeric(para[[paraSuppl]])){
      tp <- paste("NOTE: evidences in the 'supplementaryTable' should be",
                  "of numeric or integer type!")
      stop(tp, call.=FALSE)
    }
    #--- set ordering
    col_3 <- which(colnames(para)%in%paraSuppl)
    col_1_2 <- !c(1,2,3) %in% col_3
    col123 <- c(c(1,2,3)[col_1_2],col_3)
    para <- para[,col123]
    colnames(para) <- c("Regulon1","Regulon2",paraSuppl)
    para$Regulon1 <- as.character(para$Regulon1)
    para$Regulon2 <- as.character(para$Regulon2)
    #--- check duplications
    duplnms <- c(paste(para[,1],para[,2],sep="~"), paste(para[,2], para[,1],sep="~"))
    if(sum(duplicated(duplnms))>0) {
      duplnms <- duplnms[duplicated(duplnms)]
      duplnms <- unique(duplnms)
      tp <- paste("NOTE: all pairs of regulators should be unique in 'supplementaryTable'.\n", 
                  "The following pairs seem duplicated: ")
      stop(c(paste(tp, paste(duplnms,collapse=", ")), call.=FALSE))
    }
    return(para)
  }
  
  ##---
  else if(name == "tni"){
    if((class(para) != 'TNI') || is.null(para)){
      stop("NOTE: 'tni1' and 'tni2' should be TNI-class objects!", call.=FALSE)
    }
  }
  
  ##---
  else if(name=="mbrGet"){
    opts <- c("summary", "status", "results", 
              "dualRegulons", "dualsInformation",
              "TNI1", "TNI2", "testedElementsTNI1", 
              "testedElementsTNI2","para"
              )
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("NOTE: 'what' should be any one of the options:", 
                 paste(opts,collapse = ", ") ), call.=FALSE)
  }
  
  }

##------------------------------------------------------------------------------
.txtcollapse<-function(vec){
  paste("'",paste(vec[-length(vec)], collapse = "', '"),
        "'"," and '",vec[length(vec)],"'!", sep=""
  )
}

##------------------------------------------------------------------------------
is.singleNumber <- function(para){
  (is.integer(para) || is.numeric(para)) && length(para)==1L && !is.na(para)
}
is.singleInteger <- function(para){
  lg <- (is.integer(para) || is.numeric(para)) && length(para)==1L && 
    !is.na(para)
  if(lg) lg <- (para / ceiling(para)) == 1
  return(lg)
}
is.singleString <- function(para){
  is.character(para) && length(para) == 1L && !is.na(para)
}
is.singleLogical <- function(para){
  is.logical(para) && length(para) == 1L && !is.na(para)
}
all.binaryValues <- function(para){
  all( para %in% c(0, 1, NA) )
}
all.integerValues <- function(para){
  lg <- ( all(is.integer(para)) || all(is.numeric(para)) ) && 
    !any(is.na(para))
  if(lg) lg <- all ( (para / ceiling(para)) == 1 )
  return(lg)
}
all.characterValues <- function(para){
  all(is.character(para)) && !any(is.na(para))
}
is.color <- function(x){
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}

