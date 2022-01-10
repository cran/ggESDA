#' @name classic2sym
#' @title Convert classical data frame into a symbolic data.
#' @description  A function for converting a classical data, which may
#' present as a data frame or a matrix with one entry one value, into
#' a symbolic data,which is shown as a interval or a set in an entry.
#' Object after converting is ggESDA class containing interval
#' data and raw data(if it exist) and typically statistics.
#' @import rlang stats
#' @importFrom RSDA classic.to.sym
#' @importFrom stringr str_split
#' @ipmortForm tibble rownames_to_column
#' @ipmortForm tibble column_to_rownames
#' @ipmortForm prodlim row.match
#' @param data A classical data frame that you want to be converted into
#' a interval data
#' @param groupby A way to aggregate. It can be either a clustering method
#' or a variable name which exist in input data (necessary factor type)
#' . Default "kmeans".
#' @param k A number of group,which is used by clustering. Default k = 5.
#' @param minData if choose groupby parameter as 'customize',user need to define
#' which data is min data or max data.
#' @param maxData if choose groupby parameter as 'customize',user need to define
#' which data is min data or max data.
#' @usage classic2sym(data=NULL,groupby = "kmeans",k=5,minData=NULL,maxData=NULL)
#' @return classic2sym returns an object of class "ggESDA",which
#' have a interval data and others as follows.
#' \itemize{
#'   \item intervalData - The Interval data after converting also known
#'   as a RSDA object.
#'   \item rawData - Classical data that user input.
#'   \item clusterResult - Cluster results .If the groupby method is
#'   a clustering method then it will exist.
#'   \item statisticsDF - A list contains data frame including some
#'   typically statistics in each group.
#' }
#' @examples
#' #classical data to symbolic data
#' classic2sym(iris)
#' classic2sym(mtcars,groupby = "kmeans",k=10)
#' classic2sym(iris, groupby = "hclust",k=7)
#' classic2sym(iris,groupby=Species)
#'
#' x1<-runif(10,-30,-10)
#' y1<-runif(10,-10,30)
#' x2<-runif(10,-5,5)
#' y2<-runif(10,10,50)
#' x3<-runif(10,-50,30)
#' y3<-runif(10,31,60)
#'
#' d<-data.frame(min1=x1,max1=y1,min2=x2,max2=y2,min3=x3,max3=y3)
#' classic2sym(d,groupby="customize",minData=d[,c(1,3,5)],maxData=d[,c(2,4,6)])
#' classic2sym(d,groupby="customize",minData=d$min1,maxData=d$min2)
#'
#' #extract the data
#' symObj<-classic2sym(iris)
#' symObj$intervalData       #interval data
#' symObj$rawData            #raw data
#' symObj$clusterResult      #cluster result
#' symObj$statisticsDF       #statistics

#' @export
classic2sym<-function(data=NULL,groupby = "kmeans",k=5,minData=NULL,maxData=NULL){
  if("ggESDA" %in% class(data) || "symbolic_tbl" %in% class(data)){
    stop("ERROR : please using classical data frame.")
  }
  pkg.env$rawData <- data
  data <- as.data.frame(data)
  groupby<-substitute(groupby)
  numericData <- unlist(lapply(data.frame(data[1:dim(data)[2]]) ,FUN = is.numeric))
  #start debug version0825-1
  for(i in c(1:dim(data)[2])[!numericData]){
    data[[i]]<-as.factor(data[[i]])
  }
  #end debug
  numericData <- data.frame(data[,which(numericData)])
  if(isNaExist(data)[[1]]){
    warning(paste(isNaExist(data)[[2]],"NA observations will be omit."))
    numericData <- na.omit(numericData)
    data <- na.omit(data)
  }
  switch(toString(groupby),
         kmeans = {
           tryCatch({
             pkg.env$result<-stats::kmeans(numericData, centers=k)
             cluster <- as.factor(pkg.env$result$cluster)
             d<-as.data.frame(cbind(data,cluster=as.factor(pkg.env$result$cluster)))

             pkg.env$statisticsDF<-lapply(pkg.env$statistics,FUN=function(x) aggrByGroup(df=d,group = cluster,method=x))
             names(pkg.env$statisticsDF) <- pkg.env$statistics
             pkg.env$intervalData<-RSDA::classic.to.sym(d,concept = cluster)
           },warning = function(war) {
             print(paste("WARNING in kmeans:  ",war))
           },error = function(err) {
             print(paste("ERROR in kmeans:  ",err))
             print("Automatically using RSDA::classic.to.sym to improve.")
             finallySol(data,numericData)
           })
         },
         #case2 hirechical
         hclust = {
           tryCatch({
             pkg.env$result<-stats::hclust(stats::dist(numericData))
             cluster <- stats::cutree(pkg.env$result, k=k)
             cluster <- as.factor(as.integer(cluster))
             d<-as.data.frame(cbind(data,cluster=cluster))
             pkg.env$statisticsDF<-lapply(pkg.env$statistics,FUN=function(x) aggrByGroup(df=d,group = cluster,method=x))
             names(pkg.env$statisticsDF) <- pkg.env$statistics
             pkg.env$intervalData<-RSDA::classic.to.sym(d,concept = cluster)
           },warning = function(war) {
             print(paste("WARNING in hclust:  ",war))
           },error = function(err) {
             print(paste("ERROR in hclust:  ",err))
             print("Automatically using RSDA::classic.to.sym to improve.")
             finallySol(data,numericData)
           })
         },
         customize={
           idata<-customize(numericData,minData,maxData)
           pkg.env$intervalData<-idata
           pkg.env$statisticsDF <- buildStatsDf(numericData = idata)
           names(pkg.env$statisticsDF) <- c("min","median","max")

         },
         #default
         {
           if(length(as.character(groupby))>1){
             groupby <- as.character(groupby)[-1]
           }else{groupby <- as.character(groupby)[1]}
           d<-RSDA::classic.to.sym(data,groupby)
           d<-buildRowName(data,d,groupby)
           #d<-tibble::as.tibble(d)
           numericData <- unlist(lapply(data.frame(tibble::as.tibble(d)[,1:dim(d)[2]]) ,FUN = RSDA::is.sym.interval))
           numericData <- d[,numericData]

           pkg.env$statisticsDF <- buildStatsDf(numericData = numericData)
           names(pkg.env$statisticsDF) <- c("min","median","max")
           pkg.env$intervalData<-d
         }
  )
  #final
  if(!("symbolic_tbl" %in% class(data))){
    class(pkg.env$intervalData) <- append(class(pkg.env$intervalData),"symbolic_tbl")
  }
  symObj<-ggESDA$new(rawData=pkg.env$rawData,
                       statisticsDF=pkg.env$statisticsDF,
                       intervalData=pkg.env$intervalData,
                       clusterResult = pkg.env$result)

  return(symObj)
}



aggrByGroup <- function(df=NULL,group=NULL,method="mean"){
  group<-rlang::get_expr(group)#trans string to variable
  #1. select numeric data
  numericData<-unlist(lapply(df[,1:dim(df)[2]] ,FUN = is.numeric))
  nData <- as.data.frame(cbind(df[,numericData],group=group))
  #2. aggrate by method
  p<-dim(nData)[2]
  aggData <- data.frame(NULL)
  for (g in levels(group)){
    temp <- round(t(unlist(lapply(as.data.frame(nData[nData$group==g,1:p-1]),FUN=method))),2)
    aggData<-rbind(aggData,temp)
    #aggData<-rbind(aggData,cbind(temp,cluster=g))
  }
  rownames(aggData)<-levels(group)
  return(aggData)
}

isNaExist <- function(data){
  if(sum(is.na(data)*1)==0){#no na
    return(FALSE)
  }else{
    N<-dim(data)[1]
    completeN <- sum(complete.cases(data)*1)
    return(list(TRUE,N-completeN))
  }
}

finallySol <- function(data,numericData){
  pkg.env$intervalData<-RSDA::classic.to.sym(data)
  d<-RSDA::classic.to.sym(numericData)
  mind <- data.frame(matrix(sapply(d,min),nrow=dim(d)[1]))
  maxd <- data.frame(matrix(sapply(d,max),nrow=dim(d)[1]))
  colnames(mind)<-colnames(d);colnames(maxd)<-colnames(d)
  pkg.env$statisticsDF <- list(mind,maxd)
  names(pkg.env$statisticsDF) <- c("min","max")
}

buildRowName <- function(rData,iData,attrList){
  #change to factor type
  for(i in attrList){
    rData[,i]<-as.factor(rData[,i])
  }

  #all possible combination
  #1. build all attributes levels in a list
  attrLevel<-NULL
  for(i in attrList){
    attrLevel<-append(attrLevel,list(levels(rData[,i])))
  }
  #2 find combination by attributes levels list
  combination<-NULL
  if(length(attrLevel)>1){
    for(i in (length(attrLevel)-1):1){
      temp<-NULL
      for(u in 1:length(attrLevel[[i]])){
        if(i==(length(attrLevel)-1)){
          temp<-append(temp,paste(attrLevel[[i]][u],attrLevel[[i+1]],sep=","))
        }else{
          temp<-append(temp,paste(attrLevel[[i]][u],combination,sep=","))
        }
      }
      combination<-temp
    }
  }else{
    combination<-attrLevel[[1]]
    iData[,"rowname"]<-combination
    iData<-tibble::column_to_rownames(iData, var = "rowname")
    return(iData)
  }

  #data frame combination
  df<-rData[!duplicated(rData[,attrList]),attrList]

  #which combination exactly exist
  notIn <- NULL
  for(comb in combination){
    temp<-stringr::str_split(comb, ",")
    #if not exist
    if(is.na(prodlim::row.match(unlist(temp),df))){
      notIn<-append(notIn,comb)
    }
  }
  idx<-which(unlist(lapply(combination,function(x){x%in%notIn})))
  combination<-combination[-idx]

  #build rownames by combination
  rownames(iData)<-combination
  iData<-tibble::rownames_to_column(iData)
  iData<-tibble::column_to_rownames(iData, var = "rowname")
  return(iData)
}

customize <- function(nData=NULL,minList=NULL,maxList=NULL){
  if(is.null(nData)){
    stop("ERROR : Missing data input.")
  }
  nData<-as.data.frame(nData)
  if(is.null(minList)||is.null(maxList)){
    #default
    if(dim(nData)[2]==2){
      minList<-nData[,1]
      maxList<-nData[,2]
      warning("warning : Missing min data or max data.Automatically
              convert column 1 to min data,column 2 to max data.")
    }else{
      stop("ERROR : Missing min data or max data.")
    }
  }

  minList<-as.data.frame(minList)
  maxList<-as.data.frame(maxList)
  #test whether minList len == maxList
  if(dim(minList)[1]!=dim(maxList)[1] || dim(minList)[2]!=dim(maxList)[2]){
    stop("ERROR : Dimension Error,cannot match min and max.")
  }

  #test all min smaller than all max
  if(!all(minList<=maxList)){
    stop("ERROR : min values must smaller than all max values.")
  }

  p<-length(minList)
  n<-dim(minList)[1]
  allVar<-lapply(1:p,FUN=function(x) cbind(minList[x],maxList[x]))
  finalData <- matrix(NA,nrow = n)
  for(i in 1:p){
    temp<-unlist(lapply(1:n,FUN=function(x) cbind(allVar[[i]][[1]][x],allVar[[i]][[2]][x])))
    g=as.factor(rep(1:n,each=2))
    data<-data.frame(temp,g)
    data<-RSDA::classic.to.sym(data,g)
    colnames(data) <- paste0("V",i)
    finalData<-cbind(finalData,data)
  }
  finalData<-tibble::as.tibble(finalData)
  finalData<-finalData[,2:(p+1)]
  return(finalData)
}

buildStatsDf <- function(numericData = NULL){
  mind <- data.frame(matrix(sapply(numericData,min),nrow=dim(numericData)[1]))
  maxd <- data.frame(matrix(sapply(numericData,max),nrow=dim(numericData)[1]))
  mediand <- (mind+maxd)/2
  colnames(mind)<-colnames(numericData);colnames(maxd)<-colnames(numericData)
  colnames(mediand)<-colnames(numericData)
  return(list(mind,mediand,maxd))
}


pkg.env <- new.env()
pkg.env$statistics <- c("min","median","max","mean")
pkg.env$statisticsDF<-NULL
pkg.env$result<-NULL
pkg.env$intervalData<-NULL
pkg.env$rawData<-NULL

