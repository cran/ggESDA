#' @name ggInterval_radar
#' @title A interval Radar plot
#' @description  Using ggplot2 package to make a radar plot with multiple
#' variables.Each variables contains min values and max values as a
#' symbolic data.
#' @import ggplot2 ggforce
#' @importFrom gridExtra marrangeGrob
#' @importFrom RSDA is.sym.interval
#' @importFrom tidyr separate
#' @param data A ggESDA object. It can also be either RSDA object or
#' classical data frame(not recommended),which will be automatically convert to ggESDA
#' data.
#' @param layerNumber number of layer of a concentric circle,usually
#' to visuallize the reach of a observation in particularly variable.
#' @param inOneFig whether plot all observations in one figure.if not,
#' it will generate a new windows containing distinct observations.
#' @param showLegend whether show the legend.
#' @param showXYLabs whether show the x,y axis labels.
#' @param plotPartial a numeric vector,which is the row index from
#' the data.if it is not null, it will extract the row user deciding
#' to draw a radar plot from original data.Notes : the data must be
#' an interval data if the plotPartial is not null.
#' @param base_lty line type in base figure
#' @param base_circle boolean, if true, it will generate inner circle.
#' @param alpha aesthetic alpha of fill color
#' @param addText add the value of interval-valued variables in figure
#' @param type different type of radar,it can be "default","rect","quantile"
#' @param quantileNum if type == "quantile", it will provide the number of percentage
#' @param Drift The drift term, which determines the radar values beginning.
#' @param addText_modal add the factor of modal multi-valued variables in figure..
#' @param addText_modal.p add the value of modal multi-valued variables in figure..
#' @usage ggInterval_radar(data=NULL,layerNumber=3,
#' inOneFig=TRUE,showLegend=TRUE,showXYLabs=FALSE,
#' plotPartial=NULL,
#' alpha=0.5,
#' base_circle=TRUE,
#' base_lty=2,
#' addText=TRUE,
#' type="default",
#' quantileNum=4,
#' Drift=0.5,
#' addText_modal=TRUE,
#' addText_modal.p=FALSE)
#' @examples
#' # must specify plotPartial to some certain rows you want to plot
#' Environment.n <- Environment[, 5:17]
#' ggInterval_radar(Environment.n,
#'                 plotPartial = 2,
#'                 showLegend = FALSE,
#'                 base_circle = TRUE,
#'                 base_lty = 2,
#'                 addText = FALSE
#') +
#'  labs(title = "") +
#'  scale_fill_manual(values = c("gray50")) +
#'  scale_color_manual(values = c("red"))
#'
#' ggInterval_radar(Environment,
#'                 plotPartial = 2,
#'                 showLegend = FALSE,
#'                 base_circle = FALSE,
#'                 base_lty = 1,
#'                 addText = TRUE
#') +
#'  labs(title = "") +
#'  scale_fill_manual(values = c("gray50")) +
#'  scale_color_manual(values = c("gray50"))
#'
#' @export
ggInterval_radar <-function(data=NULL,layerNumber=3,
                            inOneFig=TRUE,showLegend=TRUE,showXYLabs=FALSE,
                            plotPartial=NULL,
                            alpha=0.5,
                            base_circle=TRUE,
                            base_lty=2,
                            addText=TRUE,
                            type="default",
                            quantileNum=4,
                            Drift=0.5,
                            addText_modal=TRUE,
                            addText_modal.p=FALSE){
  #globalVariables("varLevels", add = F)
  #utils::globalVariables(c("varLevels"), add=FALSE)
  varLevels <- NULL
  extendUnit <- Drift
  if(extendUnit < 0 | extendUnit > 1){
    stop("Drift must between 0 to 1.")
  }

  #extend plot original xyLimits = 1.25, extendUnit=0, colnamesDrift = 0.1
  colnamesDrift <- 0.18
  xyLimits <- extendUnit + 1.25 + colnamesDrift
  textShift<-0.065
  textSize <- 4
  textCol <- "gray20"
  #

  fillBetween=TRUE #not fix complete
  #notes
  if(dim(data)[1]<=1){
    stop("ERROR IN data : data must be full data with observations larger than 1.")
  }
  #setting parameter
  nL <- layerNumber
  if(is.null(plotPartial)){
    #test data to iData
    ggSymData <- testData(data)
    iData <- ggSymData$intervalData
    indNum <-dim(iData)[1]
    plotPartial <- c(1:indNum)
  }else{
    #checking input format
    if(!"symbolic_tbl"%in%class(data)){
      stop("ERROR : data must be symbolic data,when plotPartial is on.")
    }
    if(!is.numeric(plotPartial)){
      stop("ERROR IN plotPartial : the input data must be numeric vector as the row index in full data.")
    }
    if(!all(plotPartial<=nrow(data)) || !all(plotPartial>0)){
      stop("ERROR IN plotPartial : the input data must be numeric vector as the row index in full data.")
    }
    iData <- data
    indNum <- length(plotPartial)
  }

  if( !(type%in%c("default","rect","quantile")) ){
    warning(paste("There is no type called ",type,", automatically set default type."))
    type<-"default"
  }

  #data preprocessing
  numericData <- unlist(lapply(iData[,1:dim(data)[2]] ,FUN = RSDA::is.sym.interval))

  rawiData<-iData
  allnP <- dim(rawiData)[2]
  propData <-iData[,!numericData]
  iData <- iData[,numericData]
  if("symbolic_tbl" %in% class(iData)){
    if(dim(iData)[2] == 0){
      #all data are symbolic modal
      nP <- 0
    }else{
      #normal type , there are some symbolic interval data (len>1)
      nP <- dim(iData)[2]
    }
  }else{
    #only one symbolic interval data
    nP <- 1
    iData <- rawiData[, c(which(numericData), which(!numericData))]
  }

  #someday i will fix this bug ,maybe
  # if(nP < 2){
  #   stop("number of numerical (symbolic interval) data must greater than 2")
  # }


  #setting quantile argument
  if(type=="quantile"){
    #warning("inOneFig will be unused in type=='quantile'.")
    inOneFig<-TRUE
    indNum <-dim(iData)[1]
    if(!quantileNum%in%c(2:10)){
      stop("the number of quantile must be in the range [2,10].")
    }
    if(nP<3){
      stop("numerical variables must greater than 3 in type=='quantile'")
    }
    if(allnP!=nP){
      warning("nominal variables will be ignore in type=='quantile'")
      allnP<-nP
    }
    plotPartial<-NULL
  }

  #test ERROR
  if(type!="quantile"){
    if(indNum>50 & inOneFig==F){
      stop("Out of time limits")
    }else if(indNum>10 & inOneFig==T){
      stop("Suggest set inOneFig to FALSE when the observations are large.")
    }else if(indNum>5 & inOneFig==T){
      warning("Suggest set inOneFig to FALSE when the observations are large.")
    }
  }

  #setting original plot
  #p<-generateCircle(nL)+coord_fixed(ratio = 1)

  if(base_circle){
    p<-generateCircle(nL+1,extendUnit)+coord_fixed(ratio = 1)
  }else{
    p<-ggplot()+coord_fixed(ratio = 1)
  }

  if(allnP!=nP){
    d<-generatePoint(allnP,nL,extendUnit)
    d<-d[(allnP*nL-allnP+1):(allnP*nL),]
    rawiData<-cbind(iData,propData)
    propData<-propData[plotPartial,]
  }else{
    rawiData<-iData
    d<-generatePoint(nP,nL,extendUnit)
    d<-d[(nP*nL-nP+1):(nP*nL),]
  }
  p<-p+geom_point(data=d,aes(x=d$x,y=d$y))
  #get all variables min max data
  if(nP != 0){
    allData<-iData[,1:nP]


    if(!is.null(plotPartial)){
      iData<-iData[plotPartial,1:nP]
    }
    dataList<-lapply(1:nP,FUN=function(x) data.frame(allData[[x]]))
    iDataList<-lapply(1:nP,FUN=function(x) data.frame(iData[[x]]))
    maxList<-lapply(1:nP,FUN=function(x) max(dataList[[x]]))
    minList<-lapply(1:nP,FUN=function(x) min(dataList[[x]]))
    minminList<-minList;  maxmaxList<-maxList

    #quantile #this will over write iDataList
    if(type=="quantile"){
      indNum<-quantileNum+1
      figQuaNum<-quantileNum+3
      sepQua<-seq(0,1,1/figQuaNum)
      for(i in 1:nP){
        datatmp<-c(iDataList[[i]]$min,iDataList[[i]]$max)
        #extend adjust
        datatmp<-round(quantile(datatmp,sepQua),2)[-(figQuaNum+1)][-1] + extendUnit
        #end adjust
        tempMin<-datatmp[1:indNum]
        tempMax<-datatmp[2:(figQuaNum-1)]
        iDataList[[i]]<-data.frame(min=tempMin,max=tempMax)
      }
    }

    if(type == "quantile"){
      maxList<-lapply(1:nP,FUN=function(x) max(iDataList[[x]][, "max"]))
      minList<-lapply(1:nP,FUN=function(x) min(iDataList[[x]][, "min"]))
      minminList<-minList;  maxmaxList<-maxList
    }

    #normalize data to 0,1
    #extend adjust
    normData<-lapply(1:nP ,FUN=function(x){
      sapply(iDataList[[x]],FUN=function(elem) {
        ( (elem-minList[[x]])/(maxList[[x]]-minList[[x]]) ) + extendUnit
      })
    })




    #rescale dataframe to input form which data2Vec need
    minDF<-sapply(1:nP,FUN=function(x) matrix(normData[[x]],ncol=2)[,1])
    maxDF<-sapply(1:nP,FUN=function(x) matrix(normData[[x]],ncol=2)[,2])

    minDF<-matrix(minDF,ncol=nP)
    maxDF<-matrix(maxDF,ncol=nP)
    transMatrix<-d[1:nP,c(3,4)]



    #convert data to match radar axis scale
    if(type=="quantile"){
      setRowNameDf <- iDataList[[1]]
    }else{
      setRowNameDf<-iData
    }
    plotMin<-data2Vec(iData=setRowNameDf,data=minDF,transMat = transMatrix,type,quantileNum)
    plotMax<-data2Vec(iData=setRowNameDf,data=maxDF,transMat = transMatrix,type,quantileNum)


    #making labels
    temp<-lapply(1:nP,FUN=function(x) paste(iDataList[[x]][,1],iDataList[[x]][,2],sep=":"))
    temp<-lapply(1:nP,FUN=function(x) paste0("[",temp[[x]],"]"))
    temp<-matrix(unlist(temp),ncol=indNum,byrow=T)
    groupId<-rep(c(1:indNum),each=nP)
    myLabel=paste(groupId,paste(colnames(iData),temp),sep=" : ")


    #new a labels variable
    plotMin<-data.frame(plotMin,Variables=myLabel)
    plotMax<-data.frame(plotMax,Variables=myLabel)

    #make text in plot
    # minList<-unlist(lapply(1:nP,FUN=function(x) iDataList[[x]][,1]))
    # maxList<-unlist(lapply(1:nP,FUN=function(x) iDataList[[x]][,2]))
    #
    newDf <- data.frame(NULL)
    for(i in 1:indNum){
      for(u in 1:nP){
        newDf<-rbind(newDf,iDataList[[u]][i,])
      }
    }

    textMin<-as.data.frame(cbind(plotMin,min=newDf$min))
    textMax<-as.data.frame(cbind(plotMax,max=newDf$max))
    textMin$cos<-textMin$cos+textShift;textMin$sin<-textMin$sin+textShift
    textMax$cos<-textMax$cos+textShift;textMax$sin<-textMax$sin+textShift

    plotMin$group<-as.factor(plotMin$group)
  }#end if np != 0

  #calculate nominal variable (point)
  if(allnP!=nP){
    # preparing the nominal data set
    tmpDf<-data.frame(NULL)
    varNum <- 1
    for(var in colnames(propData)){
      propDf<-data.frame(NULL)
      counter<-1
      for(ele in 1:indNum){
        varPro<-unlist(propData[ele,var])
        propDf[counter:(counter+(length(varPro)/2)-1),"varName"] <- rep(var,length(varPro)/2)
        propDf[counter:(counter+(length(varPro)/2)-1),"varLevels"]<-paste(varPro[1:(length(varPro)/2)],round(as.numeric(varPro[(length(varPro)/2+1):length(varPro)]),2),sep=":")
        propDf[counter:(counter+(length(varPro)/2)-1),"prop"]<-as.numeric(varPro[(length(varPro)/2+1):length(varPro)])
        propDf[counter:(counter+(length(varPro)/2)-1),"groupid"] <- rep(ele,length(varPro)/2)
        #extendUnit adjust  quantile(seq(0,1,1/3),seq(0,1,1/3))[2:(3+1)]
        #pos<-seq(0,1,1/((length(varPro)+2)/2))
        pos <- as.numeric(quantile(seq(0,1+extendUnit,1/((length(varPro)+2)/2)),seq(0,1,1/((length(varPro)+2)/2)))[2:(length(varPro)/2+1)])
        #pos <- pos[-length(pos)][-1]
        #end adjust
        propDf[counter:(counter+(length(varPro)/2)-1),"x"] <- d[nP+varNum,3] * pos
        propDf[counter:(counter+(length(varPro)/2)-1),"y"] <- d[nP+varNum,4] * pos
        counter<-counter+length(varPro)/2
      }
      tmpDf<-rbind(tmpDf,propDf)
      varNum<-varNum+1
    }
    propDf<-as.data.frame(tmpDf)
    propDf$groupid <- as.factor(propDf$groupid)


    # preparing the symbolic(rectangle) to represent nominal data
    totalRectDf <- data.frame(NULL)
    obsG<-1
    #lastRect <- NULL
    heiList <- list(NULL)
    for(i in 1:indNum){
      heiList[[i]]<-0
    }
    for(ele in levels(propDf$groupid)){
      propDf.ele <- dplyr::filter(propDf,propDf$groupid==ele)
      varG <-1
      tempRectDf <- data.frame(NULL)
      for(i in unique(as.factor(propDf$varName))){
        rectDf <- data.frame(NULL)
        propDf.temp <- dplyr::filter(propDf.ele,propDf.ele$varName==i)
        varL <- 1
        heiVec <- NULL
        for(u in 1:dim(propDf.temp)[1]){
          #re range 0~1 to -0.9~0.5 (1 = 平面)
          #print(paste0(ele,i,u))
          #print(propDf.temp)
          thisHei<-reRange(-0.99,-0.8,propDf.temp[u,"prop"])
          heiVec<-cbind(heiVec,thisHei)
          rect <- build3DRect(-0.95,-0.95,thisHei,g=paste0(i,u,ele))
          #print(rect)
          distx <- mean(rect$newx-propDf.temp[u,"x"])
          disty <- mean(rect$newy[13:16]-propDf.temp[u,"y"])
          rect$newx <- rect$newx - distx
          rect$newy <- rect$newy - disty

          #add in height
          if(obsG>=2&&inOneFig){
            rect$newy <- rect$newy + unlist(heiList[[obsG-1]][varG])[u]
          }

          rect[,"varLevels"]<-varL
          rectDf<-rbind(rectDf,rect)
          varL<-varL+1
        }#inner for
        heiList[[obsG]][varG]<-list(1+heiVec)
        if(obsG>=2){
          heiList[[obsG]][varG]<-list(unlist(heiList[[obsG-1]][varG])+unlist(heiList[[obsG]][varG]))
        }
        rectDf[,"varGroup"]<-varG
        tempRectDf<-rbind(tempRectDf,rectDf)
        varG<-varG+1
      }# sec for
      tempRectDf[,"obsGroup"]<-obsG
      obsG<-obsG+1
      # if(!is.null(lastRect)){
      #   tempRectDf$newy <- lastRect$newy+tempRectDf$newy
      # }
      totalRectDf<-rbind(totalRectDf,tempRectDf)
      #lastRect<-tempRectDf
      #print(lastRect)
    }# first for


    #totalRectDf<-as.data.frame(totalRectDf)
    totalRectDf$varLevels <- as.factor(totalRectDf$varLevels)
    totalRectDf$varGroup <- as.factor(totalRectDf$varGroup)
    totalRectDf$obsGroup <- as.factor(totalRectDf$obsGroup)
    if(nP != 0){
      levels(propDf$groupid)<-levels(plotMin$group)
    }
  }

  #generate cut line for type==rect && build nominal rect
  if(type=="rect"){
    if(nP != 0){
      cutDf<-c(getCutDf(plotMin,transMatrix,nP,indNum),getCutDf(plotMax,transMatrix,nP,indNum))
    }
  }

  #plot
  if(inOneFig){
    if(type=="default"||type=="quantile"){
      #if(fillBetween){
      myPolyData<-data.frame(NULL)
      if(nP != 0){
        for(i in levels(as.factor(plotMin$group))){
          plotMin.temp <- dplyr::filter(plotMin,plotMin$group==i)
          plotMax.temp <- dplyr::filter(plotMax,plotMin$group==i)
          myPathData<-data.frame(x1=plotMin.temp$cos,y1=plotMin.temp$sin,
                                 x2=plotMax.temp$cos,y2=plotMax.temp$sin)
          if(allnP==nP){
            newTemp<-rbind(plotMin.temp,plotMin.temp[1,],plotMax.temp,plotMax.temp[1,])
            myPolyData<-rbind(myPolyData,newTemp)
          }else{ #add nominal
            tmp2<-NULL;j<-1
            # print(levels(as.factor(propDf$varName)))
            # print(propDf)
            for(u in unique(as.factor(propDf$varName))){
              # propDf.temp <- dplyr::filter(propDf,propDf$varName==u) %>%
              #   dplyr::filter(groupid==i)

              propDf.temp <- dplyr::filter(propDf,propDf$varName==u)
              propDf.temp <- dplyr::filter(propDf.temp,propDf.temp$groupid==i)
              tmp<-propDf.temp[propDf.temp$prop==max(propDf.temp$prop),]
              tmp<-tmp[1,]
              tmp2<-rbind(tmp2,plotMin.temp[1,])
              tmp2[j,c("cos","sin")]<-c(tmp$x,tmp$y)
              myPathData<-rbind(myPathData,data.frame(x1=tmp$x,y1=tmp$y,x2=tmp$x,y2=tmp$y))
              j<-j+1
            }
            #propDf.temp <- dplyr::filter(propDf,groupid==i)
            #tmp<-propDf[propDf.temp$prop==max(propDf.temp$prop),]
            #tmp<-tmp[1,] #還沒想到一次連兩個點的方法 先只取第一個
            #tmp2<-plotMin.temp[1,]
            #tmp2[,c("cos","sin")]<-c(tmp$x,tmp$y)
            #myPathData<-rbind(myPathData,data.frame(x1=tmp$x,y1=tmp$y,x2=tmp$x,y2=tmp$y))

            newTemp<-rbind(plotMin.temp,tmp2,
                           plotMin.temp[1,], plotMax.temp,tmp2,plotMax.temp[1,])
            myPolyData<-rbind(myPolyData,newTemp)
          }
          p<-p+geom_path(data=myPathData,aes(x=myPathData$x1, y=myPathData$y1),lty=0)+
            geom_path(data=myPathData,aes(x=myPathData$x2, y=myPathData$y2),lty=0)
        }
      }

      if(allnP!=nP){
        if(nP != 0){
          tmpN<-dim(myPolyData)[1]/indNum
          myPolyData[,"obsGroup"]<-as.factor(rep(1:indNum,each=tmpN))
          tempPolyDf <- myPolyData[,c(1,2,3,5)]
        }

        tempRectDf <- totalRectDf[,c(1,2,3,6)]

        colnames(tempRectDf) <- c("cos", "sin", "group", "obsGroup")
        #colnames(tempRectDf)<-colnames(tempPolyDf)
        if(nP != 0){
          polyDf<-as.data.frame(rbind(tempPolyDf,tempRectDf))
        }else{
          polyDf <- tempRectDf
        }

        p<-p+geom_polygon(data=polyDf,aes(x=polyDf$cos,y=polyDf$sin,group=polyDf$group,fill=polyDf$obsGroup,col=polyDf$obsGroup),
                          alpha=alpha)
        if(nP != 0){
          p <- p +  geom_point(data=plotMin,aes(x=plotMin$cos,y=plotMin$sin))+
            geom_point(data=plotMax,aes(x=plotMax$cos,y=plotMax$sin))
        }

      }else{
        if(type=="quantile"){
          p<-p+geom_polygon(data=myPolyData,aes(x=myPolyData$cos,y=myPolyData$sin,fill=myPolyData$group,col=myPolyData$group),alpha = alpha)
        }else{
          p<-p+geom_polygon(data=myPolyData,aes(x=myPolyData$cos,y=myPolyData$sin,fill=myPolyData$group,col=myPolyData$group),alpha = alpha,lty=0)
        }

        if(type!="quantile"){
          p<-p+geom_point(data=plotMin,aes(x=plotMin$cos,y=plotMin$sin))+
            geom_point(data=plotMax,aes(x=plotMax$cos,y=plotMax$sin))
        }

      }

    }
    #else{ #this is else for if(fillbetween)
    else if(type=="rect"){
      #add cut segment
      if(nP != 0){
        rectPolyData<-data.frame(NULL)
        gId <- 1
        for(g in levels(as.factor(cutDf[[1]]$obsGroup))){
          cutDf.temp<-lapply(1:4,FUN = function(x){as.data.frame(dplyr::filter(cutDf[[x]],cutDf[[x]]$obsGroup==g))})
          #print(cutDf.temp)
          for(i in 1:nP){
            tmpDf<-data.frame(x=c(cutDf.temp[[1]]$x2[i],cutDf.temp[[2]]$x2[i],cutDf.temp[[4]]$x2[i],cutDf.temp[[3]]$x2[i]),
                              y=c(cutDf.temp[[1]]$y2[i],cutDf.temp[[2]]$y2[i],cutDf.temp[[4]]$y2[i],cutDf.temp[[3]]$y2[i]),
                              varGroup=as.factor(rep(gId,4)),
                              obsGroup=as.factor(g))
            rectPolyData<-rbind(rectPolyData,tmpDf)
            gId<-gId+1
          }
        }
      }
      #print(cutDf)
      #print(rectPolyData)
      #numeric
      if(allnP!=nP){
        tempRectDf<-totalRectDf[,c(1,2,3,6)]
        #colnames(tempRectDf) <- colnames(rectPolyData)
        colnames(tempRectDf) <- c("x", "y", "varGroup", "obsGroup")
        if(nP != 0){
          levels(tempRectDf$obsGroup) <- levels(rectPolyData$obsGroup)
          rectPolyData<-as.data.frame(rbind(tempRectDf,rectPolyData))
        }else{
          rectPolyData <- tempRectDf
        }

      }
      p<-p+geom_polygon(data=rectPolyData,aes(x=.data$x,y=.data$y,group=.data$varGroup,
                                              fill=.data$obsGroup,col=.data$obsGroup),alpha=alpha)


    }
    #nominal
    if(allnP!=nP){
      newPropDf<-as.data.frame(tidyr::separate(propDf,varLevels,into=c("varLevels","gar"),sep=":",convert = F))
      newPropDf<-as.data.frame(dplyr::filter(newPropDf,newPropDf$groupid==levels(newPropDf$groupid)[1]))

      if(addText_modal){
        p<-p+geom_text(data=newPropDf,aes(x=newPropDf$x+textShift,y=newPropDf$y+textShift,label=newPropDf$varLevels),vjust=2.15,hjust=0.5,
                       size = textSize, col = textCol)
      }
      #print(totalRectDf)

      if(addText_modal.p){
        tempd<-data.frame(NULL)
        for(i in unique(totalRectDf$varGroup)){
          temp.varG<-dplyr::filter(totalRectDf,totalRectDf$varGroup==i)
          for(u in levels(temp.varG$obsGroup)){
            temp.obsG<-dplyr::filter(temp.varG,temp.varG$obsGroup==u)
            for(k in unique(temp.obsG$varLevels)){
              temp.varL<-dplyr::filter(temp.obsG,temp.obsG$varLevels==k)
              tempd<-rbind(tempd,data.frame(x=mean(temp.varL$newx),y=mean(temp.varL$newy)))
            }
          }
        }
        propTextDf<-data.frame(tempd,prop=propDf$prop)
        #print(propTextDf)
        p<-p+geom_text(data=propTextDf,aes(x=propTextDf$x+textShift,y=propTextDf$y,label=round(propTextDf$prop,2)),
                       size = textSize, col = textCol)
      }
    }
    #print("5")
    # if(allnP!=nP){
    #   p<-p+geom_polygon(data=totalRectDf,aes(x=newx,y=newy,group=group,fill=obsGroup),
    #                     col="black",
    #                     alpha=alpha)+
    #     geom_text(data=propDf,aes(x=propDf$x+textShift,y=propDf$y+textShift,label=varLevels),
    #               vjust=2.75)
    # }

    #colnames shift out of figure 0.1 unit (colnamesDrift)
    newD <- shift(d, unit = colnamesDrift)
    p<-p+geom_segment(data=d,aes(x=0,y=0,xend=d$x,yend=d$y),lty=base_lty,alpha=0.6, col = "gray50")
    if(nP != 0){
      p <- p + geom_point(data=plotMax,aes(x=0,y=0,alpha=plotMax$Variables))+
        geom_text(data=newD,aes(x=newD$x,y=newD$y,label=c(colnames(rawiData))))

    }



    if(addText && type!="quantile" && nP != 0){
      p<-p+geom_text(data=textMin,aes(x=textMin$cos,y=textMin$sin,label=textMin$min),
                     size = textSize, col = textCol)+
        geom_text(data=textMax,aes(x=textMax$cos,y=textMax$sin,label=textMax$max),
                  size = textSize, col = textCol)
    }
    if(!showXYLabs){
      p<-p+scale_x_continuous(labels =NULL,limits = c(-xyLimits,xyLimits))+
        scale_y_continuous(labels =NULL,limits = c(-xyLimits,xyLimits))+xlab(NULL)+ylab(NULL)
    }else{
      #adds scale y continuous label
      yLabelPos<-seq(-xyLimits,xyLimits,(xyLimits-(-xyLimits))/(nP+1))
      yLabelPos<-yLabelPos[-1]
      yLabelTitle<-paste0(colnames(iData),": [")
      yLabels<-paste0(yLabelTitle,paste(unlist(minminList),unlist(maxmaxList),sep=","),"]")
      yLabels<-c(yLabels,expression(bold("The range of each variable show as follows :")))
      p<-p+scale_x_continuous(limits = c(-xyLimits,xyLimits))+
        scale_y_continuous(breaks=yLabelPos,labels=yLabels,limits = c(-xyLimits,xyLimits))+xlab(NULL)+ylab(NULL)+
        scale_x_continuous(labels=NULL,limits = c(-xyLimits,xyLimits))
    }
    #showLegend
    if(showLegend){
      p<-p+guides(col=F)
    }else{
      p<-p+guides(fill=F,col=F,alpha=F)
    }
    if(type=="quantile"){
      myTitle<-"Radar : quantile plot"
    }else{
      #myTitle<-paste0("Radar : ",paste(sapply(rownames(iData), paste, collapse=":"), collapse=","))
      myTitle = "Radar plot"
    }
    p<-p+labs(title=myTitle)+
      scale_colour_discrete(name = "Group")+
      scale_alpha_discrete(name="Interval value of each observation")+
      scale_fill_discrete(name = "Group")+
      theme_bw()

    #make circle interpret
    if(base_circle){
      a<-cos(120*pi/180)
      b<-sin(120*pi/180)
      #extend adjust quantile(seq(0,myRange,1/nLayer),seq(0,1,1/nLayer))[2:(nLayer+1)]
      tmp<-quantile(seq(0+extendUnit,1+extendUnit,1/layerNumber),seq(0,1,1/layerNumber))
      #tmp<-seq(0,1,1/layerNumber)[-1]
      #end adjust
      tmp<-data.frame(x=a*tmp,y=b*tmp,label=paste0(as.character(round((tmp-extendUnit)*100)),"%"))
      p<-p+geom_text(data=tmp,aes(x=tmp$x,y=tmp$y,label=tmp$label))
    }
  }else{
    plotList <- NULL
    u<-1
    for(i in levels(as.factor(plotMin$group))){
      #print(i)
      if(nP == 0){
        stop("Cannot allows discrete variables only, under the parameter inOneFig = FALSE,
             please adjust inOneFig to TRUE.")
      }
      plotMin.temp <- dplyr::filter(plotMin,plotMin$group==i)
      plotMax.temp <- dplyr::filter(plotMax,plotMax$group==i)
      textMin.temp <- dplyr::filter(textMin,textMin$group==i)
      textMax.temp <- dplyr::filter(textMax,textMax$group==i)
      if(allnP!=nP){
        propDf.temp <- dplyr::filter(propDf,propDf$groupid==levels(propDf$groupid)[u])
        totalRectDf.temp <- dplyr::filter(totalRectDf,totalRectDf$obsGroup==levels(totalRectDf$obsGroup)[u])
      }else{
        propDf.temp<-0
        totalRectDf.temp<-0
      }

      if(type=="rect"){
        cutDf.temp<-list(NULL)
        for(k in 1:4){
          cutDf[[k]]$obsGroup <- as.factor(cutDf[[k]]$obsGroup)
          cutDf.temp[[k]]<-dplyr::filter(cutDf[[k]],cutDf[[k]]$obsGroup==levels(cutDf[[k]]$obsGroup)[u])
        }
      }else{
        cutDf.temp<-0
      }



      #print(paste(u,plotMin.temp))
      base<-plotFun(p,iData,plotMin.temp,plotMax.temp,d,showXYLabs,showLegend,fillBetween,base_circle,layerNumber,
                    textMin.temp,textMax.temp,rawiData,propDf,propDf.temp,allnP,nP,type,cutDf.temp,totalRectDf.temp,addText,textShift,base_lty,alpha,xyLimits,extendUnit,addText_modal,addText_modal.p)


      base<-base+labs(title=paste0("Radar : ",rownames(iData)[u]))+labs(title=paste0("Radar : ",rownames(iData)[u]))+scale_colour_discrete(name = "Group")+
        scale_alpha_discrete(name="Interval")
      plotList[[u]]<-base
      u<-u+1
    }

    #print(plotList[[1]])
    #p<-gridExtra::marrangeGrob(plotList,ncol=2,nrow=1)
    myRow <- round(sqrt(indNum))
    myCol <- ceiling(indNum/myRow)
    p <- gridExtra::marrangeGrob(plotList,
                   ncol = myCol,
                   nrow = myRow)

  }

  if(type=="quantile"){
    colShift<-4
    myValues <- rev(grey.colors(quantileNum+1+colShift))[-c(1:colShift)]
    p<-p+guides(alpha=F)+
      scale_fill_manual(name="quantile percentage",values=myValues)+
      scale_colour_manual(values=myValues)
  }
  return(p)
}




plotFun<-function(p,iData,plotMin.temp,plotMax.temp,d,showXYLabs,showLegend,fillBetween,base_circle,layerNumber,textMin.temp,textMax.temp,rawiData,propDf,propDf.temp,allnP,nP,type,cutDf,totalRectDf.temp,addText,textShift,base_lty,alpha,xyLimits,extendUnit,addText_modal,addText_modal.p){
  #if(fillBetween){
  varLevels <- NULL
  if(type=="default"){
    myPolyData<-data.frame(NULL)
    myPathData<-data.frame(x1=plotMin.temp$cos,y1=plotMin.temp$sin,
                           x2=plotMax.temp$cos,y2=plotMax.temp$sin)

    if(allnP==nP){
      myPolyData<-rbind(plotMin.temp,plotMin.temp[1,],plotMax.temp,plotMax.temp[1,])

    }else{ #add nominal
      tmp2<-NULL;j<-1
      # print(levels(as.factor(propDf$varName)))
      # print(propDf)
      for(u in unique(as.factor(propDf.temp$varName))){
        propDf.temp2 <- dplyr::filter(propDf.temp,propDf.temp$varName==u)
        tmp<-propDf.temp2[propDf.temp2$prop==max(propDf.temp2$prop),]
        tmp<-tmp[1,]
        tmp2<-rbind(tmp2,plotMin.temp[1,])
        tmp2[j,c("cos","sin")]<-c(tmp$x,tmp$y)
        myPathData<-rbind(myPathData,data.frame(x1=tmp$x,y1=tmp$y,x2=tmp$x,y2=tmp$y))
        j<-j+1
      }
      newTemp<-rbind(plotMin.temp,tmp2,
                     plotMin.temp[1,], plotMax.temp,tmp2,plotMax.temp[1,])
      myPolyData<-rbind(myPolyData,newTemp)
    }
    base<-p+geom_path(data=myPathData,aes(x=myPathData$x1, y=myPathData$y1))+
      geom_path(data=myPathData,aes(x=myPathData$x2, y=myPathData$y2))



    if(allnP!=nP){
      #tmpN<-dim(myPolyData)[1]

      this.obs <- totalRectDf.temp$obsGroup[1]
      myPolyData[,"obsGroup"]<-as.factor(this.obs)
      tempRectDf <- totalRectDf.temp[,c(1,2,3,6)]

      tempPolyDf <- myPolyData[,c(1,2,3,5)]
      colnames(tempRectDf)<-colnames(tempPolyDf)
      polyDf<-as.data.frame(rbind(tempPolyDf,tempRectDf))

      base<-base+
        geom_point(data=plotMin.temp,aes(x=plotMin.temp$cos,y=plotMin.temp$sin))+
        geom_point(data=plotMax.temp,aes(x=plotMax.temp$cos,y=plotMax.temp$sin))




      base<-base+geom_polygon(data=polyDf,aes(x=polyDf$cos,y=polyDf$sin,group=polyDf$group,
                                              fill=polyDf$obsGroup,col=polyDf$obsGroup),
                              alpha=alpha)



    }else{
      base<-base+geom_polygon(data=myPolyData,aes(x=myPolyData$cos,y=myPolyData$sin,fill=myPolyData$group,col=myPolyData$group),alpha = alpha)+
        geom_point(data=plotMin.temp,aes(x=plotMin.temp$cos,y=plotMin.temp$sin))+
        geom_point(data=plotMax.temp,aes(x=plotMax.temp$cos,y=plotMax.temp$sin))
    }

    # base<-base+geom_polygon(data=myPolyData,aes(x=myPolyData$cos,y=myPolyData$sin,fill=group),alpha = alpha,col="black")
  }else if(type=="rect"){
    #add cut segment
    rectPolyData<-data.frame(NULL)
    gId <- 1
    for(g in levels(as.factor(cutDf[[1]]$obsGroup))){
      cutDf.temp<-lapply(1:4,FUN = function(x){as.data.frame(dplyr::filter(cutDf[[x]],cutDf[[x]]$obsGroup==g))})
      #print(cutDf.temp)
      for(i in 1:nP){
        tmpDf<-data.frame(x=c(cutDf.temp[[1]]$x2[i],cutDf.temp[[2]]$x2[i],cutDf.temp[[4]]$x2[i],cutDf.temp[[3]]$x2[i]),
                          y=c(cutDf.temp[[1]]$y2[i],cutDf.temp[[2]]$y2[i],cutDf.temp[[4]]$y2[i],cutDf.temp[[3]]$y2[i]),
                          varGroup=as.factor(rep(gId,4)),
                          obsGroup=as.factor(g))
        rectPolyData<-rbind(rectPolyData,tmpDf)
        gId<-gId+1
      }

    }
    #print(cutDf)
    #print(rectPolyData)
    #numeric
    if(allnP!=nP){
      tempRectDf<-totalRectDf.temp[,c(1,2,3,6)]
      colnames(tempRectDf) <- colnames(rectPolyData)
      levels(tempRectDf$obsGroup) <- levels(rectPolyData$obsGroup)
      rectPolyData<-as.data.frame(rbind(tempRectDf,rectPolyData))
    }
    base<-p+geom_polygon(data=rectPolyData,aes(x=.data$x,y=.data$y,group=.data$varGroup,
                                               fill=.data$obsGroup,col=.data$obsGroup),alpha=alpha)

  }#end rect

  #add nominal text
  if(allnP!=nP){

    newPropDf<-as.data.frame(tidyr::separate(propDf,varLevels,into=c("varLevels","gar"),sep=":",convert = F))
    newPropDf<-as.data.frame(dplyr::filter(newPropDf,newPropDf$groupid==levels(newPropDf$groupid)[1]))
    if(addText_modal){
      base<-base+geom_text(data=newPropDf,aes(x=newPropDf$x+textShift,y=newPropDf$y+textShift,label=newPropDf$varLevels),vjust=2.75)
    }

    #print(totalRectDf.temp)
    if(addText_modal.p){
      tempd<-data.frame(NULL)
      for(i in levels(totalRectDf.temp$varGroup)){
        temp.varG<-dplyr::filter(totalRectDf.temp,totalRectDf.temp$varGroup==i)
        for(k in unique(temp.varG$varLevels)){
          temp.varL<-dplyr::filter(temp.varG,temp.varG$varLevels==k)
          tempd<-rbind(tempd,data.frame(x=mean(temp.varL$newx),y=mean(temp.varL$newy)))
        }
      }
      propTextDf<-data.frame(tempd,prop=propDf.temp$prop)
      #print(propTextDf)
      base<-base+geom_text(data=propTextDf,aes(x=propTextDf$x+textShift,y=propTextDf$y,label=round(propTextDf$prop,2)))
    }
  }

  newD <- shift(d,0.1)
  base<-base+geom_segment(data=d,aes(x=0,y=0,xend=d$x,yend=d$y),lty=base_lty,alpha=0.6,col="gray50")+
    geom_point(data=plotMax.temp,aes(x=0,y=0,alpha=plotMax.temp$Variables))+
    geom_text(data=newD,aes(x=newD$x,y=newD$y,label=c(colnames(rawiData))))

  if(addText){
    base <- base+geom_text(data=textMin.temp,aes(x=textMin.temp$cos,y=textMin.temp$sin,label=textMin.temp$min))+
      geom_text(data=textMax.temp,aes(x=textMax.temp$cos,y=textMax.temp$sin,label=textMax.temp$max))
  }

  if(!showXYLabs){
    base<-base+scale_x_continuous(labels =NULL,limits = c(-xyLimits,xyLimits))+
      scale_y_continuous(labels =NULL,limits = c(-xyLimits,xyLimits))+xlab(NULL)+ylab(NULL)
  }else{
    base<-base+scale_x_continuous(limits = c(-xyLimits,xyLimits))+
      scale_y_continuous(limits = c(-xyLimits,xyLimits))
  }#showLegend
  if(showLegend){
    base<-base+guides(fill=F,col=F)
  }else{
    base<-base+guides(fill=F,col=F,alpha=F)
  }

  if(base_circle){
    a<-cos(120*pi/180)
    b<-sin(120*pi/180)

    #extend adjust quantile(seq(0,myRange,1/nLayer),seq(0,1,1/nLayer))[2:(nLayer+1)]
    tmp<-quantile(seq(0+extendUnit,1+extendUnit,1/layerNumber),seq(0,1,1/layerNumber))
    #tmp<-seq(0,1,1/layerNumber)[-1]
    #end adjust
    tmp<-data.frame(x=a*tmp,y=b*tmp,label=paste0(as.character(round((tmp-extendUnit)*100)),"%"))
    base<-base+geom_text(data=tmp,aes(x=tmp$x,y=tmp$y,label=tmp$label))
  }
  base<-base+theme_bw()+guides(alpha = F)
  return(base)
}
data2Vec <- function(iData=NULL,data=NULL,transMat=NULL,type=NULL,quantileNum=NULL){
  if(dim(data)[2]!=dim(transMat)[1]){
    stop("ERROR : Cannot match data and transMat")
  }
  result=data.frame(NULL)
  if(type=="quantile"){
    for(i in 1:dim(data)[1]){
      result <- rbind(result,data.frame(data[i,]*transMat,
                                        group=paste(paste("Group",as.factor(letters[i])),rownames(iData)[i],sep=" : ")))
    }
  }else{
    for(i in 1:dim(data)[1]){
      result <- rbind(result,data.frame(data[i,]*transMat,
                                        group=paste(paste("Group",as.factor(i)),rownames(iData)[i],sep=" : ")))
    }
  }
  return(result)
}
generateCircle <- function(nLayer=NULL,extendUnit=0.5){
  if(is.null(nLayer)||nLayer<=0||nLayer>30){
    stop("Illegal layerNumber input. Recommended value will between 2 and 10.")
  }
  circles <- data.frame(
    x0 = rep(0,nLayer),
    y0 = rep(0,nLayer),
    #r = round(seq(0,myRange,1/nLayer)[2:(nLayer+1)],2)
    #r = quantile(seq(0,myRange,1/nLayer),seq(0,1,1/nLayer))[2:(nLayer+1)]
    r = quantile(seq(0+extendUnit,1+extendUnit,1/(nLayer-1)),seq(0,1,1/(nLayer-1)))
  )
  p<-ggplot() +
    ggforce::geom_circle(aes(x0 = circles$x0, y0 = circles$y0, r = circles$r), data = circles, col = "gray50")
  return(p)
}
generatePoint<-function(nPoly=NULL,nLayer=NULL,extendUnit=0.5){
  if(is.null(nPoly)||nPoly<=0||nPoly>360){
    stop("ERROR : Illegal parameter input in nPoly")
  }
  if(is.null(nLayer)||nLayer<=0||nLayer>30){
    stop("ERROR : Illegal parameter input in nLayer")
  }
  #rList = round(seq(0,1,1/nLayer)[2:(nLayer+1)],2)
  #rList = quantile(seq(0,myRange,1/nLayer),seq(0,1,1/nLayer))[2:(nLayer+1)]
  rList = quantile(seq(0+extendUnit,1+extendUnit,1/(nLayer-1)),seq(0,1,1/(nLayer-1)))
  degreeUnit<-round(360/nPoly)
  allPoint<-data.frame(NULL)
  for(r in rList){
    allPoint <- rbind(allPoint,data.frame(x=0,y=r,cos=0,sin=1))
    for(i in 1:(nPoly-1)){
      degree<-90+degreeUnit*i
      x1<-cos(degree*pi/180)*r
      y1<-sin(degree*pi/180)*r
      allPoint <- rbind(allPoint,data.frame(x=x1,y=y1,
                                            cos=cos(degree*pi/180),
                                            sin=sin(degree*pi/180)))
    }
  }
  return(allPoint)
}
getCutLine<-function(startPoint,len,perp){
  if(perp[1]==0&&perp[2]==0){
    stop("(0,0) vector cannot find a Perpendicular")
  }else if(perp[2]==0){
    if(perp[1]<0){
      return(c(startPoint[1]-len,startPoint[2]))
    }
    return(c(startPoint[1]+len,startPoint[2]))
  }else if(perp[1]==0){
    return(c(startPoint[1],startPoint[2]+len))
  }else{
    r<-perp[1]/perp[2]
  }
  #adjust <- abs(len/(perp[1]))
  adjust <- (len^2/((perp[1]^2)*(1+(1/r)^2)))^0.5
  result <-  c(perp[1]*adjust+startPoint[1] , perp[1]*adjust*(1/r)+startPoint[2])
  return(result)
}
getPerpendicular <-function(v=NULL,changeSide=TRUE){
  if(changeSide){
    if(is.vector(v)&&length(v)==2){
      if(v[1]==0&&v[2]==0){
        stop("cannot find perp. with a point")
      }else if(v[1]==0){
        return(c(-1,0))
      }else if(v[2]==0){
        return(c(0,-1))
      }else{
        return(c(-1,v[1]/v[2]))
      }
    }else{
      stop("generate perpendicular vector error, check data type and dimension.")
    }
  }else{
    if(is.vector(v)&&length(v)==2){
      if(v[1]==0&&v[2]==0){
        stop("cannot find perp. with a point")
      }else if(v[1]==0){
        return(c(1,0))
      }else if(v[2]==0){
        return(c(0,1))
      }else{
        return(c(1,-v[1]/v[2]))
      }
    }else{
      stop("generate perpendicular vector error, check data type and dimension.")
    }
  }
}
getCutDf <- function(tempDf,transMatrix,nP,indNum){
  startDf<-data.frame(x=tempDf$cos,y=tempDf$sin,group=tempDf$group)
  perp1<-lapply(1:nP,FUN=function(x) getPerpendicular(unlist(transMatrix[x,]),changeSide=T))
  perp2<-lapply(1:nP,FUN=function(x) getPerpendicular(unlist(transMatrix[x,]),changeSide=F))

  cutPos1<-sapply(1:(nP*indNum),FUN=function(x) {
    u <- x %% nP
    if(u==0){
      u<-nP
    }
    getCutLine(startDf[x,],0.025,perp1[[u]])
  })
  cutPos2<-sapply(1:(nP*indNum),FUN=function(x) {
    u <- x %% nP
    if(u==0){
      u<-nP
    }
    getCutLine(startDf[x,],0.025,perp2[[u]])
  })
  #print(as.data.frame(t(cutPos1)))
  cutDf1<-cbind(startDf,as.data.frame(t(cutPos1)))
  cutDf2<-cbind(startDf,as.data.frame(t(cutPos2)))
  colnames(cutDf1)<-c("x1","y1","obsGroup","x2","y2");colnames(cutDf2)<-c("x1","y1","obsGroup","x2","y2");
  cutDf1$x2<-as.numeric(cutDf1$x2);cutDf1$y2<-as.numeric(cutDf1$y2)
  cutDf2$x2<-as.numeric(cutDf2$x2);cutDf2$y2<-as.numeric(cutDf2$y2)
  return(list(data.frame(cutDf1),data.frame(cutDf2)))
}
build3DRect <- function(len=NULL,wid=NULL,hei=NULL,g=NULL){
  if(is.null(g)){
    g<-as.factor("g")
  }
  vertice<-get_vertice(3)
  vertice[c(3,4,7:8),2]<-vertice[c(3,4,7:8),2]+hei
  vertice[c(5:8),1]<-vertice[c(5:8),1]+wid
  vertice[c(2,4,6,8),3]<-vertice[c(2,4,6,8),3]+len
  vertice<-as.data.frame(vertice)
  vertice <- point3Dto2D(vertice)
  #rearrange
  vertice[,6:11]<-0
  vertice[c(1,2,4,3),6]<-paste0(g,1)
  vertice[c(3,4,8,7),7]<-paste0(g,2)
  vertice[c(5,6,8,7),8]<-paste0(g,3)
  vertice[c(1,2,6,5),9]<-paste0(g,4)
  vertice[c(1,3,7,5),10]<-paste0(g,5)
  vertice[c(2,4,8,6),11]<-paste0(g,6)
  d<-vertice[,c(4,5)]
  d<-rbind(d[c(1,2,4,3),],
           d[c(3,4,8,7),],
           d[c(5,6,8,7),],
           d[c(1,2,6,5),],
           d[c(1,3,7,5),],
           d[c(2,4,8,6),])
  group<-matrix(c(vertice[,6],vertice[,7],
                  vertice[,8],vertice[,9],
                  vertice[,10],vertice[,11]),ncol=1)
  group<-group[group!=0]
  d<-cbind(d,group)
  d$group <- as.factor(d$group)
  return(d)
}
reRange <- function(min,max,data){
  dist<-max-min
  return(data*dist-abs(min))
}
shift <-function(data,unit){
  unit <- 1+unit
  data$x<-round(data$x,10) ; data$y<-round(data$y,10)
  data$x <- data$x * unit
  data$y <- data$y * unit
  return(data)
}
