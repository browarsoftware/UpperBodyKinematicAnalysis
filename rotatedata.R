generateroty <- function(angley)
{
  ry <- matrix(nrow=3, ncol=3)
  recangley <- angley * 2 * pi / 360
  ry[1,1] <- cos(recangley)
  ry[1,2] <- 0
  ry[1,3] <- sin(recangley)
  
  ry[2,1] <- 0
  ry[2,2] <- 1
  ry[2,3] <- 0
  
  ry[3,1] <- -sin(recangley)
  ry[3,2] <- 0
  ry[3,3] <- cos(recangley)
  
  return (ry)  
}

require('subplex')

#mydata <- read.csv("e:\\Publikacje\\artyku³y kine\\sample3.csv")
#referencedata <- read.csv("e:\\Publikacje\\artyku³y kine\\sample7.csv")

#mydata <- read.csv("f:\\mocap\\karate\\2016-12-22 ShoriunRiu MP\\mae_geri_right\\segmented\\sample3.csv")
#referencedata <- read.csv("f:\\mocap\\karate\\2016-12-22 ShoriunRiu MP\\mae_geri_right\\segmented\\sample7.csv")

rotatedata <- function(mydata, referencedata, referencebodypartnamex, referencebodypartnamez,
                       hipsbodypartnamex, hipsbodypartnamez)
{
  columnsToFind <- colnames(mydata)[grep("[[:alnum:]]+\\.Dx", colnames(mydata), ignore.case = TRUE)]
  columnsToFind <- c(columnsToFind, colnames(mydata)[grep("[[:alnum:]]+\\.Dy", colnames(mydata), ignore.case = TRUE)])
  columnsToFind <- c(columnsToFind, colnames(mydata)[grep("[[:alnum:]]+\\.Dz", colnames(mydata), ignore.case = TRUE)])
  
  columnsToFindNoEnds <- list()
  for (a in 1:length(columnsToFind))
  {
    columnsToFindNoEnds[a] <- substring(text = columnsToFind[a], first = 1, last = nchar(columnsToFind[a]) - 3)
  }
  columnsToFindNoEnds <- unique(unlist(columnsToFindNoEnds))
  
  allResults <- list()
  allResultsColnames <- list()
  Time <- (1:length(mydata[,hipsbodypartnamex])) / 100
  
  df <- data.frame(Time = Time)
  
  
  #v1 <- c(mydata$LeftFoot.Dx[1] - mydata$Hips.Dx[1], 0, mydata$LeftFoot.Dz[1] - mydata$Hips.Dz[1])
  #v2 <- c(referencedata$LeftFoot.Dx[1] - referencedata$Hips.Dx[1], 0, referencedata$LeftFoot.Dz[1] - referencedata$Hips.Dz[1])
  
  v1 <- c(mydata[1,referencebodypartnamex] - mydata[1,hipsbodypartnamex], 
          0, 
          mydata[1,referencebodypartnamez] - mydata[1,hipsbodypartnamez])
  v2 <- c(referencedata[1,referencebodypartnamex] - referencedata[1,hipsbodypartnamex], 
          0, 
          referencedata[1,referencebodypartnamez] - referencedata[1,hipsbodypartnamez])
  
  
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  optimizeangle <- function(x)
  {
    my <- generateroty(x)
    v11 <- my %*% v1
    #v11 <- rotres[1,]
    return (euc.dist(v11, v2))
  }
  response <- subplex(par=c(0),fn=optimizeangle)
  response$par
  
  #euc.dist(v1, v2)
  #euc.dist(my %*% matrix(v1), v2)
  
  #my <- generateroty(response$par)
  #v1
  #my %*% matrix(v1)
  #v2
  
  #vec1 <- rotres[1,]
  
  #v1
  #v2
  
  #columnsToFind
  
  xx <- grep("[[:alnum:]]+\\.Dx", colnames(mydata), ignore.case = TRUE)[1]
  yy <- grep("[[:alnum:]]+\\.Dy", colnames(mydata), ignore.case = TRUE)[1]
  zz <- grep("[[:alnum:]]+\\.Dz", colnames(mydata), ignore.case = TRUE)[1]
  
  tt <- colnames(mydata[xx])
  xx <- substring(tt,nchar(tt)-2,nchar(tt))
  tt <- colnames(mydata[yy])
  yy <- substring(tt,nchar(tt)-2,nchar(tt))
  tt <- colnames(mydata[zz])
  zz <- substring(tt,nchar(tt)-2,nchar(tt))
  
  for (a in 1:length(columnsToFindNoEnds))
  {
    #a = 2
    xx1 <- paste(columnsToFindNoEnds[a] , xx, sep = "")
    yy1 <- paste(columnsToFindNoEnds[a] , yy, sep = "")
    zz1 <- paste(columnsToFindNoEnds[a] , zz, sep = "")
    
    xxx <- list()
    yyy <- list()
    zzz <- list()
    
    for (b in 1:length(Time))
    {
      vec1 <- c(mydata[b,xx1], mydata[b,yy1], mydata[b,zz1])
      my <- generateroty(response$par)
      #rotres <- my * vec1
      vv <- as.vector(my %*% matrix(vec1))
      
      xxx[[b]] <- vv[1]
      yyy[[b]] <- vv[2]
      zzz[[b]] <- vv[3]
    }
    df[xx1] <- unlist(xxx)
    df[yy1] <- unlist(yyy)
    df[zz1] <- unlist(zzz)
  }
  return (df)
}