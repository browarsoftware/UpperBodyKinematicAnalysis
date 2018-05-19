euc.dist1d <- function(x1, x2) abs(x1[1] - x2[1])
library(compiler)
euc.dist1dCmp <- cmpfun(euc.dist1d)

vec.length <- function(x) sqrt(sum(x * x))

euc.dist1dangle <- function(x1, x2)
{
  angle <- abs(x1[1] - x2[1])
  if (angle > pi)
    angle <- (2 * pi) - angle
  return (angle)
}
library(compiler)
euc.dist1dangleCmp <- cmpfun(euc.dist1dangle)

calculatenormalizeddistance <- function(path1, path2, signal1, signal2, FUN)
{
  distanceHelper <- 0
  for (a in 1:length(path1))
  {
    distanceHelper <- distanceHelper + FUN(signal1[[path1[a]]], signal2[[path2[a]]])
  }
  return (distanceHelper / (length(signal1) + length(signal2)))
}

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
library(compiler)
euc.distCmp <- cmpfun(euc.dist)

angle.between <- function(a,b) acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
angle.betweenCmp <- cmpfun(angle.between)

vector.dot <- function(a,b) sum(a*b)
vector.dotCmp <- cmpfun(vector.dot)

vector.norm <- function(x) sqrt(sum(x^2))
vector.norm <- cmpfun(vector.norm)

calc.angle <- function(a,b,vn) 
{
  angle <- angle.betweenCmp(a,b) * 180 / pi
  if (sign(vector.dot(vn, vector.cross(a, b))) > 0)
    return (angle)
  else
    return ((360) - angle)
}


vector.cross <- function(a, b) {
  if(length(a)!=3 || length(b)!=3){
    stop("Cross product is only defined for 3D vectors.");
  }
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (a[i1]*b[i2] - a[i2]*b[i1])
}

ArgMin3 <- function(a,b,c)
{
  if (a<b)
  {
    if (a<c)
    {
      return (0)
    }
    else
    {
      return (2)
    }
  }
  else
  {
    if (b<c)
    {
      return (1)
    }
    else
    {
      return (2);
    }
  }
  return (0)
}

library(compiler)
ArgMin3Cmp <- cmpfun(ArgMin3)


distanceTo <- function(a,b)
{
  dist <- (a-b)*(a-b)
  return (dist)
}



myDTW <- function(FUN,averageS,sequence)
{
  #averageS <- rs
  #sequence <- s1
  
  tupleAssociation <- list();
  #for (t in 1:length(averageS))
  #  tupleAssociation[[t]] <- list();
  for (t in 1:length(averageS))
    tupleAssociation[[t]] <- data.frame(v1 = numeric(), 
                                        v2 = numeric(),
                                        v3 = numeric(),
                                        v4 = numeric(),
                                        stringsAsFactors = FALSE);
  
  #for t=1:size(averageS,2)
  #tupleAssociation{t}=[];
  #end
  
  sl <- length(averageS) * length(sequence)
  seq1 <- rep(0, sl)
  #mat1 <- matrix(seq1, length(sequences[[1]]))
  
  costMatrix <- matrix(seq1, length(averageS))
  pathMatrix <- matrix(seq1, length(averageS))
  
  
  
  #costMatrix[1,1] <- quat_similarityCmp(unlist(averageS[1]),unlist(sequence[1]))
  costMatrix[1,1] <- FUN(unlist(averageS[1]),unlist(sequence[1]))
  
  
  pathMatrix[1,1] <- -1;
  
  
  for (i in 2:length(averageS))
  {
    #costMatrix[i,1] <- costMatrix[i-1,1] + quat_similarityCmp(unlist(averageS[i]),unlist(sequence[1]));
    costMatrix[i,1] <- costMatrix[i-1,1] + FUN(unlist(averageS[i]),unlist(sequence[1]));
    
    pathMatrix[i,1] <- 2;
  }
  
  for (j in 2:length(sequence))
  {
    #costMatrix[1,j] <- costMatrix[1,j-1] + quat_similarityCmp(unlist(sequence[j]),unlist(averageS[1]));
    costMatrix[1,j] <- costMatrix[1,j-1] + FUN(unlist(sequence[j]),unlist(averageS[1]));
    pathMatrix[1,j] <- 1;
  }
  
  for (i in 2:length(averageS))
  {
    for (j in 2:length(sequence))
    {
      #indiceRes <- ArgMin3(costMatrix[i-1,j-1],costMatrix[i,j-1],costMatrix[i-1,j]);
      indiceRes <- ArgMin3Cmp(costMatrix[i-1,j-1],costMatrix[i,j-1],costMatrix[i-1,j]);
      pathMatrix[i,j] <- indiceRes;
      
      if (indiceRes==0)
      {
        res <- costMatrix[i-1,j-1];
      }
      else if (indiceRes==1)
      {
        res <- costMatrix[i,j-1]
      }
      else if (indiceRes==2)
      {
        res <- costMatrix[i-1,j]
      }
      #costMatrix[i,j] <- res + distanceTo(averageS[i],sequence[j])
      #costMatrix[i,j] <- res + quat_similarity(unlist(averageS[i]),unlist(sequence[j]))
      #costMatrix[i,j] <- res + quat_similarityCmp(unlist(averageS[i]),unlist(sequence[j]))
      costMatrix[i,j] <- res + FUN(unlist(averageS[i]),unlist(sequence[j]))
      
    }
  }
  
  i <- length(averageS)
  j <- length(sequence)
  
  distance <- 0
  
  path1 <- list()
  path2 <- list()
  
  a <- 1
  while(TRUE)
  {
    path1[[a]] <- i
    path2[[a]] <- j
    a <- a + 1
    
    ttt <- tupleAssociation[[i]]
    nr <- nrow(ttt) + 1
    ttt[nr,1:4] <- unlist(sequence[j])[1:4]
    
    tupleAssociation[[i]] <- ttt
    #print(costMatrix[i,j])
    distance <- costMatrix[i,j] + distance
    #print(costMatrix[i,j])
    if (pathMatrix[i,j]==0)
    {
      i=i-1;
      j=j-1;
      #print(paste('0)', i))
    } else if (pathMatrix[i,j]==1)
    {
      j=j-1;
      #print(paste('1)', i))
    } else if (pathMatrix[i,j]==2)
    {
      i=i-1;          
      #print(paste('2)', i))
    } else
    {
      break
    }
  }
  
  path1 <- rev(unlist(path1))
  path2 <- rev(unlist(path2))
  
  #plot(path1, path2)
  #end
  
  #return(distance / (length(averageS) + length(sequence)))
  #return (costMatrix[length(averageS), length(sequence)] / (length(averageS) + length(sequence)))
  normalized_distance = costMatrix[length(averageS), length(sequence)] / (length(averageS) + length(sequence))
  newList <- list("path1" = path1, "path2" = path2, 'normalized_distance' = normalized_distance)
  return (newList)
}
myDTWCmp <- cmpfun(myDTW)


generateFeatures <- function(dataToCalculate)
{
  dx_vector <- names(dataToCalculate)[grepl(".Dx",names(dataToCalculate))]

  mysubstring <- function(x)
  {
    return (substring(x, first = 1, last = nchar(x) - 3))
  }
  
  featuresNames <- unlist(lapply(dx_vector, mysubstring))
  
  listhelper <- list()
  df <- data.frame(a=1:length(dataToCalculate[,1]))
  
  
  
  
  dataToCalculate[paste(featuresNames[1],".Dx",sep = "")]
  
  for (a in 1:length(featuresNames))
  {
    listhelper <- list()
    for (b in 1:length(dataToCalculate[,paste(featuresNames[a],".Dx",sep = "")]))
      listhelper[[b]] <- c(dataToCalculate[b,paste(featuresNames[a],".Dx",sep = "")], 
                           dataToCalculate[b,paste(featuresNames[a],".Dy",sep = "")], 
                           dataToCalculate[b,paste(featuresNames[a],".Dz",sep = "")])
    df[[featuresNames[a]]] <- listhelper
  }
  
  ############################################
  #vectors
  dataRightKnee <- list()
  dataLeftKnee <- list()

  dataListRightThighX <- list()
  dataListRightThighY <- list()
  dataListRightThighZ <- list()
  
  dataListLeftThighX <- list()
  dataListLeftThighY <- list()
  dataListLeftThighZ <- list()
  

  dataListHipsX <- list()
  dataListHipsY <- list()
  dataListHipsZ <- list()
  
  dataRightElbow <- list()
  dataLeftElbow <- list()
  
  dataListRightShoulderX <- list()
  dataListRightShoulderY <- list()
  dataListRightShoulderZ <- list()
  
  dataListLeftShoulderX <- list()
  dataListLeftShoulderY <- list()
  dataListLeftShoulderZ <- list()
  
  
  dataListLeftArmX <- list()
  dataListLeftArmY <- list()
  dataListLeftArmZ <- list()
  
  dataListRightArmX <- list()
  dataListRightArmY <- list()
  dataListRightArmZ <- list()
  
  for (fposition in 1:length(dataToCalculate$Hips.Dx))
  {
    ################################
    #Legs
    
    #Knees
    dataRightKnee[[fposition]] <- angle.between(c(dataToCalculate$RightLeg.Dx[fposition] - dataToCalculate$RightThigh.Dx[fposition],
                                                  dataToCalculate$RightLeg.Dy[fposition] - dataToCalculate$RightThigh.Dy[fposition],
                                                  dataToCalculate$RightLeg.Dz[fposition] - dataToCalculate$RightThigh.Dz[fposition]),
                                                c(dataToCalculate$RightLeg.Dx[fposition] - dataToCalculate$RightFoot.Dx[fposition], 
                                                  dataToCalculate$RightLeg.Dy[fposition] - dataToCalculate$RightFoot.Dy[fposition],
                                                  dataToCalculate$RightLeg.Dz[fposition] - dataToCalculate$RightFoot.Dz[fposition]))
                                                
    
    dataLeftKnee[[fposition]] <- angle.between(c(dataToCalculate$LeftLeg.Dx[fposition] - dataToCalculate$LeftThigh.Dx[fposition],
                                                  dataToCalculate$LeftLeg.Dy[fposition] - dataToCalculate$LeftThigh.Dy[fposition],
                                                  dataToCalculate$LeftLeg.Dz[fposition] - dataToCalculate$LeftThigh.Dz[fposition]),
                                                c(dataToCalculate$LeftLeg.Dx[fposition] - dataToCalculate$LeftFoot.Dx[fposition], 
                                                  dataToCalculate$LeftLeg.Dy[fposition] - dataToCalculate$LeftFoot.Dy[fposition],
                                                  dataToCalculate$LeftLeg.Dz[fposition] - dataToCalculate$LeftFoot.Dz[fposition]))
    
    #TIGHT JOINTS    
    xt <- c(dataToCalculate$RightThigh.Dx[fposition], dataToCalculate$RightThigh.Dy[fposition], dataToCalculate$RightThigh.Dz[fposition]) - c(dataToCalculate$LeftThigh.Dx[fposition], dataToCalculate$LeftThigh.Dy[fposition], dataToCalculate$LeftThigh.Dz[fposition])
    xxt <- xt / vector.norm(xt)
    yt <- c(0,1,0)
    zt <- vector.cross(xxt, yt)
    zzt <-  zt / vector.norm(zt)
    yyt <- vector.cross(xt, zzt)
    yyt <- yyt / vector.norm(yyt)
    
    dataListRightLeg <- c(dataToCalculate$RightThigh.Dx[fposition] - dataToCalculate$RightLeg.Dx[fposition], 
                         dataToCalculate$RightThigh.Dy[fposition] - dataToCalculate$RightLeg.Dy[fposition],
                         dataToCalculate$RightThigh.Dz[fposition] - dataToCalculate$RightLeg.Dz[fposition])

    dataListRightThighX[[fposition]] <- angle.between(yyt,dataListRightLeg)
    dataListRightThighY[[fposition]] <- angle.between(zzt,dataListRightLeg)
    dataListRightThighZ[[fposition]] <- angle.between(xxt,dataListRightLeg)
    

    dataListLeftLeg <- c(dataToCalculate$LeftThigh.Dx[fposition] - dataToCalculate$LeftLeg.Dx[fposition], 
                                      dataToCalculate$LeftThigh.Dy[fposition] - dataToCalculate$LeftLeg.Dy[fposition],
                                      dataToCalculate$LeftThigh.Dz[fposition] - dataToCalculate$LeftLeg.Dz[fposition])
    
    dataListLeftThighX[[fposition]] <- angle.between(yyt,dataListLeftLeg)
    dataListLeftThighY[[fposition]] <- angle.between(zzt,dataListLeftLeg)
    dataListLeftThighZ[[fposition]] <- angle.between(xxt,dataListLeftLeg)
    
    #HIPS
    
    hipschest <- c(dataToCalculate$Chest.Dx[fposition] - dataToCalculate$Hips.Dx[fposition], 
                   dataToCalculate$Chest.Dy[fposition] - dataToCalculate$Hips.Dy[fposition],
                   dataToCalculate$Chest.Dz[fposition] - dataToCalculate$Hips.Dz[fposition])
    
    dataListHipsX[[fposition]] <- angle.between(c(0,1,0), hipschest)
    dataListHipsY[[fposition]] <- angle.between(c(0,0,1),xt)
    dataListHipsZ[[fposition]] <- angle.between(c(1,0,0), hipschest)
    
    ####################
    #UPPER BODY
    
    #Elbows
    
    dataRightElbow[[fposition]] <- angle.between(c(dataToCalculate$RightArm.Dx[fposition] - dataToCalculate$RightShoulder.Dx[fposition],
                                                  dataToCalculate$RightArm.Dy[fposition] - dataToCalculate$RightShoulder.Dy[fposition],
                                                  dataToCalculate$RightArm.Dz[fposition] - dataToCalculate$RightShoulder.Dz[fposition]),
                                                c(dataToCalculate$RightArm.Dx[fposition] - dataToCalculate$RightForearm.Dx[fposition], 
                                                  dataToCalculate$RightArm.Dy[fposition] - dataToCalculate$RightForearm.Dy[fposition],
                                                  dataToCalculate$RightArm.Dz[fposition] - dataToCalculate$RightForearm.Dz[fposition]))
    
    
    dataLeftElbow[[fposition]] <- angle.between(c(dataToCalculate$LeftArm.Dx[fposition] - dataToCalculate$LeftShoulder.Dx[fposition],
                                                 dataToCalculate$LeftArm.Dy[fposition] - dataToCalculate$LeftShoulder.Dy[fposition],
                                                 dataToCalculate$LeftArm.Dz[fposition] - dataToCalculate$LeftShoulder.Dz[fposition]),
                                               c(dataToCalculate$LeftArm.Dx[fposition] - dataToCalculate$LeftForearm.Dx[fposition], 
                                                 dataToCalculate$LeftArm.Dy[fposition] - dataToCalculate$LeftForearm.Dy[fposition],
                                                 dataToCalculate$LeftArm.Dz[fposition] - dataToCalculate$LeftForearm.Dz[fposition]))
    
    #Shoulders
    
    xt <- c(dataToCalculate$RightShoulder.Dx[fposition], dataToCalculate$RightShoulder.Dy[fposition], dataToCalculate$RightShoulder.Dz[fposition]) - c(dataToCalculate$LeftShoulder.Dx[fposition], dataToCalculate$LeftShoulder.Dy[fposition], dataToCalculate$LeftShoulder.Dz[fposition])
    xxt <- xt / vector.norm(xt)
    yt <- c(0,1,0)
    zt <- vector.cross(xxt, yt)
    zzt <-  zt / vector.norm(zt)
    yyt <- vector.cross(xt, zzt)
    yyt <- yyt / vector.norm(yyt)
    
    
    dataListLeftShoulder <- c(dataToCalculate$LeftShoulder.Dx[fposition] - dataToCalculate$LeftArm.Dx[fposition], 
                              dataToCalculate$LeftShoulder.Dy[fposition] - dataToCalculate$LeftArm.Dy[fposition],
                              dataToCalculate$LeftShoulder.Dz[fposition] - dataToCalculate$LeftArm.Dz[fposition])
    
    #dataListLeftShoulderX[[fposition]] <- angle.between(yyt,dataListLeftShoulder)
    #dataListLeftShoulderY[[fposition]] <- angle.between(zzt,dataListLeftShoulder)
    #dataListLeftShoulderZ[[fposition]] <- angle.between(xxt,dataListLeftShoulder)
    
    v <- dataListLeftShoulder
    n <- vector.cross(xxt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListLeftShoulderY[[fposition]] <- angle.between(proj3, xxt)
    
    n <- vector.cross(yyt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListLeftShoulderZ[[fposition]] <- angle.between(proj3, yyt)
    
    
    n <- vector.cross(zzt, xxt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListLeftShoulderX[[fposition]] <- angle.between(proj3, zzt)
    
    dataListRightShoulder <- c(dataToCalculate$RightShoulder.Dx[fposition] - dataToCalculate$RightArm.Dx[fposition], 
                              dataToCalculate$RightShoulder.Dy[fposition] - dataToCalculate$RightArm.Dy[fposition],
                              dataToCalculate$RightShoulder.Dz[fposition] - dataToCalculate$RightArm.Dz[fposition])
    
    #dataListRightShoulderX[[fposition]] <- angle.between(yyt,dataListRightShoulder)
    #dataListRightShoulderY[[fposition]] <- angle.between(zzt,dataListRightShoulder)
    #dataListRightShoulderZ[[fposition]] <- angle.between(xxt,dataListRightShoulder)
    
    v <- dataListRightShoulder
    n <- vector.cross(xxt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListRightShoulderY[[fposition]] <- angle.between(proj3, xxt)
    
    n <- vector.cross(yyt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListRightShoulderZ[[fposition]] <- angle.between(proj3, yyt)
    
    
    n <- vector.cross(zzt, xxt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListRightShoulderX[[fposition]] <- angle.between(proj3, zzt)
    
    #############################

    dataListLeftArm <- c(dataToCalculate$LeftArm.Dx[fposition] - dataToCalculate$LeftForearm.Dx[fposition], 
                              dataToCalculate$LeftArm.Dy[fposition] - dataToCalculate$LeftForearm.Dy[fposition],
                              dataToCalculate$LeftArm.Dz[fposition] - dataToCalculate$LeftForearm.Dz[fposition])
    
    #dataListLeftArmX[[fposition]] <- angle.between(yyt,dataListLeftArm)
    #dataListLeftArmY[[fposition]] <- angle.between(zzt,dataListLeftArm)
    #dataListLeftArmZ[[fposition]] <- angle.between(xxt,dataListLeftArm)
    
    v <- dataListLeftArm
    n <- vector.cross(xxt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListLeftArmY[[fposition]] <- angle.between(proj3, xxt)
    
    n <- vector.cross(yyt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListLeftArmZ[[fposition]] <- angle.between(proj3, yyt)
    
    
    n <- vector.cross(zzt, xxt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListLeftArmX[[fposition]] <- angle.between(proj3, zzt)
    
    
    dataListRightArm <- c(dataToCalculate$RightArm.Dx[fposition] - dataToCalculate$RightForearm.Dx[fposition], 
                               dataToCalculate$RightArm.Dy[fposition] - dataToCalculate$RightForearm.Dy[fposition],
                               dataToCalculate$RightArm.Dz[fposition] - dataToCalculate$RightForearm.Dz[fposition])
    
    #dataListRightArmX[[fposition]] <- angle.between(yyt,dataListRightArm)
    #dataListRightArmY[[fposition]] <- angle.between(zzt,dataListRightArm)
    #dataListRightArmZ[[fposition]] <- angle.between(xxt,dataListRightArm)
    
    
    
    v <- dataListRightArm
    n <- vector.cross(xxt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListRightArmY[[fposition]] <- angle.between(proj3, xxt)
    
    n <- vector.cross(yyt, zzt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListRightArmZ[[fposition]] <- angle.between(proj3, yyt)
    
    
    n <- vector.cross(zzt, xxt)
    a <- vector.cross(n, vector.cross(v, n))
    #proj <- v - ((vector.dot(a, v) / vector.dot(a, a)) * a)
    proj3 <- v - (n * (vector.dot(v, n)))
    dataListRightArmX[[fposition]] <- angle.between(proj3, zzt)
  }
  
  df$ListRightKnee <- dataRightKnee
  df$ListLeftKnee <- dataLeftKnee
  
  df$ListLeftThighX <- dataListLeftThighX
  df$ListLeftThighY <- dataListLeftThighY
  df$ListLeftThighZ <- dataListLeftThighZ
  
  df$ListRightThighX <- dataListRightThighX
  df$ListRightThighY <- dataListRightThighY
  df$ListRightThighZ <- dataListRightThighZ
  
  df$ListHipsX <- dataListHipsX
  df$ListHipsY <- dataListHipsY
  df$ListHipsZ <- dataListHipsZ
  
  df$ListRightElbow <- dataRightElbow
  df$ListLeftElbow <- dataLeftElbow
  
  df$ListLeftShoulderX <- dataListLeftShoulderX
  df$ListLeftShoulderY <- dataListLeftShoulderY
  df$ListLeftShoulderZ <- dataListLeftShoulderZ
  
  df$ListRightShoulderX <- dataListRightShoulderX
  df$ListRightShoulderY <- dataListRightShoulderY
  df$ListRightShoulderZ <- dataListRightShoulderZ
  
  
  df$ListLeftArmX <- dataListLeftArmX
  df$ListLeftArmY <- dataListLeftArmY
  df$ListLeftArmZ <- dataListLeftArmZ
  
  df$ListRightArmX <- dataListRightArmX
  df$ListRightArmY <- dataListRightArmY
  df$ListRightArmZ <- dataListRightArmZ
  
  
  df$a <- NULL
  

  return(df)
}



plotsmoothingresults <- function(smoothingresults, plottitle, plotifnoextreams = TRUE, plotsmoothed = FALSE, ylab = "Distance [cm]", legenPosition = "topright")
{
  #smoothingresults <- footddf
  #plottitle <- "footddf"
  #plotifnoextreams <- TRUE
  #plotsmoothed <- TRUE
  
  #plot smoothed data
  if (plotsmoothed)
  {
    plot(smoothingresults$smoothdata , col = "black")
    title(main = plottitle)
    lines(smoothingresults$smoothdata , col = "black")
    
    for (a in 1:length(smoothingresults$extremumbool))
    {
      if (smoothingresults$extremumbool[a])
      {
        points(a, smoothingresults$smoothdata[a], col = "red",pch = 4)
        idhelp <- a - 1
        end <- FALSE
        #footddf$extremumtreshold
        while (smoothingresults$derivative[idhelp] > 0 && !end)
        {
          points(idhelp, smoothingresults$smoothdata[idhelp], col = "blue",pch = 4)
          if (idhelp > 1)
            idhelp <- idhelp - 1
          else
            end <- TRUE
          
        }
      }
    }
  }
  if (plotifnoextreams || length(smoothingresults$resultsList) > 0)
  {
    #plot(smoothingresults$data , col = "black")
    plot(smoothingresults$data, xlab = "Time [10^-100 s]", ylab = ylab, col = 'black', type='l',
	ylim = c(min(smoothingresults$data), max(smoothingresults$data) * 1.5))
    
    
    title(main = plottitle)
    #lines(smoothingresults$data , col = "black")
    if (length(smoothingresults$resultsList) > 0)
      for (a in 1:length(smoothingresults$resultsList))
      {
        vec <- smoothingresults$resultsList[[a]]
        for (b in 1:length(vec))
        {
          points(vec[b], smoothingresults$data[vec[b]], col = "red",pch = 4, lwd=3)
          idhelp <- vec[b] - 1
          end <- FALSE
          #footddf$extremumtreshold
          while (smoothingresults$derivative[idhelp] > 0 && !end)
          {
            points(idhelp, smoothingresults$data[idhelp], col = "blue",pch = 4)
            if (idhelp > 1)
              idhelp <- idhelp - 1
            else
              end <- TRUE
            
          }
        }
      }
    legend(x= legenPosition, y=max(smoothingresults$data) * 1.5, legend=c("Original", "Maxima over treshold", "ROI"), col=c("black", 'red', "blue"), lty=c(1,1,1), cex=0.8)
  }
}


rglplotanalyzedata <- function(refdatakinematic, inputdataalignmentkinematic, xx1, xx2, path1, path2, resultdata, whattodraw = "LeftFoot")
{
  #xx1 <- refdatakinematicf$dataRightKnee
  #xx2 <- inputdataalignmentkinematicf$dataRightKnee
  #path1 <- footddf$path1
  #path2 <- footddf$path2
  #resultdata <- kneeadf
  #whattodraw <- "LeftLeg"
  library("rgl")
  rgl.open() # Open a new RGL device
  #rgl.bg(col="white")
  #pointscolors <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  pointscolors <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  alpha <- 0.1
  
  idref <- -1
  idinput <- -1
  
  #for (a in seq(1, length(path1), 20))
  #a = 1
  #for (a in length(path1):length(path1))
  for (a in 1:length(path1))
  {
    if (!(path2[a] %in% resultdata$extremumid))
    {
      if (idref != path1[a])
      {
        idref = path1[a]
        #renderactor(refdatakinematic, idref, pointscolors, "green", "green", a / (4 * length(path1)), showspheres = FALSE)
        renderactor(refdatakinematic, idref, pointscolors, "green", "green", 16 / length(path1), showspheres = FALSE,linewidth =  2)
      }
      if (idinput != path2[a])
      {
        idinput = path2[a]
        #renderactor(inputdataalignmentkinematic, idinput, pointscolors, "red", "red", a / (4 * length(path2)), showspheres = FALSE)
        renderactor(inputdataalignmentkinematic, idinput, pointscolors, "red", "red", 16 / length(path2), showspheres = FALSE, 2)
      }
      
      if (path2[a] %in% resultdata$extremumid)
      {
        lines3d(c(refdatakinematic[path1[a],paste(whattodraw, ".Dx",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dx",sep = "")]),
                c(refdatakinematic[path1[a],paste(whattodraw, ".Dy",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dy",sep = "")]),
                c(refdatakinematic[path1[a],paste(whattodraw, ".Dz",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dz",sep = "")])
                , color = "yellow", alpha = 1, lwd=5)
        #rgl.texts((xx2[[path2[a]]][1] + xx1[[path1[a]]][1]) / 2, 
        #          (xx2[[path2[a]]][2] + xx1[[path1[a]]][2]) / 2, 
        #          (xx2[[path2[a]]][3] + xx1[[path1[a]]][3]) / 2, 
        #          a, lwd = 5)
      }
      else
      {
        lines3d(c(refdatakinematic[path1[a],paste(whattodraw, ".Dx",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dx",sep = "")]),
                c(refdatakinematic[path1[a],paste(whattodraw, ".Dy",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dy",sep = "")]),
                c(refdatakinematic[path1[a],paste(whattodraw, ".Dz",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dz",sep = "")])
                , color = "blue", alpha = 1, lwd=1)
        #rgl.texts((xx2[[path2[a]]][1] + xx1[[path1[a]]][1]) / 2, 
        #          (xx2[[path2[a]]][2] + xx1[[path1[a]]][2]) / 2, 
        #          (xx2[[path2[a]]][3] + xx1[[path1[a]]][3]) / 2, 
        #        a, lwd = 1)
      }
    }
  }
  for (b in resultdata$extremumid)
  {
    for (a in 1:length(path2))
    {
      if (path2[a] == b)
      {
        renderactor(refdatakinematic, path1[a], pointscolors, "green", "green", 1, showspheres = FALSE, 5)
        renderactor(inputdataalignmentkinematic, path2[a], pointscolors, "red", "red", 1, showspheres = FALSE, 5)
        
        lines3d(c(refdatakinematic[path1[a],paste(whattodraw, ".Dx",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dx",sep = "")]),
                c(refdatakinematic[path1[a],paste(whattodraw, ".Dy",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dy",sep = "")]),
                c(refdatakinematic[path1[a],paste(whattodraw, ".Dz",sep = "")],inputdataalignmentkinematic[path2[a],paste(whattodraw, ".Dz",sep = "")])
                , color = "yellow", alpha = 1.0, lwd=50)
        #rgl.texts((xx2[[path2[a]]][1] + xx1[[path1[a]]][1]) / 2, 
        #          (xx2[[path2[a]]][2] + xx1[[path1[a]]][2]) / 2, 
        #          (xx2[[path2[a]]][3] + xx1[[path1[a]]][3]) / 2, 
        #          a, lwd = 5)
      }
    }
  }
  
  borderx <- c(-100, -100, 100, 100, -100, -100, 100, 100) * 1
  bordery <- c(200, -20, 200, -20, 200, -20, 200, -20) * 1
  borderz <- c(-100, -100, -100, -100, 100, 100, 100, 100) * 1
  
  spheres3d(borderx, bordery, borderz, col = "red", r=0.01)
  
  
  y <- 0#min(refdata[1,featurename])
  planes3d(0,-1,0,y, col = 'orange', alpha = 0.3)
  #axes3d(col="black", alpha = 1, lwd = 1)
  axes3d(col="white", alpha = 1, lwd = 1)  
  
  
  
	fposition <- 1
	xt <- c(refdatakinematic$RightShoulder.Dx[fposition], refdatakinematic$RightShoulder.Dy[fposition], refdatakinematic$RightShoulder.Dz[fposition]) - c(refdatakinematic$LeftShoulder.Dx[fposition], refdatakinematic$LeftShoulder.Dy[fposition], refdatakinematic$LeftShoulder.Dz[fposition])
	xxt <- xt / vector.norm(xt)
	yt <- c(0,1,0)
	zt <- vector.cross(xxt, yt)
	zzt <-  zt / vector.norm(zt)
	yyt <- vector.cross(xt, zzt)
	yyt <- yyt / vector.norm(yyt)

	xh <- refdatakinematic$Hips.Dx[fposition] + (xxt[1] * 50)
	yh <- refdatakinematic$Hips.Dy[fposition]
	zh <- refdatakinematic$Hips.Dz[fposition] + (xxt[3] * 50)
	
	lines3d(c(0 + xh, xxt[1] * 50 + xh),
			c(0 + yh, xxt[2] * 50 + yh),
			c(0 + zh, xxt[3] * 50 + zh),alpha = 1, lwd=15, color = "orange")
	lines3d(c(0 + xh, yyt[1] * 50 + xh),
			c(0 + yh, yyt[2] * 50 + yh),
			c(0 + zh, yyt[3] * 50 + zh),alpha = 1, lwd=15, color = "white")
	lines3d(c(0 + xh, zzt[1] * 50 + xh),
			c(0 + yh, zzt[2] * 50 + yh),
			c(0 + zh, zzt[3] * 50 + zh),alpha = 1, lwd=15, color = "black")
}

scatterplotanalyzedta <- function(x1, y1, z1, x2, y2, z2, path1, path2)
{
  
  require('scatterplot3d')
  color <- rep("green", length(x1))
  
  #s3d <- scatterplot3d(x1, z1, y1, color, type='l', pch=20, main="DTW signal mapping of RightThigh roation", angle=90, xlim=c(-1, 1), ylim=c(-1, 0), zlim=c(-1, -0.3),
  #s3d <- scatterplot3d(x1, z1, y1, color, type='l', pch=20, main="DTW signal mapping of Hips roation", angle=45, xlim=c(-1, 1), ylim=c(-1, 1), zlim=c(-1, 1),
  #                     s3d <- scatterplot3d(x1, z1, y1, color, type='l', pch=20, main=paste("DTW signal mapping of", "Hips", "roation"), angle=90, 
  
  #at the begining it was angle=45
  s3d <- scatterplot3d(x1, z1, y1, color, type='l', pch=20, main=paste("DTW signal mapping of", "Hips", "rotation"),# angle=90,
                       #xlim=c(-0.5, 0.6), ylim=c(-1, 0.1), zlim=c(-.1, 0.9),
                       xlab = 'X',ylab = 'Z',zlab = 'Y')
  s3d$points3d(x1, z1, y1, col='green', type="p", pch=20)
  s3d$points3d(x2, z2, y2, col='blue', type="l", pch=20)
  s3d$points3d(x2, z2, y2, col='blue', type="p", pch=20)
  text(s3d$xyz.convert(x  = x1, y = z1, z = y1), labels = 1:length(x1), col='green')
  text(s3d$xyz.convert(x  = x2, y = z2, z = y2), labels = 1:length(x2), col='blue')
  
  for (a in 1:length(path1))
  {
    x3 <- c(x1[path1[a]], x2[path2[a]])
    y3 <- c(y1[path1[a]], y2[path2[a]])
    z3 <- c(z1[path1[a]], z2[path2[a]])
    s3d$points3d(x3, z3, y3, col='red', type="l", pch=20)
  }
  
}

smothdata2 <- function(dd, smoothSize = 0.2, extremumtreshold = -1, smoothSizeHelper = 0.1)
{
  require(smoother)
  #smoothSize <- 0.1
  #smoothSize <- 0.2
  #smoothSizeHelper <- 0.1
  
  options('smoother.window' = smoothSize)
  ddsmooth <- smth.gaussian(dd, tails = TRUE)
  minvalue = min(ddsmooth)
  maxvalue = max(ddsmooth)
  diffvalue <- maxvalue - minvalue
  
  #findAllMinimums(x)
  
  derivative <- rep(0, length(ddsmooth))
  for (a in 2:(length(ddsmooth)-1))
  {
    derivative[a] <- (ddsmooth[a+1] - ddsmooth[a-1]) / 2
  }
  
  extremum <- rep(0, length(ddsmooth))
  
  
  for (a in 2:(length(derivative)-1))
  {
    if (derivative[a] < 0 &&  derivative[a+1] > 0)
      extremum[a] <- -1
    else if (derivative[a] > 0 &&  derivative[a+1] < 0)
      extremum[a] <- 1
  }
  for (a in 1:(length(derivative)*smoothSizeHelper))
  {
    extremum[a] <- 0
  }
  
  for (a in (length(derivative)*(1-smoothSizeHelper)):length(derivative))
  {
    extremum[a] <- 0
  }
  
  if (FALSE)
  {
    results <- rep(FALSE, length(extremum))
    for (a in 1:length(extremum))
    {
      if (extremum[a] == 1)
      {
        for (b in a:1)
        {
          if (ddsmooth[b] > ddsmooth[1] && derivative[b] > 0)
            results[b] <- TRUE
          else
            b <- -1
        }
        for (b in a:length(extremum))
        {
          if (ddsmooth[b] > ddsmooth[1] && derivative[b] < 0)
            results[b] <- TRUE
          else
            b <- length(extremum) + 1
        }
      }
    }
  } else
  {
    results <- rep(FALSE, length(extremum))
    for (a in 1:length(extremum))
    {
      if (extremum[a] == 1)
      {
        results[a] = TRUE
      }
    }
  }
  if (extremumtreshold < 0)
  {
    extremumtreshold <- mean(ddsmooth) + sd(ddsmooth)
    for (a in 1:length(extremum))
    {
      if (extremum[a] == 1 &&  ddsmooth[a] >= extremumtreshold) #&& (ddsmooth[a] - minvalue) >= (diffvalue * extremumtreshold))
      {
        results[a] = TRUE
      }
      else
      {
        results[a] = FALSE
        extremum[a] = 0
      }
    }
  }
  else {
    for (a in 1:length(extremum))
    {
      if (extremum[a] == 1 &&  (ddsmooth[a] - minvalue) >= (diffvalue * extremumtreshold))
      {
        results[a] = TRUE
      }
      else
      {
        results[a] = FALSE
        extremum[a] = 0
      }
    }
  }
  
  vec <- -1
  resultsList <- list()
  countList <- 1
  for (a in 1:length(results))
  {
    if (results[a] == TRUE)
    {
      if (vec[1] == -1)
        vec <- a
      else
        vec <- c(vec,a)
    }
    else
    {
      if (vec[1] != -1)
      {
        resultsList[[countList]] <- vec
        countList <- countList + 1
        vec <- -1
      }
    }
  }
  
  extremumid <- which(extremum %in% 1)
  resultdata<-list(data = dd, smoothdata = ddsmooth, derivative = derivative, extremum = extremum, extremumbool = results, resultsList = resultsList, extremumid = extremumid, extremumtreshold = extremumtreshold)
  return(resultdata)
}




smothdata <- function(dd, smoothSize = 0.2, smoothSizeHelper = 0.1)
{
  require(smoother)
  #smoothSize <- 0.1
  #smoothSize <- 0.2
  #smoothSizeHelper <- 0.1
  
  options('smoother.window' = smoothSize)
  ddsmooth <- smth.gaussian(dd, tails = TRUE)
  minvalue = min(ddsmooth)
  maxvalue = max(ddsmooth)
  diffvalue <- maxvalue - minvalue
  
  #findAllMinimums(x)
  
  derivative <- rep(0, length(ddsmooth))
  for (a in 2:(length(ddsmooth)-1))
  {
    derivative[a] <- (ddsmooth[a+1] - ddsmooth[a-1]) / 2
  }
  
  extremum <- rep(0, length(ddsmooth))
  
  
  for (a in 2:(length(derivative)-1))
  {
    if (derivative[a] < 0 &&  derivative[a+1] > 0)
      extremum[a] <- -1
    else if (derivative[a] > 0 &&  derivative[a+1] < 0)
      extremum[a] <- 1
  }
  for (a in 1:(length(derivative)*smoothSizeHelper))
  {
    extremum[a] <- 0
  }
  
  for (a in (length(derivative)*(1-smoothSizeHelper)):length(derivative))
  {
    extremum[a] <- 0
  }
  
  if (FALSE)
  {
    results <- rep(FALSE, length(extremum))
    for (a in 1:length(extremum))
    {
      if (extremum[a] == 1)
      {
        for (b in a:1)
        {
          if (ddsmooth[b] > ddsmooth[1] && derivative[b] > 0)
            results[b] <- TRUE
          else
            b <- -1
        }
        for (b in a:length(extremum))
        {
          if (ddsmooth[b] > ddsmooth[1] && derivative[b] < 0)
            results[b] <- TRUE
          else
            b <- length(extremum) + 1
        }
      }
    }
  } else
  {
    results <- rep(FALSE, length(extremum))
    for (a in 1:length(extremum))
    {
      if (extremum[a] == 1)
      {
        results[a] = TRUE
      }
    }
  }
  
  vec <- -1
  resultsList <- list()
  countList <- 1
  for (a in 1:length(results))
  {
    if (results[a] == TRUE)
    {
      if (vec[1] == -1)
        vec <- a
      else
        vec <- c(vec,a)
    }
    else
    {
      if (vec[1] != -1)
      {
        resultsList[[countList]] <- vec
        countList <- countList + 1
        vec <- -1
      }
    }
  }
  
  extremumid <- which(extremum %in% 1)
  resultdata<-list(data = dd, smoothdata = ddsmooth, derivative = derivative, extremum = extremum, extremumbool = results, resultsList = resultsList, extremumid = extremumid, 
                   minvalue = minvalue, maxvalue = maxvalue, diffvalue = diffvalue)
  return(resultdata)
}



analyzedta <- function(refdatakinematic, inputdataalignmentkinematic, xx1, xx2, FUN=euc.dist, smoothSize = 0.2, pointscolors=c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                       plottitle = "", path1 = NA, path2 = NA, extremumtreshold = -1, smoothSizeHelper = 0.1, whattodraw = "LeftFoot", plotrgl = FALSE)
{
  #xx1 = refdatakinematicf$RightFoot
  #xx2 = inputdataalignmentkinematicf$RightFoot
  #FUN=euc.dist
  
  #xx1 = refdatakinematicf$dataListRightThigh
  #xx2 = inputdataalignmentkinematicf$dataListRightThigh
  #FUN=angle.betweenCmp
  
  #xx1 = refdatakinematicf$dataListRightLegX
  #xx2 = inputdataalignmentkinematicf$dataListRightLegX
  #FUN=euc.dist
  
  #pointscolors = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  dd <- list()
  if (is.na(path1) || is.na(path2))
  {
    ff <- myDTWCmp(FUN,xx1, xx2)
    path1 <- ff$path1
    path2 <- ff$path2
  }
  #plot(path1, path2)
  
  x1 <- list()
  y1 <- list()
  z1 <- list()
  for (a in 1:length(xx1))
  {
    x1[[a]] <- xx1[[a]][1]
    y1[[a]] <- xx1[[a]][2]
    z1[[a]] <- xx1[[a]][3]
  }
  x1 <- unlist(x1)
  y1 <- unlist(y1)
  z1 <- unlist(z1)
  
  x2 <- list()
  y2 <- list()
  z2 <- list()
  for (a in 1:length(xx2))
  {
    x2[[a]] <- xx2[[a]][1]
    y2[[a]] <- xx2[[a]][2]
    z2[[a]] <- xx2[[a]][3]
  }
  x2 <- unlist(x2)
  y2 <- unlist(y2)
  z2 <- unlist(z2)
  
  #scatterplotanalyzedta(x1, y1, z1, x2, y2, z2, path1, path2)
  
  for (a in 1:length(path1))
  {
    xx <- c(x1[path1[a]], y1[path1[a]], z1[path1[a]])
    yy <- c(x2[path2[a]], y2[path2[a]], z2[path2[a]])
    dd[[a]] <- FUN(xx, yy)
  }
  
  #resultdata <- smothdata(dd, smoothSize = smoothSize)
  #plotsmoothingresults(resultdata)
  
  dd2 <- generatedtwalignment(FUN, x1, y1, z1, x2, y2, z2, path1, path2)
  
  resultdata <- smothdata(dd2, smoothSize = smoothSize, smoothSizeHelper = smoothSizeHelper)
  #resultdata <- smothdata(dd2, smoothSize = smoothSize, extremumtreshold = extremumtreshold, smoothSizeHelper = smoothSizeHelper)
  
  #if (plotrgl)
  #  rglplotanalyzedata(refdatakinematic, inputdataalignmentkinematic, xx1, xx2, path1, path2, resultdata, whattodraw = whattodraw)
  
  
  
  #plot(dd2, col = "blue")
  #lines(dd2, col = "blue")
  
  
  #resultdata <- smothdata(dd2, smoothSize = smoothSize)
  #plotsmoothingresults(resultdata, plottitle)
  if (exists("ff"))
  {
    ff <- c(ff, resultdata)
    return(ff)
  }
  return (resultdata)
}

generatedtwalignment <- function(FUN, x1, y1, z1, x2, y2, z2, path1, path2)
{
  dd2 <- list()
  idd2 <- 1
  idd2prev <- -1
  for (a in 1:length(path1))
  {
    idd2 <- path2[a]
    xx <- c(x1[path1[a]], y1[path1[a]], z1[path1[a]])
    yy <- c(x2[path2[a]], y2[path2[a]], z2[path2[a]])
    
    if (idd2 != idd2prev)
    {
      dd2[[idd2]] <- FUN(xx, yy)
    }
    else
    {
      dd2[[idd2]] <- max(dd2[[idd2prev]], FUN(xx, yy))
    }
    idd2prev <- idd2
  }
  
  dd2 <- unlist(dd2)
  return(dd2)
}

aligninputandrefdata <- function(inputdata, refdata, limbname)
{
  dx3 <- grep("[[:alnum:]]+\\.Dx", colnames(inputdata), ignore.case = TRUE)
  dy3 <- grep("[[:alnum:]]+\\.Dy", colnames(inputdata), ignore.case = TRUE)
  dz3 <- grep("[[:alnum:]]+\\.Dz", colnames(inputdata), ignore.case = TRUE)
  
  
  
  dx <- inputdata[1,paste(limbname, ".Dx", sep = "")] - refdata[1,paste(limbname, ".Dx", sep = "")]
  dy <- inputdata[1,paste(limbname, ".Dy", sep = "")] - refdata[1,paste(limbname, ".Dy", sep = "")]
  dz <- inputdata[1,paste(limbname, ".Dz", sep = "")] - refdata[1,paste(limbname, ".Dz", sep = "")]
  
  for (a in 1:length(dx3))
  {
    inputdata[dx3[a]] <- inputdata[dx3[a]] - dx
    inputdata[dy3[a]] <- inputdata[dy3[a]] - dy
    inputdata[dz3[a]] <- inputdata[dz3[a]] - dz
  }
  return (inputdata)
}

gradient <-function(x)
{
  resultlist <- list()
  for (a in 1:length(x))
  {
    if (a == 1)
      resultlist[[a]] <-x[[a + 1]] - x[[a]]
    else if (a == length(x))
      resultlist[[a]] <- x[[a]] - x[[a - 1]]
    else 
      resultlist[[a]] <- (x[[a + 1]] - x[[a - 1]]) / 2
  }
  return (resultlist)
}

compute.trajectory <- function(x)
{
  tlength <- 0
  for (a in 2:length(x))
  {
    tlength <- tlength + vec.length(x[[a]] - x[[a - 1]])
  }
  return (tlength)
}

##
calculateparameters <- function(dftoanalyze)
{
  #dftoanalyze <- inputdataalignmentkinematicf
  #in seconds
  MovementTime <- length(dftoanalyze$Hips) / 100
  
  #acceleration and velocity m/s
  LeftFootV <- gradient(dftoanalyze$LeftFoot)
  LeftFootA <- gradient(LeftFootV)
  
  #trajectory - kicking and returning
  LeftFootTrajectory <- compute.trajectory(dftoanalyze$LeftFoot) / 100
  
  LeftFootVMax <- max(unlist(lapply(X = LeftFootV, FUN = vec.length)))
  LeftFootAMax <- max(unlist(lapply(X = LeftFootA, FUN = vec.length)))
  
  #Radians
  LeftKneeMax <- max(unlist(dftoanalyze$ListLeftKnee))
  LeftKneeMin <- min(unlist(dftoanalyze$ListLeftKnee))
  
  LeftKneeROM <- LeftKneeMax - LeftKneeMin
  
  LeftKneeV <- as.list(unlist(gradient(dftoanalyze$ListLeftKnee)) * 100)
  LeftKneeVMax <- max(unlist(LeftKneeV))
  
  #acceleration and velocity
  RightFootV <- gradient(dftoanalyze$RightFoot)
  RightFootA <- gradient(RightFootV)
  
  #trajectory - kicking and returning
  RightFootTrajectory <- compute.trajectory(dftoanalyze$RightFoot) / 100
  
  RightFootVMax <- max(unlist(lapply(X = RightFootV, FUN = vec.length)))
  RightFootAMax <- max(unlist(lapply(X = RightFootA, FUN = vec.length)))
  
  RightKneeMax <- max(unlist(dftoanalyze$ListRightKnee))
  RightKneeMin <- min(unlist(dftoanalyze$ListRightKnee))
  
  RightKneeV <- as.list(unlist(gradient(dftoanalyze$ListRightKnee)) * 100)
  RightKneeVMax <- max(unlist(RightKneeV))
  
  
  
  #Radians
  RightKneeROM <- RightKneeMax - RightKneeMin
  resultlist <- list(MovementTime = MovementTime,
                     LeftFootV = LeftFootV,
                     LeftFootA = LeftFootA,
                     LeftFootTrajectory = LeftFootTrajectory,
                     LeftFootVMax = LeftFootVMax,
                     LeftFootAMax = LeftFootAMax,
                     LeftKneeMax = LeftKneeMax,
                     LeftKneeMin = LeftKneeMin,
                     LeftKneeROM = LeftKneeROM,
                     LeftKneeV = LeftKneeV,
                     LeftKneeVMax = LeftKneeVMax,
                     RightFootV = RightFootV,
                     RightFootA = RightFootA,
                     RightFootTrajectory = RightFootTrajectory,
                     RightFootVMax = RightFootVMax,
                     RightFootAMax = RightFootAMax,
                     RightKneeMax = RightKneeMax,
                     RightKneeMin = RightKneeMin,
                     RightKneeROM = RightKneeROM,
                     RightKneeV = RightKneeV,
                     RightKneeVMax = RightKneeVMax)
  return(resultlist)
}  



alignextremum <- function(footddf, kneeadf, analyzerange = 10)
{
  if (length(kneeadf$extremumid) == 0)
    return (kneeadf)
  if (length(footddf$extremumid) == 0)
    return (footddf)
  #analyzerange <- 10
  for (a in 1:length(kneeadf$extremumid))
  {
    found <- FALSE
    for (b in 1:length(footddf$extremumid))
    {
      if (footddf$extremumid[b] - analyzerange <= kneeadf$extremumid[a] &&  kneeadf$extremumid[a] <= footddf$extremumid[b] + analyzerange)
        found <- TRUE
    }
    if (!found)
    {
      kneeadf$resultsList[[a]] <- NA
      kneeadf$extremumbool[kneeadf$extremumid[a]] <- FALSE
      kneeadf$extremumid[a] <- NA
      kneeadf$extremum[a] <- 0
    }
    #kneeadf$extremumbool[a] <- TRUE
  }
  kneeadf$extremumid <- kneeadf$extremumid[!is.na(kneeadf$extremumid)]
  a <- 1
  while (a <= length(kneeadf$resultsList))
  {
    if (is.na((kneeadf$resultsList[[a]])))
    {
      kneeadf$resultsList <- kneeadf$resultsList[-a]
    }
    else {
      a <- a + 1
    }
  }
  return (kneeadf)
}

tresholdresults <- function(extreamdf, extremumtreshold = -1)
{
  if (extremumtreshold < 0)
  {
    extremumtreshold <- mean(extreamdf$smoothdata) + sd(extreamdf$smoothdata)
    for (a in 1:length(extreamdf$extremum))
    {
      if (extreamdf$extremum[a] == 1 &&  extreamdf$smoothdata[a] >= extremumtreshold) #&& (ddsmooth[a] - minvalue) >= (diffvalue * extremumtreshold))
      {
        extreamdf$extremumbool[a] = TRUE
      }
      else
      {
        extreamdf$extremumbool[a] = FALSE
        extreamdf$extremum[a] = 0
      }
    }
  }
  else {
    for (a in 1:length(extreamdf$extremum))
    {
      if (extreamdf$extremum[a] == 1 &&  (extreamdf$smoothdata[a] - extreamdf$minvalue) >= (extreamdf$diffvalue * extremumtreshold))
      {
        extreamdf$extremumbool[a] = TRUE
      }
      else
      {
        extreamdf$extremumbool[a] = FALSE
        extreamdf$extremum[a] = 0
      }
    }
  }
  
  vec <- -1
  resultsList <- list()
  countList <- 1
  for (a in 1:length(extreamdf$extremumbool))
  {
    if (extreamdf$extremumbool[a] == TRUE)
    {
      if (vec[1] == -1)
        vec <- a
      else
        vec <- c(vec,a)
    }
    else
    {
      if (vec[1] != -1)
      {
        resultsList[[countList]] <- vec
        countList <- countList + 1
        vec <- -1
      }
    }
  }
  
  extreamdf$extremumid <- unlist(resultsList)
  
  extreamdf$resultsList <- resultsList
  return (extreamdf)
  
}
