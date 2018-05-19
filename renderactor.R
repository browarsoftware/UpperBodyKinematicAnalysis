library("rgl")

sphere.f <- function(x0 = 0, y0 = 0, z0 = 0, r = 1, n = 101, alpha, ...){
  f <- function(s, t) cbind(r * cos(s) * cos(t) + x0,
                            r * sin(s) * cos(t) + y0, 
                            r * sin(t) + z0)
  persp3d(f, slim = c(0, pi), tlim = c(0, 2*pi), n = n, add = T, alpha = alpha, ...)
}

renderactor <- function(mocapdata, mocapindex, pointscolors, primarycolor, secondarycolor, alpha = 1.0, showspheres = TRUE, linewidth = 50)
{
  
  xdata <- c(mocapdata$Hips.Dx[mocapindex], mocapdata$LeftThigh.Dx[mocapindex], mocapdata$LeftLeg.Dx[mocapindex], mocapdata$LeftFoot.Dx[mocapindex],
             mocapdata$RightThigh.Dx[mocapindex], mocapdata$RightLeg.Dx[mocapindex], mocapdata$RightFoot.Dx[mocapindex],
             mocapdata$SpineLow.Dx[mocapindex], mocapdata$SpineMid.Dx[mocapindex], mocapdata$Chest.Dx[mocapindex],
             mocapdata$LeftShoulder.Dx[mocapindex], mocapdata$LeftArm.Dx[mocapindex], mocapdata$LeftForearm.Dx[mocapindex], mocapdata$LeftHand.Dx[mocapindex],
             mocapdata$RightShoulder.Dx[mocapindex], mocapdata$RightArm.Dx[mocapindex], mocapdata$RightForearm.Dx[mocapindex], mocapdata$RightHand.Dx[mocapindex],
             mocapdata$Neck.Dx[mocapindex], mocapdata$Head.Dx[mocapindex]
  )
  
  ydata <- c(mocapdata$Hips.Dy[mocapindex], mocapdata$LeftThigh.Dy[mocapindex], mocapdata$LeftLeg.Dy[mocapindex], mocapdata$LeftFoot.Dy[mocapindex],
             mocapdata$RightThigh.Dy[mocapindex], mocapdata$RightLeg.Dy[mocapindex], mocapdata$RightFoot.Dy[mocapindex],
             mocapdata$SpineLow.Dy[mocapindex], mocapdata$SpineMid.Dy[mocapindex], mocapdata$Chest.Dy[mocapindex],
             mocapdata$LeftShoulder.Dy[mocapindex], mocapdata$LeftArm.Dy[mocapindex], mocapdata$LeftForearm.Dy[mocapindex], mocapdata$LeftHand.Dy[mocapindex],
             mocapdata$RightShoulder.Dy[mocapindex], mocapdata$RightArm.Dy[mocapindex], mocapdata$RightForearm.Dy[mocapindex], mocapdata$RightHand.Dy[mocapindex],
             mocapdata$Neck.Dy[mocapindex], mocapdata$Head.Dy[mocapindex]
  )
  
  zdata <- c(mocapdata$Hips.Dz[mocapindex], mocapdata$LeftThigh.Dz[mocapindex], mocapdata$LeftLeg.Dz[mocapindex], mocapdata$LeftFoot.Dz[mocapindex],
             mocapdata$RightThigh.Dz[mocapindex], mocapdata$RightLeg.Dz[mocapindex], mocapdata$RightFoot.Dz[mocapindex],
             mocapdata$SpineLow.Dz[mocapindex], mocapdata$SpineMid.Dz[mocapindex], mocapdata$Chest.Dz[mocapindex],
             mocapdata$LeftShoulder.Dz[mocapindex], mocapdata$LeftArm.Dz[mocapindex], mocapdata$LeftForearm.Dz[mocapindex], mocapdata$LeftHand.Dz[mocapindex],
             mocapdata$RightShoulder.Dz[mocapindex], mocapdata$RightArm.Dz[mocapindex], mocapdata$RightForearm.Dz[mocapindex], mocapdata$RightHand.Dz[mocapindex],
             mocapdata$Neck.Dz[mocapindex], mocapdata$Head.Dz[mocapindex]
  )
  pointscolorsHelper <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  segmentList <- c(1,2,2,3,3,4,
                   1,5,5,6,6,7,
                   1,8,8,9,9,10,10,11,11,12,12,13,13,14,
                   10,15,15,16,16,17,17,18,
                   10,19,19,20)
  for (a in seq(1,length(segmentList),2))
  {
    if (pointscolors[segmentList[a+1]] == 0)
      colorHelp = primarycolor
    else
    {
      colorHelp = secondarycolor
      pointscolorsHelper[segmentList[a+1]] = 1
      pointscolorsHelper[segmentList[a]] = 1
    }
    segments3d(c(xdata[segmentList[a]], xdata[segmentList[a+1]]),
               c(ydata[segmentList[a]], ydata[segmentList[a+1]]),
               c(zdata[segmentList[a]], zdata[segmentList[a+1]]), color = colorHelp, lwd=linewidth, alpha = alpha)
  }
  if (showspheres)
  {
    for (a in 1:length(xdata))
    {
      if (pointscolorsHelper[a] == 0)
        colorHelp = primarycolor
      else
        colorHelp = secondarycolor
      sphere.f(xdata[a], ydata[a], zdata[a], color = colorHelp, r=5, n=31, alpha = alpha)
      #rgl::spheres3d(xdata[a], ydata[a], zdata[a], color = colorHelp, r=5)
    }
  }
  returndata <- list()
  returndata[[1]] <- xdata
  returndata[[2]] <- ydata
  returndata[[3]] <- zdata
  returndata[[4]] <- pointscolorsHelper
  return (returndata)
}


drawscene <- function(inputdata, refdata, inputdataid, refdataid, featurename, pointscolors = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), alpha = 1)
{
  #dx3 <- grep("[[:alnum:]]+\\.Dx", colnames(inputdata), ignore.case = TRUE)
  #dy3 <- grep("[[:alnum:]]+\\.Dy", colnames(inputdata), ignore.case = TRUE)
  #dz3 <- grep("[[:alnum:]]+\\.Dz", colnames(inputdata), ignore.case = TRUE)
  
  
  
  #dx <- inputdata[1,"LeftFoot.Dx"] - refdata[1,"LeftFoot.Dx"]
  #dy <- inputdata[1,"LeftFoot.Dy"] - refdata[1,"LeftFoot.Dy"]
  #dz <- inputdata[1,"LeftFoot.Dz"] - refdata[1,"LeftFoot.Dz"]
  
  #for (a in 1:length(dx3))
  #{
  #  inputdata[dx3[a]] <- inputdata[dx3[a]] - dx
  #  inputdata[dy3[a]] <- inputdata[dy3[a]] - dy
  #  inputdata[dz3[a]] <- inputdata[dz3[a]] - dz
  #}
  
  
  library("rgl")
  rgl.open() # Open a new RGL device
  #pointscolors <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  returndata1 <- renderactor(inputdata, inputdataid, pointscolors, "lightgray", "green", alpha)
  returndata2 <- renderactor(refdata, refdataid, pointscolors, "darkgray", "red", alpha)
  
  borderx <- c(-100, -100, 100, 100, -100, -100, 100, 100) * 1
  bordery <- c(200, -20, 200, -20, 200, -20, 200, -20) * 1
  borderz <- c(-100, -100, -100, -100, 100, 100, 100, 100) * 1
  
  #rgl.material(alpha = 0.5)
  spheres3d(borderx, bordery, borderz, col = "red", r=0.01)
  
  for (a in 1:length(pointscolors))
  {
    if (pointscolors[a] != 0)
      lines3d(c(returndata1[[1]][a],returndata2[[1]][a]),
              c(returndata1[[2]][a],returndata2[[2]][a]),
              c(returndata1[[3]][a],returndata2[[3]][a]), color = "yellow", lwd = 20)
  }
  
  y <- min(refdata[1,featurename])
  planes3d(0,-1,0,y, col = 'orange', alpha = 0.3)
  axes3d(col="white")
  
}