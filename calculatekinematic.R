correctMe3D <- function(startC, stopC, signal1, signal2, signal3, all, dx_vector, dy_vector, dz_vector)
{
  diff2 <- all[stopC, signal1] - all[startC, signal1] 
  all[stopC:length(all$Hips.Dy),dz_vector] = all[stopC:length(all$Hips.Dy),dz_vector] - diff2
  
  diff2 <- all[stopC, signal2] - all[startC, signal2] 
  all[stopC:length(all$Hips.Dy),dx_vector] = all[stopC:length(all$Hips.Dy),dx_vector] - diff2
  
  
  diff2 <- all[stopC, signal3] - all[startC, signal3] 
  all[stopC:length(all$Hips.Dy),dy_vector] = all[stopC:length(all$Hips.Dy),dy_vector] - diff2
  
  return (all)
}



calculatekinematic <- function(dd, bodypartname)
{
  #dd <- mocapdata1
  #bodypartname<-"LeftFoot"
  #window_size <- 10 / length(dd$RightFoot.Dx)
  window_size <- 0.05
  dx_vector <- names(dd)[grepl(".Dx",names(dd))]
  dy_vector <- names(dd)[grepl(".Dy",names(dd))]
  dz_vector <- names(dd)[grepl(".Dz",names(dd))]
  

  require(smoother)
  #zzR <- smth(sqrt(dd$RightFoot.Dx ^ 2 + dd$RightFoot.Dy ^ 2 + dd$RightFoot.Dz ^ 2),window = window_size,method = "gaussian") #SMOOTHING
  #zzL <- smth(sqrt(dd$LeftFoot.Dx ^ 2 + dd$LeftFoot.Dy ^ 2 + dd$LeftFoot.Dz ^ 2),window = window_size,method = "gaussian") #SMOOTHING
  #return(dd)
  #plot(zzR, type = 'n', xlab = "Acceleration [g]", ylab="Y [cm]")
  #title("Heian Shodan")
  #lines(zzR, col="red")
  #lines(zzL, col="blue")
  #return(dd)
  #zzR[is.na(zzR)] <- 0
  #zzL[is.na(zzL)] <- 0  
  
  for (a in 1:(length(dd$RightFoot.Dx) - 1))
  {
    
    dd <- correctMe3D(a, a + 1, 
                      paste(bodypartname, ".Dz", sep = ""),
                      paste(bodypartname, ".Dx", sep = ""),
                      paste(bodypartname, ".Dy", sep = ""), dd,
                      dx_vector, dy_vector, dz_vector)
    
    #dd <- correctMe3D(a, a + 1, "LeftFoot.Dz", "LeftFoot.Dx", "LeftFoot.Dy", dd)
    
    #if (zzR[a] > zzL[a] && zzR[a + 1] > zzL[a + 1])
    #{
    #  dd <- correctMe3D(a, a + 1, "LeftFoot.Dz", "LeftFoot.Dx", "LeftFoot.Dy", dd)
    #}
    #else if (zzR[a] < zzL[a] && zzR[a + 1] < zzL[a + 1])
    #{
    #  dd <- correctMe3D(a, a + 1, "RightFoot.Dz", "RightFoot.Dx", "RightFoot.Dy", dd)
    #}
  }
  #dd <- calculatekinematiccorrection(dd)
  return (dd)
}

MoveToTheGround <- function (all, a, groundPosition,  signal1)
{
  diff2 <- all[a, signal1] - groundPosition
  all[a,dy_vector] = all[a,dy_vector] - diff2
  return (all)
}


dyEps <- 5

calculatekinematiccorrection <- function(dd)
{
  window_size <- 100 / length(dd$RightFoot.ax)
  
  require(smoother)
  zzR <- smth(dd$RightFoot.Dy,window = window_size,method = "gaussian") #SMOOTHING
  #plot(zzR, type = 'n', xlab = "Time [ms]", ylab="Y [cm]")
  #title("Heian Shodan")
  #lines(zzR, col="red")
  zzL <- smth(dd$LeftFoot.Dy,window = window_size,method = "gaussian") #SMOOTHING
  #lines(zzL, col="blue")
  
  zzR[is.na(zzR)] <- 0
  zzL[is.na(zzL)] <- 0
  
  for (a in 1:(length(dd$RightFoot.ax) - 1))
  {
    if (abs(zzR[a] - zzL[a]) > dyEps)
    {
      if (zzR[a] > zzL[a] && zzR[a + 1] > zzL[a + 1])
      {
        dd <- correctMe3D(a, a + 1, "LeftFoot.Dz", "LeftFoot.Dx", "LeftFoot.Dy", dd)
      }
      else if (zzR[a] < zzL[a] && zzR[a + 1] < zzL[a + 1])
      {
        dd <- correctMe3D(a, a + 1, "RightFoot.Dz", "RightFoot.Dx", "RightFoot.Dy", dd)
      }
    }
  }
  
  groundPosition <- dd$RightFoot.Dy[1]
  for (a in 1:(length(dd$RightFoot.ax)))
  {
    if (zzR[a] >= zzL[a])
    {
      #przesuñ ca³¹ sylwetkê, aby Y stopy dotyka³o ziemi
      dd <- MoveToTheGround(dd, a, groundPosition, "LeftFoot.Dy")
    }
    else if (zzR[a] < zzL[a])
    {
      dd <- MoveToTheGround(dd, a, groundPosition, "RightFoot.Dy")
    }
  }
  
  
  zzR <- smth(dd$RightFoot.Dy,window = window_size,method = "gaussian") #SMOOTHING
  zzL <- smth(dd$LeftFoot.Dy,window = window_size,method = "gaussian") #SMOOTHING
  #plot(zzL, type = 'n', xlab = "Time [ms]", ylab="Y [cm]")
  #title("Heian Shodan")
  #lines(zzR, col="red")
  zzL <- smth(dd$LeftFoot.Dy,window = window_size,method = "gaussian") #SMOOTHING
  #lines(zzL, col="blue")
  
  return (dd)
}

