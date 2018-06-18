save_file <- "~/Documents/threestate/results/"
data_file <- "~/Documents/threestate/ReplciationFiles/data/"

modifyResults <- function(yearOrRegion){
  
  # yearOrRegion <- "Europe+S&P"
  if(is.numeric(yearOrRegion)){year      <- yearOrRegion
  filename  <- paste0(data_file, "madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else if(yearOrRegion=="Maddison"){
    filename  <- paste0(data_file,"madisonFromMaxNA2-1950.csv")
    z         <- read.table(filename,header = T,sep = " ")
    z         <- data.matrix(z)
  } else {fname1 <- paste0(data_file,"mds_G7-1950.csv")
  fname2 <- paste0(data_file,"mds_Europe-1950.csv")
  fname3 <- paste0(data_file,"mds_S&P-1950.csv")
  z_g7         <- read.table(fname1,header = T,sep = ";"); z_g7 <- data.matrix(z_g7)
  z_eu         <- read.table(fname2,header = T,sep = ";"); z_eu <- data.matrix(z_eu)
  z_sp         <- read.table(fname3,header = T,sep = ";"); z_sp <- data.matrix(z_sp)
  
  
  if(yearOrRegion=="Europe+G7"|yearOrRegion=="G7+Europe")
  {z<- matrix(c(z_g7,z_eu),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_eu))} else
    if(yearOrRegion=="Europe+S&P"|yearOrRegion=="S&P+Europe")
    {z<- matrix(c(z_eu,z_sp),dim(z_eu)[1],); colnames(z)<- c(colnames(z_eu),colnames(z_sp))} else
      if(yearOrRegion=="G7+S&P"|yearOrRegion=="S&P+G7")
      {z<- matrix(c(z_g7,z_sp),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_sp))} else {stop("Unknown Country List")}
  if(sum(duplicated(t(z)))!=0){z<- z[,-which(duplicated(t(z)))]}
  }
  
  
  pdat   <- apply(combn(ncol(z),2),2, function(x) log(z[,x[1]]) - log(z[,x[2]]))
  pnames <- apply(combn(colnames(z),2), 2, function(n) paste0(n,collapse = " - ")) 
  ress   <- get(load(paste0(save_file,yearOrRegion,"_ress.rda")))

  inds   <- which(sapply(ress, function(r) r$convergence) == 0)
  
  lowerV <- c( -2,-2,  -2, .6 , .6 , .6 , 0, 0, 0,0.001,-5,-5,-5)
  inits2 <- c(1.7, 1.5, 1, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1) 
  upperV <- c(  2, 2,   2,.999,.999,.999,.3,.3,.3,    5, 5, 5, 5)
  
  const_mat <- matrix(0,13,13)
  diag(const_mat) <- 1
  const_mat <- rbind(const_mat,-const_mat)
  dmmvec <- -c(0,0,0,1,0,0,1,0,0,0,0,0,0)
  dmmvec <- c(dmmvec,0,dmmvec,0,dmmvec)[1:39]
  const_mat <- rbind(const_mat,matrix(dmmvec,3,,T))
  const_mat <- cbind(const_mat,c(lowerV,-upperV,rep(-1,3)))
  
  # ress     <- lapply(1:ncol(pdat), function(i) constrOptim(inits, function(p) -lnviD3(p,pdat[,i]), NULL, ui = const_mat[,-14], const_mat[,14]))
  newress <- lapply(inds, function(i) {cat(paste0(i,"\n")); constrOptim(inits2, function(p) -lnviD3(p,pdat[,i]), NULL, ui = const_mat[,-14], const_mat[,14])})
  
  modress       <- ress
  modress[inds] <- newress 
  save(newress, file = paste0(save_file,"../modified_results/",yearOrRegion,"_modifyingress.rda"))
  save(modress, file = paste0(save_file,"../modified_results/",yearOrRegion,"_modress.rda"))
  
  return(ress) 
}
# 
# ress <- get(load("Documents/threestate/results/1930ress.rda"))
# inds <- sapply(ress, function(r) r$convergence) == 0
# 
# dat    <- read.csv("~/Documents/threestate/ReplciationFiles/data/madisonFrom-1930.csv",sep=";")
# pdat   <- apply(combn(ncol(dat),2),2, function(x) log(dat[,x[1]]) - log(dat[,x[2]]))
# pnames <- apply(combn(names(dat),2), 2, function(n) paste0(n,collapse = " - "))
# 
# tempres <- ress[inds][[1]]
# tempdat <- pdat[,inds][,1]
# 
# lowerV <- c( -2,-2,  -2, .6 , .6 , .6 , 0, 0, 0,0.001,-5,-5,-5)
# inits  <- c(1.2, 1, 0.7, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1)
# inits2 <- c(1.7, 1.5, 1, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1)
# upperV <- c(  2, 2,   2,.999,.999,.999,.3,.3,.3,    5, 5, 5, 5)
# 
# const_mat <- matrix(0,13,13)
# diag(const_mat) <- 1
# const_mat <- rbind(const_mat,-const_mat)
# dmmvec <- -c(0,0,0,1,0,0,1,0,0,0,0,0,0)
# dmmvec <- c(dmmvec,0,dmmvec,0,dmmvec)[1:39]
# const_mat <- rbind(const_mat,matrix(dmmvec,3,,T))
# const_mat <- cbind(const_mat,c(lowerV,-upperV,rep(-1,3)))
# 
# tempres2 <- constrOptim(inits2, function(p) -lnviD3(p,tempdat), NULL, ui = const_mat[,-14], const_mat[,14])
# ds <- sapply(1:length(ress), function(x) matrix(ress[[x]]$par[1:3],1,) %*% return_path(ress[[x]]$par,pdat[,x]))

# 
# selected <- c("Finland - Greece", "Germany - Italy", "Canada - USA", "Ecuador - Guatemala")
# 
# sdat <- pdat[,pnames %in% selected]
# 
# 
# lowerV <- c( -2,-2,  -2, .6 , .6 , .6 , 0, 0, 0,0.001,-5,-5,-5)
# inits  <- c(1.2, 1, 0.7, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1) 
# upperV <- c(  2, 2,   2,.999,.999,.999,.3,.3,.3,    5, 5, 5, 5)
# 
# 
# const_mat <- matrix(0,13,13)
# diag(const_mat) <- 1
# const_mat <- rbind(const_mat,-const_mat)
# dmmvec <- -c(0,0,0,1,0,0,1,0,0,0,0,0,0)
# dmmvec <- c(dmmvec,0,dmmvec,0,dmmvec)[1:39]
# const_mat <- rbind(const_mat,matrix(dmmvec,3,,T))
# const_mat <- cbind(const_mat,c(lowerV,-upperV,rep(-1,3)))
# 
# ress    <- lapply(1:ncol(pdat), function(i) constrOptim(inits, function(p) -lnviD3(p,pdat[,i]), NULL, ui = const_mat[,-14], const_mat[,14]))
# save(ress, file = "~/threestate/g7ress.rda")
# thepath <- return_path(res$par,ser)
# plot(apply(res$par[1:3]*thepath,2,sum),type="l")
# 
# res4 <- constrOptim(inits, function(p) -lnviD3(p,sdat[,4]), NULL, ui = const_mat[,-14], const_mat[,14])
# 
# res2 <- res
# plot(1930:2010,matrix(res2$par[1:3],1,) %*% return_path(res2$par,sdat[,2]), type="l")
# 
# 
# ds <- sapply(1:length(ress), function(x) matrix(ress[[x]]$par[1:3],1,) %*% return_path(ress[[x]]$par,pdat[,x]))
# 
# x <- 17
# pnames[x]
# plot(1950:2010,matrix(ress[[x]]$par[1:3],1,) %*% return_path(ress[[x]]$par,pdat[,x]), type="l")
# 
# plot(1930:2010,matrix(res4$par[1:3],1,) %*% return_path(res4$par,sdat[,4]), type="l")
