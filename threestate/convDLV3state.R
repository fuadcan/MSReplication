save_file <- "~/Documents/threestate/results/"
data_file <- "~/Documents/threestate/ReplciationFiles/data/"
source("~/Documents/threestate/lnviDM3.R")
source("~/Documents/threestate/return_path.R")
# yearOrRegion <- "Europe+S&P"
convDLV3state <- function(yearOrRegion){

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
  
  
  lowerV <- c( -2,-2,  -2, .6 , .6 , .6 , 0, 0, 0,0.001,-5,-5,-5)
  inits  <- c(1.2, 1, 0.7, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1) 
  upperV <- c(  2, 2,   2,.999,.999,.999,.3,.3,.3,    5, 5, 5, 5)
  
  const_mat <- matrix(0,13,13)
  diag(const_mat) <- 1
  const_mat <- rbind(const_mat,-const_mat)
  dmmvec <- -c(0,0,0,1,0,0,1,0,0,0,0,0,0)
  dmmvec <- c(dmmvec,0,dmmvec,0,dmmvec)[1:39]
  const_mat <- rbind(const_mat,matrix(dmmvec,3,,T))
  const_mat <- cbind(const_mat,c(lowerV,-upperV,rep(-1,3)))
  
  # ress    <- lapply(1:ncol(pdat), function(i) constrOptim(inits, function(p) -lnviD3(p,pdat[,i]), NULL, ui = const_mat[,-14], const_mat[,14]))
  ress    <- lapply(1:ncol(pdat), function(i) {cat(paste0(i,"\n")); constrOptim(inits, function(p) -lnviD3(p,pdat[,i]), NULL, ui = const_mat[,-14], const_mat[,14])})

  save(ress, file = paste0(save_file,yearOrRegion,"_ress.rda"))
  
  return(ress) 
}

