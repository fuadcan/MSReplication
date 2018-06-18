save_file <- "~/Documents/threestate/results/"
data_file <- "~/Documents/threestate/ReplciationFiles/data/"
source("~/Documents/threestate/lnviDM3.R")
source("~/Documents/threestate/return_path.R")
# yearOrRegion <- "G7+S&P"

whatwehave <- dir("Documents/threestate/results/")


fname1 <- paste0(data_file,"mds_G7-1950.csv")
fname2 <- paste0(data_file,"mds_Europe-1950.csv")
fname3 <- paste0(data_file,"mds_S&P-1950.csv")
z_g7   <- read.table(fname1,header = T,sep = ";"); z_g7 <- data.matrix(z_g7)
z_eu   <- read.table(fname2,header = T,sep = ";"); z_eu <- data.matrix(z_eu)
z_sp   <- read.table(fname3,header = T,sep = ";"); z_sp <- data.matrix(z_sp)

pnames_g7 <- colnames(z_g7)
pnames_eu <- colnames(z_eu)
pnames_sp <- colnames(z_sp)

c("Europe+S&P","G7+Europe")
whatwehave  <- list(c(pnames_eu,pnames_sp), c(pnames_g7,pnames_eu))
whatwehave  <- lapply(whatwehave, function(c) c[!duplicated(c)])
pairswehave <- lapply(whatwehave, function(c) apply(combn(c,2), 2, function(n) paste0(sort(n),collapse = " - ")))
pairsneeded <- c(pnames_g7,pnames_sp)
pairsneeded <- pairsneeded[!duplicated(pairsneeded)]
pairsneeded <- apply(combn(pairsneeded,2), 2, function(n) paste0(sort(n),collapse = " - "))

pairsneeded <- setdiff(pairsneeded, unlist(pairswehave))
cneeded     <- unlist(strsplit(pairsneeded, " - "))
cneeded     <- cneeded[!duplicated(cneeded)]

dat    <- cbind(z_sp,z_eu,z_g7)
dat    <- dat[,colnames(dat) %in% cneeded]
pdat   <- apply(combn(ncol(dat),2),2, function(x) log(dat[,x[1]]) - log(dat[,x[2]]))
pnames <- apply(combn(colnames(dat),2), 2, function(n) paste0(sort(n),collapse = " - ")) 
pdat   <- pdat[,pnames %in% pairsneeded]

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
  
  dat <- cbind(z_g7,z_sp)
  dat <- dat[,!duplicated(colnames(dat))]
  pdat   <- apply(combn(ncol(dat),2),2, function(x) log(dat[,x[1]]) - log(dat[,x[2]]))
  pnames <- apply(combn(colnames(dat),2), 2, function(n) paste0(sort(n),collapse = " - ")) 
  pdat   <- pdat[,pnames %in% pairsneeded]
  
  
  return_pairs <- function(yearOrRegion){
    fname1 <- paste0(data_file,"mds_G7-1950.csv")
    fname2 <- paste0(data_file,"mds_Europe-1950.csv")
    fname3 <- paste0(data_file,"mds_S&P-1950.csv")
    z_g7   <- read.table(fname1,header = T,sep = ";"); z_g7 <- data.matrix(z_g7)
    z_eu   <- read.table(fname2,header = T,sep = ";"); z_eu <- data.matrix(z_eu)
    z_sp   <- read.table(fname3,header = T,sep = ";"); z_sp <- data.matrix(z_sp)
    
    pnames_g7 <- colnames(z_g7)
    pnames_eu <- colnames(z_eu)
    pnames_sp <- colnames(z_sp)
    
    if(yearOrRegion=="Europe+G7"|yearOrRegion=="G7+Europe")
    {z<- matrix(c(z_g7,z_eu),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_eu))} else
      if(yearOrRegion=="Europe+S&P"|yearOrRegion=="S&P+Europe")
      {z<- matrix(c(z_eu,z_sp),dim(z_eu)[1],); colnames(z)<- c(colnames(z_eu),colnames(z_sp))} else
        if(yearOrRegion=="G7+S&P"|yearOrRegion=="S&P+G7")
        {z<- matrix(c(z_g7,z_sp),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_sp))} else {stop("Unknown Country List")}
    if(sum(duplicated(t(z)))!=0){z<- z[,-which(duplicated(t(z)))]}
    pnames <- apply(combn(colnames(z),2), 2, function(n) paste0(sort(n),collapse = " - ")) 
    
    return(pnames)
    }
    
  ress_eusp <- get(load("Documents/threestate/results/Europe+S&P_ress.rda"))
  ress_g7eu <- get(load("Documents/threestate/results/G7+Europe_ress.rda"))
  ress_rest <- get(load("Documents/threestate/results/restofG7nSP_ress.rda"))
  
  names(ress_g7eu) <- return_pairs("Europe+G7")
  names(ress_eusp) <- return_pairs("Europe+S&P")
  names(ress_rest) <- pairsneeded
  
  poolofress <- c(ress_eusp,ress_g7eu,ress_rest)
  poolofress <- poolofress[!duplicated(names(poolofress))]
  ress_g7sp  <- poolofress[names(poolofress) %in% return_pairs("G7+S&P")]
  
  save(ress_g7sp, file = "Documents/threestate/results/G7+S&P_ress.rda")
  dat    <- cbind(z_sp,z_eu); dat <- dat[,!duplicated(colnames(dat))]
  pnames <- apply(combn(colnames(dat),2), 2, function(n) paste0(sort(n),collapse = " - ")) 
  names(ress_eusp) <- pnames
  
  
  
  ress <- lapply(1:length(pnames), function(x) NULL)
  names(ress) <- pnames
  ress[match(pairsneeded, pnames)] <- ress_rest
  
  outress <- c(ress, )
  
  save(ress, file = paste0(save_file,"restofG7nSP_ress.rda"))
  
  return(ress) 
}
