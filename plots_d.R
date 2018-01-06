# Reading outputs
ress_d   <- lapply(dir("output","_D_.*resALL.rda"), function(d) get(load(paste0("output/",d))))
# Reformating results
ress_d <- lapply(ress_d, function(r) correctRes(r[,1:7]))
# data names and pair panels
dnames     <- gsub("d_|_D_resALL.rda","",dir("output","_D_.*resALL.rda"))
pdats      <- lapply(dnames, gen_pdat)

# d estimations vs year for each pair 
dss     <- lapply(1:length(ress_d), function(i) sapply(1:nrow(ress_d[[i]]), function(x) {matrix(ress_d[[i]][x,1:2],1,) %*% dlvPath_d(ress_d[[i]][x,1:6],pdats[[i]][,x])}))

# proportions of d < 1 at each time point
d_plots     <- lapply(dss, function(ds) apply(ds<1,1,mean)) 

# Plots of prop(d<1) for each dataset
for(i in 1:length(dnames)){
  jpeg(paste0("rejplots/d_lt_1_",dnames[i],"_D.jpg"))
  plot((2010 - length(d_plots[[i]]) + 1):2010, d_plots[[i]],ylab="proportion (d < 1)",xlab="year",main=paste0("Proportions of d < 1, ",dnames[i]))
  dev.off()
}

# Code for plotting d estimations vs year of given pair for given dataset
plot_specific <- function(n1,n2,datname){
  pairname <- paste0(n1," - ",n2)
  datind   <- which(datname==dnames)
  pnames   <- colnames(gen_pdat(datname))
  serind   <- which(pairname == pnames)
  ser      <- dss[[datind]][,serind]; maintext <- paste0(pnames[serind], " ", dnames[datind])
  jpeg(paste0("pairplots/",dnames[datind],"_",pnames[serind],"_D.jpg"))
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
  dev.off()
  plot((2010 - length(ser)+ 1):2010, ser,type="l",xlab="year",ylab="d",main=maintext)
}

# Europe Case
plot_specific("Finland","Greece",1930)
plot_specific("France","Sweden",1930)
plot_specific("Germany","Italy",1930)
plot_specific("Germany","Netherlands",1930)

# Developed Case
plot_specific("Belgium","Netherlands",1930)
plot_specific("Sweden","USA",1930)
plot_specific("Netherlands","UK",1930)
plot_specific("Canada","USA",1930)

# Developing Case
plot_specific("Chile","Morocco",1950)
plot_specific("Ecuador","Guatemala",1930)
plot_specific("Peru","Philippines",1950)
plot_specific("Peru","South.Africa",1950)

