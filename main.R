if(Sys.info()["sysname"]=="Linux"){setwd("~/Documents/MSReplication/")} else {setwd("~/MSReplication/")}

source("convDLV_d.R")
source("convDLV_dm!!!!.R") # remove it !!!!
source("utils.R")
source("dlvPath.R")
source("plotAll.R")
source("plotRejs.R")
source("lnviD2.R")
if(Sys.info()["sysname"]=="Linux")
  {library("parallel");mclapply.hack <- mclapply} else {source("mclapplyhack.R")}
dir.create("output",showWarnings = F)
dir.create("results",showWarnings = F)
dir.create("rejplots",showWarnings = F)
dir.create("pairplots",showWarnings = F)

# Estimating d's for only d switching analysis
res1930   <- convDLV_d(1930)
res1940   <- convDLV_d(1940)
res1950   <- convDLV_d(1950)
reseu     <- convDLV_d("Europe")
resg7eu   <- convDLV_d("Europe+G7")
resg7sp   <- convDLV_d("G7+S&P")
resspeu   <- convDLV_d("Europe+S&P")

# Estimating d's for both d and mu switching analysis
res1930   <- convDLV_dm(1930)
res1940   <- convDLV_dm(1940)
res1950   <- convDLV_dm(1950)
reseu     <- convDLV_dm("Europe")
resg7eu   <- convDLV_dm("Europe+G7")
resg7sp   <- convDLV_dm("G7+S&P")
resspeu   <- convDLV_dm("Europe+S&P")

# Reading outputs
ress_d   <- lapply(dir("output","_D_.*resALL.rda"), function(d) get(load(paste0("output/",d))))
ress_dm  <- lapply(dir("output","_DM_.*resALL.rda"), function(d) get(load(paste0("output/",d))))

# Converting to tables
ress_dm <- lapply(ress_dm, totable)
ress_dm <- lapply(ress_dm, function(r) r[,-7])

# Reformating results
ress_d  <- lapply(ress_d, function(r) correctRes(r[,1:7]))
ress_dm <- lapply(ress_dm, correctRes)


# data names and pair panels
dnames     <- gsub("d_|_D_resALL.rda","",dir("output","_D_.*resALL.rda"))
pdats      <- lapply(dnames, gen_pdat)

# State switching series for dm and d
pathss_DM   <- lapply(1:length(ress_dm), function(i) lapply(1:nrow(ress_dm[[i]]), function(x) dlvPath_dm(ress_dm[[i]][x,-8],pdats[[i]][,x])))
pathss_D    <- lapply(1:length(ress_d), function(i) lapply(1:nrow(ress_d[[i]]), function(x) dlvPath_d(ress_d[[i]][x,-7],pdats[[i]][,x])))

# Correcting results
ischange  <- lapply(pathss_DM, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,temp); return(temp)})))
parss_dm  <- lapply(1:length(ress_dm),  function(i)  ress_dm[[i]][,-8] * ischange[[i]])

ischange  <- lapply(pathss_D, function(paths) t(sapply(paths, function(p) {temp <- apply(p,1,sum) > 0; temp <- c(temp,temp,T,T); return(temp)})))
parss_d   <- lapply(1:length(ress_d),  function(i)  ress_d[[i]][,-7] * ischange[[i]])

# Reporting outputs
rep_dm <- do.call(rbind,lapply(parss_dm,report))
rep_dm <- data.frame(c(sapply(dnames, function(n) rep(n,5))),rownames(rep_dm),rep_dm)

rep_d  <- do.call(rbind,lapply(parss_d,report))
rep_d  <- data.frame(c(sapply(dnames, function(n) rep(n,5))),rownames(rep_d),rep_d)

colnames(rep_dm)  <- c("Data","Conv?",colnames(rep_dm)[-(1:2)])
colnames(rep_d)   <- c("Data","Conv?",colnames(rep_d)[-(1:2)])

# Writing reports
write.csv(rep_dm, "results/report_dm.csv")
write.csv(rep_d , "results/report_d.csv")

# Plotting all path graphs
lapply(c(1930, 1940,1950), plotAll)
lapply(c("G7+S&P","Europe+G7","Europe+S&P"),plotAll)
lapply(c(1930, 1940), plotRejs)
lapply(c("G7+S&P","Europe+G7","Europe+S&P"), plotRejs)
