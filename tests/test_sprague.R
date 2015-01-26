## test sprague index:
library(simPop)
x <- data.frame(age=as.factor(c(
  "0-4",
  "5-9","10-14","15-19", "20-24",
  "25-29","30-34","35-39","40-44","45-49",
  "50-54","55-59","60-64","65-69","77-74","75-79","80+"
    )),
  pop=c(1971990, 2095820,2157190, 2094110,2116580,   2003840, 1785690,
        1502990, 1214170, 796934,  627551,  530305, 488014,
        364498, 259029,158047,  125941)
)

# debugonce(sprague)
s  <- sprague(x[,2])


if(!all.equal(sum(s), sum(x[,2]))) stop("not identical")

################ debug sa
# 
# rm(list=ls())
# library(Rcpp)
# library(simPop)
# data(eusilcP)
# data(eusilcS)
# 
# ## synth pop:
# pop <- eusilcP
# colnames(pop)[3] <- "hhsize"
# 
# ## donor data:
# donors <- sampHH(pop, strata="region", hsize="hhsize")
# 
# ## microdata (1) and donor (0)
# pop <- rbind(pop, donors)
# pop$weights <- rep(c(1,0), c(nrow(eusilcP), nrow(donors)))
# pop$new.weights <- pop$weights
# 
# index <- pop$weights == 1
# tab <- table(pop$region[index],pop$gender[index], pop$hhsize[index])
# 
# ## create target marginals
# totals <- tableWt(eusilcS[,c("db040", "rb090","hsize")], weights=eusilcS$rb050)
# totals <- nrow(pop[index])/sum(totals) * totals
# totals <- ceiling(totals)
# totals <- as.data.table(totals)
# setnames(totals, c("region","gender","hhsize","N"))
# 
# data <- data.table(pop, keep.rownames=TRUE)
# totals <- totals
# hid <- "hid"
# parameter <- c("gender","hhsize")
# split <- "region"
# temp <- 10
# eps.factor <-  0.02
# maxiter <- 250
# temp.cooldown <- 0.90
# factor.cooldown <- 0.95
# min.temp <- 10^-2
# sample <- TRUE
# parallel <- TRUE
# verbose <- FALSE
# 
# sourceCpp("src/calibPop.cpp") 
# source("R/calibPop.R")
# 
# data <- calibPop_cpp(data, totals, hid, parameter, split, temp = temp, eps.factor = eps.factor, 
#                      maxiter=maxiter, temp.cooldown = temp.cooldown, factor.cooldown = factor.cooldown, min.temp = min.temp, verbose=verbose, parallel=parallel)
# 
# # compare with fixed margins
# index <- data$weights == 1
# tab_old <- as.data.table(table(data$region[index],data$gender[index], data$hhsize[index]))
# setnames(tab_old, c("region",parameter,"before"))
# setkeyv(tab_old,c("region",parameter))
# 
# index <- data$new.weights == 1
# tab <- as.data.table(table(data$region[index],data$gender[index], data$hhsize[index]))
# setnames(tab, c("region",parameter,"after"))
# setkeyv(tab,c("region",parameter))
# 
# setnames(totals, c("region",parameter,"goal"))
# setkeyv(totals,c("region",parameter))
# 
# res <- totals[tab_old]
# res <- res[tab]
# res

############### debug ipu
# 
# ## goal: calculate household weights that should also fulfil person-type constraints
# rm(list=ls())
# library(simPop)
# # load sample and population data
# data(eusilcS)
# data(eusilcP)
# 
# # variable generation and preparation
# eusilcS$hsize <- factor(eusilcS$hsize)
# 
# # make sure, factor levels in sample and population match
# eusilcP$region <- factor(eusilcP$region, levels = levels(eusilcS$db040))
# eusilcP$gender <- factor(eusilcP$gender, levels = levels(eusilcS$rb090))
# eusilcP$hsize  <- factor(eusilcP$hsize , levels = levels(eusilcS$hsize))
# 
# # generate input matrix
# # we want to adjust to variable "db040" (region) as household variables and
# # variable "rb090" (gender) as individual information
# samp <- data.table(eusilcS)
# pop <-  data.table(eusilcP)
# setkeyv(samp, "db030")
# hh <- samp[!duplicated(samp$db030),]
# hhpop <- pop[!duplicated(pop$hid),]
# 
# # reg contains for each region the number of households
# reg <- data.table(model.matrix(~db040 +0, data=hh))
# # hsize contains for each household size the number of households
# hsize <- data.table(model.matrix(~factor(hsize) +0, data=hh))
# 
# # aggregate persons-level characteristics per household
# # gender contains for each household the number of males and females
# gender <- data.table(model.matrix(~db030+rb090 +0, data=samp))
# setkeyv(gender, "db030")
# gender <- gender[, lapply(.SD, sum), by = key(gender)]
# 
# # bind together and use it as input
# inp <- cbind(reg,
#              hsize,
#              gender)
# # the totals we want to calibrate to
# con <- c(
#   as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$region)),
#   as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$hsize)),
#   as.list(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$gender))
# )
# # we need to have the same names as in 'inp'
# names(con) <- setdiff(names(inp), "db030")
# 
# # run ipu und check results
# res <- ipu(inp=inp, hid="db030", con=con, verbose=TRUE)
# 
# is <- sapply(2:(ncol(res)-1), function(x) { 
#   sum(res[,x]*res$weights)
# }) 
# data.frame(required=unlist(con), is=is)
# 
# 
# 
# 
# ######Basic not converting example
# # basic example
# require(simPop)
# inp <- as.data.frame(matrix(0, nrow=8, ncol=7))
# colnames(inp) <- c("hhid","hh1","hh2","p1","p2","p3","p4")
# inp$hhid <- 1:8
# inp$hh1[1:3] <- 1
# inp$hh2[4:8] <- 1
# inp$p1 <- c(1,1,2,1,0,1,2,1)
# inp$p2 <- c(1,0,1,0,2,1,1,1)
# inp$p3 <- c(1,1,0,2,1,0,2,0)
# inp$p4 <- c(0,0,0,0,0,0,0,0)
# con <- list(hh1=35, hh2=65, p1=91, p2=65, p3=104,p4=5)
# res <- ipu(inp=inp, hid="hhid", con=con, verbose=TRUE)
# 

######## debug mosaic
# 
# rm(list=ls())
# load("eusilcP.rdata")
# load("test0.rdata")
# 
# # original population
# index <- eusilcP$weights == 1
# eusilcP$ageG <- cut(eusilcP$age, 10)
# datO <- eusilcP[index,]
# tab1a <- table(datO$region,datO$gender, datO$hsize) 
# tab1b <- table(datO$region,datO$gender, datO$ecoStat) 
# tab1c <- table(datO$region,datO$gender, datO$ageG) 
# 
# # pop after sim-annealing
# index2 <- unlist(test0$weights) == 1
# datN <- eusilcP[index2,]
# tab2a <- table(datN$region,datN$gender, datN$hsize)
# tab2b <- table(datN$region,datN$gender, datN$ecoStat) 
# tab2c <- table(datN$region,datN$gender, datN$ageG)
# 
# # region x sex x hsize
# pdf(file="male_region_hsize.pdf", width=12, height=6)
# par(mfrow=c(1,2))
# mosaicplot(tab1a[,1,], main="males x region x hsize (orig)") 
# mosaicplot(tab2a[,1,], main="males x region x hsize (new)")
# dev.off()
# pdf(file="females_region_hsize.pdf", width=12, height=6)
# par(mfrow=c(1,2))
# mosaicplot(tab1a[,2,], main="females x region x hsize (orig)") 
# mosaicplot(tab2a[,2,], main="females x region x hsize (new)")
# dev.off()
# 
# # region x sex x ecoStat
# pdf(file="males_region_ecostat.pdf", width=12, height=6)
# par(mfrow=c(1,2))
# mosaicplot(tab1b[,1,], main="males x region x ecoStat (orig)")
# mosaicplot(tab2b[,1,], main="males x region x ecoStat (new)")
# dev.off()
# pdf(file="females_region_ecostat.pdf", width=12, height=6)
# par(mfrow=c(1,2))
# mosaicplot(tab1b[,2,], main="females x region x ecoStat (orig)")
# mosaicplot(tab2b[,2,], main="females x region x ecoStat (new)")
# dev.off()
# 
# # region x sex x ageGroups
# pdf(file="males_region_age.pdf", width=12, height=6)
# par(mfrow=c(1,2))
# mosaicplot(tab1c[,1,], main="males x region x ageG (orig)") # males X region x ecoStat
# mosaicplot(tab2c[,1,], main="males x region x ageG (new)") # males X region x ecoStat
# dev.off()
# pdf(file="females_region_age.pdf", width=12, height=6)
# par(mfrow=c(1,2))
# mosaicplot(tab1c[,2,], main="females x region x ageG (orig)") # females X region x ecoStat
# mosaicplot(tab2c[,2,], main="females x region x ageG (new)") # females X region x ecoStat
# dev.off()


############## debug simCategorical
# 
# library(parallel)
# library(simPopulation)
# library(party)
# library(LiblineaR)
# library(stringr)
# library(microbenchmark)
# library(e1071)
# 
# rm(list=ls())
# seed <- 1234
# data(eusilcS)   # load sample data
# 
# # 0.1s
# system.time({
#   eusilcP <- simStructure(eusilcS)
# })
# 
# sapply(list.files("R", full.names=TRUE), source)
# 
# 
# # non-parallel ~ 24secs
# system.time({
#   eusilcP2 <- simCategoricalOld(eusilcS, eusilcP)
# })
# 
# # conditional probabilities
# system.time({
#   eusilcP2_0 <- simCategorical(eusilcS, eusilcP, method="distribution", parallel=TRUE)
# })
# # multinomial regression from package multinom
# system.time({
#   eusilcP2_1 <- simCategorical(eusilcS, eusilcP, method="multinom", parallel=TRUE)
# })
# # recursive partitioning from package party
# #system.time({
# #  eusilcP2_2 <- simCategorical(eusilcS, eusilcP, method="ctree", parallel=FALSE)
# #})
# # naivebayes from package e1071
# system.time({
#   eusilcP2_3 <- simCategorical(eusilcS, eusilcP, method="naivebayes", parallel=TRUE)
# })
# # liblinear from package LiblineaR
# #system.time({
# #  eusilcP2_4 <- simCategorical(eusilcS, eusilcP, method="liblinear", parallel=FALSE)
# #})
# 
# # results
# tab <- rbind(
#   table(eusilcP2_0$pl030), # distr
#   table(eusilcP2_1$pl030), # multinom
#   #table(eusilcP2_2$pl030), # ctree
#   table(eusilcP2_3$pl030) # nbayes
#   #table(eusilcP2_4$pl030) # liblinear
# )
# barplot(tab, beside=TRUE, legend.text=c("distr","multinom","nbayes"), args.legend=list(x="topright"))
# 
# # benchmarking
# nrruns <- 5
# microbenchmark(
#   simCategorical(eusilcS, eusilcP, method="distribution", parallel=TRUE),
#   simCategorical(eusilcS, eusilcP, method="multinom", parallel=TRUE),
#   #simCategorical(eusilcS, eusilcP, method="ctree", parallel=TRUE),
#   simCategorical(eusilcS, eusilcP, method="naivebayes", parallel=TRUE),
#   #simCategorical(eusilcS, eusilcP, method="liblinear", parallel=TRUE),
#   times = 5
# )
# 
# 
# 
