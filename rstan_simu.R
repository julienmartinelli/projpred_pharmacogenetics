#set the working directory

#load the packages
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
source("stanTools.R")
source("functions.R")
library(reshape2)
rstan_options(auto_write = TRUE)
set.seed(11191951)

runStan <- TRUE
modelName <- "Simu"
namedir <- "/output_globlasso_corrected"
savefile <- "StanFit_globlasso.Rsave"
# modelfile <- "horseshoe.stan"
modelfile <- "global_local_laplace.stan"
projectDir <- getwd()
figDir <- tabDir <- outDir <- paste(projectDir, namedir, sep = "")
dir.create(file.path(projectDir, namedir), showWarnings = FALSE)
dir.create(file.path(figDir, "/pkparams"), showWarnings = FALSE)
dir.create(file.path(figDir, "/betaparams"), showWarnings = FALSE)
dir.create(file.path(figDir, "/cond"), showWarnings = FALSE)
dir.create(file.path(figDir, "/pred"), showWarnings = FALSE)

#load the data set
datainf <- read.csv('simu1_corrected.csv',header=TRUE)

#pk data only and removing na
xdata <- datainf %>%
  dplyr::select(ID, time, dv) %>%
  filter(!is.na(dv))

#SNP data only
tabSNP<-datainf[!duplicated(datainf$ID),grep("rs",names(datainf))]  
sum(lengths(apply(tabSNP,2,table))==1) #make sure its 0 i.e. all 134 SNPs are polymorphic

#sample size
listind<-unique(datainf$ID);(n.ind<-length(listind)) #400 patients

#the start row for each subject vector of dv
(off.data<-c(which(!duplicated(datainf$ID))))
p0 <- 20
scale_global <- 5/(dim(tabSNP)[2]-5)/sqrt(n.ind) # 5 is prior guess on relevant variables
## create data set
data <- list(
  ntot=dim(datainf)[1],
  N=n.ind,
  off_data=off.data,
  time=datainf$time,
  cObs=datainf$dv,
  dose= 200,
  nSNP=dim(tabSNP)[2], #number of tested SNPs
  tabSNP=tabSNP, #array of SNPs
  scale_global=scale_global,
  nu_global=10,
  nu_local=5,
  slab_scale=10,
  slab_df=5
)

nChains <- 8
nPost <- 1500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1500 ## Number of burn-in samples per chain after thinning
nThin <- 4

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

init <- function(){
  list(kaHat = 1,
       CLHat = 2.5,
       V1Hat = 200,
       omegaka = .1,
       omegaCL = .1,
       omegaV1 = .1,
       sigma1 = .01,
       etaka = rep(0,n.ind),
       etaCL = rep(0,n.ind),
       etaV1 = rep(0,n.ind),
       #lambda=3,
       #beta=rep(0,data["nSNP"])
       lambda_square=rep(1, data["nSNP"]),
       beta=rep(0,data["nSNP"]),
       z = rep(0, dim(tabSNP)[2])
  )
}

### Specify the variables for which you want history and density plots

# parametersToPlot <- c("kaHat", "CLHat", "V1Hat",
#                       "omegaka", "omegaCL", "omegaV1","sigma1",
#                       "beta", "tau", "logp")#, "lambda")

parametersToPlot <- c("kaHat", "CLHat", "V1Hat",
                      "omegaka", "omegaCL", "omegaV1", "sigma1",
                      "beta", "tau")


## Additional variables to monitor
otherRVs <- c("cObsPred", "cObsCond", "ka", "CL", "V1")#

parameters <- c(parametersToPlot, otherRVs)

if(runStan){
fit <- stan(#file = "one_comp_ka_CLV_LaplacePriorOnbeta.stan",
            file = modelfile,
              data = data,
              pars = parameters,
              iter = nIter,
              warmup = nBurnin,
              thin = nThin,
              init = init,
              chains = nChains,
              refresh = 50,
              control = list(adapt_delta = 0.8, stepsize = 0.01),
              cores = min(nChains, parallel::detectCores()))
  save(fit, file = savefile)
}else{
    load(file.path(projectDir, paste(savefile, sep = "")))
}
fit@stanmodel

load(savefile)
#dev.off()

pkdir <- paste(figDir, "/pkparams", sep = "")
pdf(file = file.path(pkdir, paste(modelName, "pkparams_plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

################################################################################################
### Posterior distributions of parameters

 mcmcDensity(fit, pars=c("kaHat", "CLHat", "V1Hat",
                         "omegaka", "omegaCL", "omegaV1","sigma1", "tau"),
                          byChain = TRUE)
#mcmcDensity(fit, pars=c("kaHat", "CLHat", "V1Hat",
#                        "omegaka", "omegaCL", "omegaV1","sigma1"),#, "lambda"),
#                         byChain = TRUE)
pairs(fit, pars=c("kaHat", "CLHat", "V1Hat",
                  "omegaka", "omegaCL", "omegaV1","sigma1", "tau")#, "lambda")
      )
#pairs(fit, pars=c("kaHat", "CLHat", "V1Hat",
#                  "omegaka", "omegaCL", "omegaV1","sigma1")#, "lambda")
#      )


################################################################################################
### Posterior distributions of parameters

ptable <- parameterTable(fit, parametersToPlot)

betadir <- paste(figDir, "/betaparams", sep = "")

dev.off()

pdf(file = file.path(betadir, paste(modelName, "betaparams_plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

beta_CL_estim <- parameterTable(fit, "beta")

################################################################################################
### Posterior predictive distributions
par(mar=c(7, 4, 0, 2) + 0.1)
boxplot(t(rbind(ptable[grepl("beta",row.names(ptable)),4:8]
                #,ptable[grepl("betars37",row.names(ptable)),4:8]
                ))
        #,names=c(names(tabSNP),"rs37")
        ,ylab="Genetic effect size",las=3)  
abline(h=0,col="red")

mcmcDensity(fit, parametersToPlot[grepl("beta", parametersToPlot)])
mcmcHistory(fit, parametersToPlot[grepl("beta", parametersToPlot)])
mcmcDensity(fit, parametersToPlot[grepl("beta", parametersToPlot)], byChain=TRUE)

write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
### Posterior predictive distributions

#### Prediction of future observations in the same studies, i.e., posterior predictions conditioned
#### on observed data from the same study

dev.off()

conddir <- paste(figDir, "/cond", sep = "")
pdf(file = file.path(conddir, paste(modelName, "cond_plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

pred <- as.data.frame(fit, pars = "cObsCond") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

nPerPage = 16
IDs <- sort(unique(pred$ID))
nID <- length(IDs)
nPages <- ceiling(nID / nPerPage) * 2
IDs <- data.frame(ID = IDs,
                  page = sort(rep(1:nPages, length = nID)),
                  stringsAsFactors = FALSE)
pred <- pred %>% left_join(IDs)

for(i in 1:nPages){
  xplot <- subset(pred, page == i)
  p1 <- ggplot(xplot, aes(x = time, y = dv))
  p1 <- p1 + geom_point() +
    labs(title = "individual predictions",
         x = "time (h)",
         y = "plasma concentration") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ ID)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

#### Prediction of future observations in a new study, i.e., posterior predictive distributions

dev.off()

preddir <- paste(figDir, "/pred", sep = "")
pdf(file = file.path(preddir, paste(modelName, "pkparams_plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)
  
nPerPage = 16
IDs <- sort(unique(pred$ID))
nID <- length(IDs)
nPages <- ceiling(nID / nPerPage) * 2
IDs <- data.frame(ID = IDs,
                  page = sort(rep(1:nPages, length = nID)),
                  stringsAsFactors = FALSE)
pred <- pred %>% left_join(IDs)

for(i in 1:nPages){
  xplot <- subset(pred, page == i)
  p1 <- ggplot(xplot, aes(x = time, y = dv))
  p1 <- p1 + geom_point() +
    labs(title = "individual predictions",
         x = "time (h)",
         y = "plasma concentration") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ ID)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}