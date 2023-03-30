#set the working directory
# setwd("C:/Users/Julie Bertrand/ownCloud/Documents/Encadrement/Ibtissem_Rebai/code")

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

modelName <- "Simu"

scriptDir <- projectDir <-dataDir <- modelDir <- getwd()
figDir <- tabDir <- outDir <- paste(projectDir, "/output", sep = "")
toolsDir <- getwd()

#load the data set
datainf <- read.csv('simu1.csv',header=TRUE)

#pk data only and removing na
xdata <- datainf %>%
  dplyr::select(ID, time, dv) %>%
  filter(!is.na(dv))

#SNP data only
tabSNP<-datainf[!duplicated(datainf$ID),grep("rs",names(datainf))]  
sum(lengths(apply(tabSNP,2,table))==1)#make sure its 0 i.e. all 134 SNPs are polymorphic

#sample size
listind<-unique(datainf$ID);(n.ind<-length(listind)) #400 patients

#the start row for each subject vector of dv
(off.data<-c(which(!duplicated(datainf$ID)))) 

## create data set
data <- list(
  ntot=dim(datainf)[1],
  N=n.ind,
  off_data=off.data,
  time=datainf$time,
  cObs=datainf$dv,
  dose= 200,
  nSNP=dim(tabSNP)[2], #number of tested SNPs
  tabSNP=tabSNP #array of SNPs
)

nChains <- 8
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

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
       etaV1 = rep(0,n.ind)
       ,lambda=40,  
       beta=rep(0,data["nSNP"])
  )
}

### Specify the variables for which you want history and density plots

parametersToPlot <- c("kaHat", "CLHat", "V1Hat",
                      "omegaka", "omegaCL", "omegaV1","sigma1"
                    ,"lambda","beta")#

## Additional variables to monitor
otherRVs <- c("cObsPred","cObsCond","ka", "CL", "V1")#  

parameters <- c(parametersToPlot, otherRVs)

fit <- stan(file = "one_comp_ka_CLV_LaplacePriorOnbeta.stan",
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
  save(fit, file = "StanFit_Laplace_lambda.Rsave")

  
load("StanFit_Laplace_lambda.Rsave")

fit
#goodness of mixing plots
plot(fit,pars=c("kaHat", "CLHat", "V1Hat",
                "omegaka", "omegaCL", "omegaV1","sigma1"
                ,"lambda"),plotfun="trace"
)

################################################################################################
### Posterior distributions of parameters

mcmcDensity(fit, pars=c("kaHat", "CLHat", "V1Hat",
                        "omegaka", "omegaCL", "omegaV1","sigma1"
                        ,"lambda"), byChain = TRUE)

pairs(fit, pars=c("kaHat", "CLHat", "V1Hat",
                  "omegaka", "omegaCL", "omegaV1","sigma1"
                  ,"lambda")
      )

ptable <- parameterTable(fit, parametersToPlot)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
### Posterior predictive distributions

#### Prediction of future observations in the same studies, i.e., posterior predictions conditioned
#### on observed data from the same study

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
nPages <- ceiling(nID / nPerPage)
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

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)
