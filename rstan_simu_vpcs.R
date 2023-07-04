#set the working directory
setwd("C:/Users/Julie Bertrand/ownCloud/Documents/Encadrement/Ibtissem_Rebai/code")

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

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
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
       #,lambda=40
       ,beta=rep(0,data["nSNP"])
  )
}

### Specify the variables for which you want history and density plots

parametersToPlot <- c("kaHat", "CLHat", "V1Hat",
                      "omegaka", "omegaCL", "omegaV1","sigma1"
                    #,"lambda"
                    #,"beta"
                    )#

## Additional variables to monitor
otherRVs <- c("cObsPred","ka", "CL", "V1")#  

parameters <- c(parametersToPlot, otherRVs)

fit <- stan(file = "one_comp_ka_CLV_LaplacePriorOnbeta.stan",
              data = data,
              pars = parameters,
              iter = nIter,
              warmup = nBurnin,
              thin = nThin, 
              init = init,
              chains = nChains, 
              refresh = 10,
              control = list(adapt_delta = 0.8, stepsize = 0.01),
              cores = min(nChains, parallel::detectCores()))
  save(fit, file = "StanFit_nobeta_lognormalpriorsonmus.Rsave")

  
load("StanFit_nobeta_normalpriorsonmus.Rsave")

fit
#goodness of mixing plots
plot(fit,pars=c("kaHat", "CLHat", "V1Hat",
                "omegaka", "omegaCL", "omegaV1","sigma1"
                #,"lambda"
                )
     ,plotfun="trace"
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

###############################################################################################
###boxplots of the betas

NewBoxPlot <- function(x,y){
  y <- as.character(y)
  nbox <- length(unique(y))
  f <- unique(y)
  lims <- c(0 - .5, nbox -.5)
  cents <- 0:(nbox-1)
  #labs.y <- quantile(x, c(0,.025,.1,.25,.5,.75,.9,.975,1)) # quantiles
  plot(0, xaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',xlim = lims, ylim = quantile(summm, probs = c(.025,.975), na.rm = TRUE)*3)
  for (i in 1:nbox){#i<-30
    xi <- x[y==f[i]]
    cent = cents[i]
    b <- quantile(xi, c(.025,.1,.25,.5,.75,.9,.975), na.rm = TRUE) # quantiles
    mima <- range(xi)                        # max and min
    rect(cent - .1, b[3], cent + .1, b[5])
    points(cent + c(-.1,.1), c(b[4], b[4]), type = 'l', lwd = 3)
    arrows(cent, b[5], cent, b[7],  angle = 90, length = 0.05)
    arrows(cent, b[3], cent, b[1],  angle = 90, length = 0.05)
    #points(rep(cent, 2), b[c(1,7)], pch = '-', cex = 1.5)
  }
 axis(side = 1, at = cents, labels = f, las = 2, cex.axis = 0.5)
 #axis(side = 2, at = labs.y, labels = round(labs.y,1))
}

x <- unlist(fit@sim$samples[[1]][9:(134+8)])#indexes depending on number of chains (1st []) and iterations (below)
y <- rep(names(tabSNP),each = 600)
                   
NewBoxPlot(x,y)
abline(h=0, col = "red")
################################################################################################
### Posterior predictive distributions

#### Prediction of future observations in the same studies, i.e., posterior predictions conditioned
#### on observed data from the same study

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  bind_cols(id = rep(rep(1:400), each = nPost * nChains *12)) %>%
  bind_cols(time = rep(rep(1:12,each = nPost * nChains),400)) %>%
  bind_cols(D = rep(1:(nPost * nChains),400 * 12))


pred_perc <- NULL
for(it in 1:(nPost * nChains)){# it <- 1
  pred_percit <- pred %>%
    filter(D==it) %>%
    group_by(time) %>%
    summarize(lb_D = quantile(value, probs = 0.05, na.rm = TRUE),
              median_D = quantile(value, probs = 0.5, na.rm = TRUE),
              ub_D = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
    bind_cols(D=rep(it,12))
  
  pred_perc <- rbind(pred_perc, pred_percit)
}

l5 <- pred_perc %>%
  group_by(time) %>%
  summarize(lb = quantile(lb_D, probs = 0.05, na.rm = TRUE),
            median = quantile(lb_D, probs = 0.5, na.rm = TRUE),
            ub = quantile(lb_D, probs = 0.95, na.rm = TRUE))

l50 <- pred_perc %>%
  group_by(time) %>%
  summarize(lb = quantile(median_D, probs = 0.05, na.rm = TRUE),
            median = quantile(median_D, probs = 0.5, na.rm = TRUE),
            ub = quantile(median_D, probs = 0.95, na.rm = TRUE))
  
  
l95 <- pred_perc %>%
  group_by(time) %>%
  summarize(lb = quantile(ub_D, probs = 0.05, na.rm = TRUE),
            median = quantile(ub_D, probs = 0.5, na.rm = TRUE),
            ub = quantile(ub_D, probs = 0.95, na.rm = TRUE))

obs_perc <- datainf %>%
    group_by(time) %>%
    summarize(lb = quantile(dv, probs = 0.05, na.rm = TRUE),
              median = quantile(dv, probs = 0.5, na.rm = TRUE),
              ub = quantile(dv, probs = 0.95, na.rm = TRUE))


plot(NA, NA, xlim=c(0,12), ylim=c(0,20), xlab = "Time", ylab ="Concentrations")
polygon(x = c(l5$time,rev(l5$time)), y = c(l5$lb,rev(l5$ub)), border = F, col= "lightblue")
lines(obs_perc$time,obs_perc$lb);points(obs_perc$time,obs_perc$lb)

polygon(x = c(l50$time,rev(l50$time)), y = c(l50$lb,rev(l50$ub)), border = F, col= "lightblue")
lines(obs_perc$time,obs_perc$median);points(obs_perc$time,obs_perc$median)

polygon(x = c(l95$time,rev(l95$time)), y = c(l95$lb,rev(l95$ub)), border = F, col= "lightblue")
lines(obs_perc$time,obs_perc$ub);points(obs_perc$time,obs_perc$ub)


