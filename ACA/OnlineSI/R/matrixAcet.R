## Front Matters

packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos, ...)
    require(x,character.only=TRUE)
  }
}

base <- getwd()
dataDIR <- paste(base, "Data", sep="/")
plotDIR <- paste(base, "figures", sep="/")

packages(arm)
packages(lattice)
packages(tikzDevice)

packages(rv)
packages(rstan, repos="https://cloud.r-project.org/", dependencies=T)
packages(readxl)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 200000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Importing data

### Raw data stored in separate Excel files are imported. Data from "Day1"
### were used in this analysis. After importing, the original data files
### are not touched.

## reading raw data
### 1. bovine
Chem_bovine <-
    read_excel(paste(dataDIR,
                     "D&P_cal curves_Bovine1 plasma_Dr_Song_Data.xlsx",
                     sep="/"), sheet="Acetachlor")

## extracting relevant rows and columns
bovineChem <- Chem_bovine[c(3:22,27:38),-c(3,5,7)]
names(bovineChem) <- c("Conc", "Day1", "Day3", "Day14", "Day28")

## reformating data into a data frame

bovineChem <- as.data.frame(sapply(bovineChem, as.numeric))
bovineChem$experiment <- c(rep("calib", 20), paste("QA", rep(1:3, each=4)))
bovineChem$medium <- "bovine"

### 2. human
Chem_human <-
    read_excel(paste(dataDIR,
                     "D&P_cal curves_Human plasma_Dr_Song_Data.xlsx",
                     sep="/"), sheet="Acetochlor")

humanChem <- Chem_human[c(3:22,28:39),-c(3,5,7)]
names(humanChem) <- c("Conc", "Day1", "Day3", "Day14", "Day28")

humanChem <- as.data.frame(sapply(humanChem, as.numeric))
humanChem$experiment <- c(rep("calib", 20), paste("QA", rep(1:3, each=4)))
humanChem$medium <- "human"

## 3. PBS
Chem_pbs <-
    read_excel(paste(dataDIR,
                     "D&P_cal curves_PBS_Dr_Song_Data.xlsx",
                     sep="/"), sheet="Acetochlor")

pbsChem <- Chem_pbs[c(3:17,45:53),c(1,2)]
names(pbsChem) <- c("Conc", "Day1")

pbsChem <- as.data.frame(sapply(pbsChem, as.numeric))
pbsChem$experiment <- c(rep("calib", 15), paste("QA", rep(1:3, each=3)))
pbsChem$medium <- "pbs"

## 4. Rabit
Chem_rabbit<-
    read_excel(paste(dataDIR,
                     "D&P_cal curves_Rabbit plasma_Dr_Song_Data.xlsx",
                     sep="/"), sheet="Acetochlor")

rabbitChem <- Chem_rabbit[c(5:24,31:42),-c(1,4,6,8)]
names(rabbitChem) <- c("Conc", "Day1", "Day3", "Day14", "Day28")

rabbitChem <- as.data.frame(sapply(rabbitChem, as.numeric))
rabbitChem$experiment <- c(rep("calib", 20), paste("QA", rep(1:3, each=4)))
rabbitChem$medium <- "rabbit"

## 5. rat
Chem_rat <-
    read_excel(paste(dataDIR,
                     "D&P_cal curves_Rat plasma_Dr_Song_Data.xlsx",
                     sep="/"), sheet="Acetochlor")
ratChem <- Chem_rat[c(3:22,27:38),-c(3,5,7)]
names(ratChem) <- c("Conc", "Day1", "Day3", "Day14", "Day28")

ratChem <- as.data.frame(sapply(ratChem, as.numeric))
ratChem$experiment <- c(rep("calib", 20), paste("QA", rep(1:3, each=4)))
ratChem$medium <- "rat"

## Combined data frame for multilevel model
ChemData <- data.frame(conc=c(humanChem$Conc, bovineChem$Conc, pbsChem$Conc,
                             rabbitChem$Conc, ratChem$Conc),
                      respns=c(humanChem$Day1, bovineChem$Day1, pbsChem$Day1,
                               rabbitChem$Day1, ratChem$Day1),
                      medium=c(humanChem$medium, bovineChem$medium,
                               pbsChem$medium,
                               rabbitChem$medium, ratChem$medium),
                      experiment=c(paste(humanChem$experiment, "H"),
                                   paste(bovineChem$experiment, "B"),
                                   paste(pbsChem$experiment, "P"),
                                   paste(rabbitChem$experiment, "Rb"),
                                   paste(ratChem$experiment, "Rt")))

ChemData$medium <- ordered(ChemData$medium,
                           levels=c("human","bovine","rabbit","rat","pbs"),
                           labels=c("Human","Bovine","Rabbit","Rat","PBS"))

## simple linear model
### a function to perform simple linear regression and prediction
lmPred <- function(Medium="Human"){
    tmpH <- ChemData$medium==Medium & substring(ChemData$experiment, 1,1)=="c"
    humanM <- lm(log(respns) ~ log(conc), data=ChemData, sub=tmpH)
    tmpH_pred <- ChemData$medium==Medium &
        substring(ChemData$experiment, 1,1)=="Q"

    lmPred_sim <- posterior(humanM, 4000)
    return(list(Pred=rvnorm(1,
    (log(ChemData$respns[tmpH_pred])-lmPred_sim$beta[1])/lmPred_sim$beta[2],
    lmPred_sim$sigma/lmPred_sim$beta[2]), Real=ChemData$conc[tmpH_pred]))
}

## plotting simple linear fit (not used in the paper)
par(mfrow=c(2,3))
plot(lmPred(Medium="Bovine")$Pred, ylim=log(c(0.5,100)))
abline(h=log(unique(ChemData$conc[ChemData$medium=="Bovine" &
                                  substring(ChemData$experiment, 1,1)=="Q"])))
plot(lmPred()$Pred, ylim=log(c(0.5,100)))
abline(h=log(unique(ChemData$conc[ChemData$medium=="Human" &
                                  substring(ChemData$experiment, 1,1)=="Q"])))
plot(lmPred(Medium="PBS")$Pred, ylim=log(c(0.5,100)))
abline(h=log(unique(ChemData$conc[ChemData$medium=="PBS" &
                                  substring(ChemData$experiment, 1,1)=="Q"])))
plot(lmPred(Medium="Rabbit")$Pred, ylim=log(c(0.5,100)))
abline(h=log(unique(ChemData$conc[ChemData$medium=="Rabbit" &
                                  substring(ChemData$experiment, 1,1)=="Q"])))
plot(lmPred(Medium="Rat")$Pred, ylim=log(c(0.5,100)))
abline(h=log(unique(ChemData$conc[ChemData$medium=="Rat" &
                                  substring(ChemData$experiment, 1,1)=="Q"])))

## multilevel model for fitting standard curves, a MLE alternative to BHM
## (not used in the paper)

tmp <- substring(ChemData$experiment, 1,1)=="c"

## tikz is to generate figures for LaTeX.
## dot plot of the adta (not used in the paper)
tikz(file=paste(plotDIR, "dataplot1_Acet.tex", sep="/"),
     height=3, width=4, standAlone=F)
##png("dataplot1_Acet.png", height=3*120, width=4*120)
key <- simpleKey(levels(ordered(ChemData$medium)), space = "right")
xyplot(respns~jitter(conc),
       data=ChemData, group=medium, sub=tmp,
       ylab="response",xlab="concentration",
       key=key)
dev.off()

tikz(file=paste(plotDIR, "dataplot2_Acet.tex", sep="/"),
     height=3, width=4,standAlone=T)

trellis.par.set(theme = canonical.theme("postscript", col=FALSE))

##png("dataplot2_Acet.png", height=3*120, width=4*120)
key <- simpleKey(levels(ordered(ChemData$medium)), space = "right")
##key <- list(space="right", text=levels(ordered(ChemData$medium)), points=list(pch=1:5, col="black"))
xyplot(log(respns)~jitter(log(conc)),
       data=ChemData, group=medium, sub=tmp,
       ylab="log response",xlab="log concentration",
       key=key, cex=0.5)
dev.off()

## plotting the above two figures into one figure

plotdata <- ChemData[, c("respns", "conc", "medium")]
plotdata$Log <- "Original Scale"
plotdata_log <- data.frame(respns=log(plotdata$respns), conc=log(plotdata$conc), medium=plotdata$medium)
plotdata_log$Log <- "Log Scale"

plotdata0 <- rbind(plotdata, plotdata_log)

tikz(file=paste(plotDIR, "dataplotBoth_Acet.tex", sep="/"),
     height=3.5, width=6.5,standAlone=F)
key <- simpleKey(levels(ordered(plotdata0$medium)), space = "top", columns=5, cex=0.7)
xyplot(respns~jitter(conc)|Log,
       data=plotdata0, group=medium,
       ylab="Response",xlab="Concentration ($\\mu$g/L)",
       scale=list(x="free", y="free"), key=key, aspect=1, cex=0.75)
dev.off()

## used in the paper (Figure 8)
tikz(file=paste(plotDIR, "dataplotBoth_Acet_bw.tex", sep="/"),
     height=3.5, width=6.5,standAlone=F)
trellis.par.temp <- trellis.par.get()
trellis.par.temp$add.text$cex=0.5 ## change strip text font size
trellis.par.temp$par.xlab.text$cex=0.75 ## xlab font size
trellis.par.temp$par.ylab.text$cex=0.75
trellis.par.set(theme = canonical.theme("postscript", col=FALSE),
                trellis.par.temp)
key <- simpleKey(levels(ordered(plotdata0$medium)), space = "top", columns=5, cex=0.7)
xyplot(respns~jitter(conc)|Log,
       data=plotdata0, group=medium,
       ylab="Response",xlab="Concentration ($\\mu$g/L)",
       scale=list(x="free", y="free"), key=key, aspect=1)
dev.off()

## Multilevel model (not used in the paper)
M1 <- lmer(respns~I(conc/10)+(1+I(conc/10)|medium),
           data=ChemData, sub=tmp)
summary(M1)
dotplot(ranef(M1))
coef(M1)

M2 <- lmer(log(respns)~log(conc)+(1+log(conc)|medium),
           data=ChemData, sub=tmp)
summary(M2)
dotplot(ranef(M2))
coef(M2)

## Bayesian calibration and estimation

## The stan model
BH_calib <- "
data{
  int<lower=1> N; //calibration sample size
  int<lower=1> Ngrp; //number of mediums
  int<lower=0> M; //number of observed samples
  int<lower=0> MM; //number of unique samples
  vector[N] y;
  vector[N] x;
  vector[M] y_smp;
  int ws[M];
  int<lower=1,upper=Ngrp> grp[N];
  int<lower=1,upper=Ngrp> grp_smp[M];
  int<lower=1,upper=Ngrp> grp_test[MM];
  real xlower;
  real xupper;
  real xmid;
}
parameters{
  vector[Ngrp] zb0;
  vector[Ngrp] zb1;
  real mu0;
  real mu1;
  vector<lower=xlower,upper=xupper>[MM] x_smp;
  real<lower=0> sigmaY;
  real<lower=0> sigma0;
  real<lower=0> sigma1;
  vector[Ngrp] x0mu;
  vector<lower=0>[Ngrp] x0sig;
}
transformed parameters{
  vector[N] yhat;
  vector[M] yhat_smp;
  vector[Ngrp] beta0;
  vector[Ngrp] beta1;
  beta0 = mu0+zb0*sigma0;
  beta1 = mu1+zb1*sigma1;
  for (i in 1:N){
    yhat[i] = beta0[grp[i]] + beta1[grp[i]] * (x[i]-xmid);
  }
  for (j in 1:M){
    yhat_smp[j] = beta0[grp_smp[j]] + beta1[grp_smp[j]] * (x_smp[ws[j]]-xmid);
  }
}
model{
  sigma0 ~ normal(0, 2);
  sigma1 ~ normal(0, 2);
  sigmaY ~ normal(0, 2);
  zb0 ~ std_normal();
  zb1 ~ std_normal();
  x0mu ~ normal(2.5, 5);
  x0sig ~ normal(0,5);
  for (j in 1:MM)
    x_smp[j] ~ normal(x0mu[grp_test[j]], x0sig[grp_test[j]]);
  target += normal_lpdf(y|yhat, sigmaY);
  target += normal_lpdf(y_smp| yhat_smp, sigmaY);
}
"

## compiling the stan model
fit1  <- stan_model(model_code=BH_calib)

## processing input data for the stan model
metdata <- function(data=ChemData, chains=nchains, Log=T){
    tmp <- !is.na(data$respns)
    data <- data[tmp,]
    calib <-  substring(data$experiment, 1,1)=="c"
    n <- sum(calib)
    m <- sum(!calib)
    data$grp <- as.numeric(ordered(data$medium))
    calibdata <- data[calib,]
    testdata <- data[!calib,]
    gr <- calibdata$grp
    ngr <- max(gr)
    grp0 <- testdata$grp
    wsM <- as.numeric(ordered(testdata$medium))
    ws <- as.numeric(ordered(testdata$experiment))
    oo <- order(ws)
    grptest <- wsM[oo][cumsum(table(ws[oo]))]
    ininits <- list()

    if (Log){
        indata <- list(N=n, M=m, Ngrp=ngr,MM=max(ws),
                       grp=gr, grp_smp=grp0, grp_test=grptest,  ws=ws,
                       y=log(calibdata$respns), x=log(calibdata$conc),
                       y_smp=log(testdata$respns),
                       xlower=min(log(data$conc)), xupper=max(log(data$conc)),
                       xmid=log(50))
        for (i in 1:chains)
            ininits[[i]] <- list(zb=rnorm(ngr), zb1=rnorm(ngr),
                                 mu0=rnorm(1), mu1=rnorm(1),
                                 x_smp=runif(max(ws),min(log(data$conc)),
                                             max(log(data$conc))),
                                 x0mu=runif(ngr), x0sig=runif(ngr),
                                 sigmaY=runif(1), sigma0=runif(1),
                                 sigma1=runif(1))
    } else {
        indata <- list(N=n, M=m, Ngrp=ngr,MM=max(ws),
                   grp=gr, grp_smp=grp0, grp_test=grptest,  ws=ws,
                   y=(calibdata$respns), x=(calibdata$conc),
                   y_smp=(testdata$respns),
                   xlower=min((data$conc)), xupper=max((data$conc)),
                   xmid=(50))
        for (i in 1:chains)
            ininits[[i]] <- list(zb=rnorm(ngr), zb1=rnorm(ngr),
                                 mu0=rnorm(1), mu1=rnorm(1),
                                 x_smp=runif(max(ws),min((data$conc)),
                                             max((data$conc))),
                                 x0mu=runif(ngr), x0sig=runif(ngr),
                                 sigmaY=runif(1), sigma0=runif(1),
                                 sigma1=runif(1))
    }
    para=c("beta0","beta1","x_smp","sigma0","sigma1","sigmaY","mu0","mu1")
    return(list(data=indata, inits=ininits, pars=para, nchains=chains))
}

### running the stan model
input.to.stan <- metdata()

fit2keep <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)#,
save(fit2keep, file="AcetBH_reults_log.RData")
print(fit2keep)

## load("AcetBH_reults_log.RData")

## processing results

BHest <- rvsims(extract(fit2keep, pars="x_smp")[[1]])

tmp <- substring(ChemData$experiment, 1, 2)=="QA"
realC <- ChemData$conc[tmp]
smpID_C <- ordered(ChemData$experiment[tmp])
smpID <- as.numeric(smpID_C)

## QA 1:
qa1 <- summary(qa1rv <- c(lmPred(Medium="Bovine")$Pred[1:4],
                 lmPred(Medium="Human")$Pred[1:4],
                 lmPred(Medium="PBS")$Pred[1:3],
                 lmPred(Medium="Rabbit")$Pred[1:4],
                 lmPred(Medium="Rat")$Pred[1:4]))
ind1 <- rep(1:5, c(4,4,3,4,4))
qa1real <- c(lmPred(Medium="Bovine")$Real[1:4],
                 lmPred(Medium="Human")$Real[1:4],
                 lmPred(Medium="PBS")$Real[1:3],
                 lmPred(Medium="Rabbit")$Real[1:4],
                 lmPred(Medium="Rat")$Real[1:4])
## QA 2
qa2 <- summary(qa2rv <- c(lmPred(Medium="Bovine")$Pred[5:8],
                 lmPred(Medium="Human")$Pred[5:8],
                 lmPred(Medium="PBS")$Pred[4:6],
                 lmPred(Medium="Rabbit")$Pred[5:8],
                 lmPred(Medium="Rat")$Pred[5:8]))
ind2 <- rep(6:10, c(4,4,3,4,4))
qa2real <- c(lmPred(Medium="Bovine")$Real[5:8],
                 lmPred(Medium="Human")$Real[5:8],
                 lmPred(Medium="PBS")$Real[4:6],
                 lmPred(Medium="Rabbit")$Real[5:8],
                 lmPred(Medium="Rat")$Real[5:8])
## QA 3
qa3 <- summary(qa3rv <- c(lmPred(Medium="Bovine")$Pred[9:12],
                 lmPred(Medium="Human")$Pred[9:12],
                 lmPred(Medium="PBS")$Pred[7:9],
                 lmPred(Medium="Rabbit")$Pred[9:12],
                 lmPred(Medium="Rat")$Pred[9:12]))
ind3 <- rep(11:15, c(4,4,3,4,4))
qa3real <- c(lmPred(Medium="Bovine")$Real[9:12],
                 lmPred(Medium="Human")$Real[9:12],
                 lmPred(Medium="PBS")$Real[7:9],
                 lmPred(Medium="Rabbit")$Real[9:12],
                 lmPred(Medium="Rat")$Real[9:12])


BHest_sum <- summary(BHest)

tikz(file=paste(plotDIR, "comparisons.tex", sep="/"),
     height=4, width=5, standAlone=T) ## not used in the paper
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1,tck=0.01)
plot(smpID, log(realC), pch=16, col="red", cex=0.5, ylim=c(0.5,4.5),
     ylab="log concentration", xlab="test sample index")
segments(x0=(1:15)+0.1, x1=(1:15)+0.1,
         y0=BHest_sum$"2.5%", y1=BHest_sum$"97.5%", col="blue")
segments(x0=(1:15)+0.1, x1=(1:15)+0.1,
         y0=BHest_sum$"25%", y1=BHest_sum$"75%", lwd=3, col="blue")
points((1:15)+0.1, BHest_sum$mean, cex=0.5)

## qa1
segments(x0=(1:5)-0.1, x1=(1:5)-0.1,
         y0=qa1[c(1, 5, 9, 12, 16),]$"2.5%",
         y1=qa1[c(1, 5, 9, 12, 16),]$"97.5%",col="gray")
segments(x0=(1:5)-0.1, x1=(1:5)-0.1,
         y0=qa1[c(1, 5, 9, 12, 16),]$"25%",
         y1=qa1[c(1, 5, 9, 12, 16),]$"75%", lwd=3, col="gray")
points((1:5)-0.1, qa1[c(1, 5, 9, 12, 16),]$mean, cex=0.5)

segments(x0=(1:5)-0.2, x1=(1:5)-0.2,
         y0=qa1[c(1, 5, 9, 12, 16)+1,]$"2.5%",
         y1=qa1[c(1, 5, 9, 12, 16)+1,]$"97.5%", col="gray")
segments(x0=(1:5)-0.2, x1=(1:5)-0.2,
         y0=qa1[c(1, 5, 9, 12, 16)+1,]$"25%",
         y1=qa1[c(1, 5, 9, 12, 16)+1,]$"75%", lwd=3, col="gray")
points((1:5)-0.2, qa1[c(1, 5, 9, 12, 16)+1,]$mean, cex=0.5)

segments(x0=(1:5)-0.3, x1=(1:5)-0.3,
         y0=qa1[c(1, 5, 9, 12, 16)+2,]$"2.5%",
         y1=qa1[c(1, 5, 9, 12, 16)+2,]$"97.5%", col="gray")
segments(x0=(1:5)-0.3, x1=(1:5)-0.3,
         y0=qa1[c(1, 5, 9, 12, 16)+2,]$"25%",
         y1=qa1[c(1, 5, 9, 12, 16)+2,]$"75%", lwd=3, col="gray")
points((1:5)-0.3, qa1[c(1, 5, 9, 12, 16)+2,]$mean, cex=0.5)

segments(x0=c(1,2,4,5)-0.4, x1=c(1,2,4,5)-0.4,
         y0=qa1[c(1, 5, 12, 16)+3,]$"2.5%",
         y1=qa1[c(1, 5, 12, 16)+3,]$"97.5%", col="gray")
segments(x0=c(1,2,4,5)-0.4, x1=c(1,2,4,5)-0.4,
         y0=qa1[c(1, 5, 12, 16)+3,]$"25%",
         y1=qa1[c(1, 5, 12, 16)+3,]$"75%", lwd=3, col="gray")
points(c(1,2,4,5)-0.4, qa1[c(1, 5, 12, 16)+3,]$mean, cex=0.5)

## qa2
segments(x0=(6:10)-0.1, x1=(6:10)-0.1,
         y0=qa2[c(1, 5, 9, 12, 16),]$"2.5%",
         y1=qa2[c(1, 5, 9, 12, 16),]$"97.5%", col="gray")
segments(x0=(6:10)-0.1, x1=(6:10)-0.1,
         y0=qa2[c(1, 5, 9, 12, 16),]$"25%",
         y1=qa2[c(1, 5, 9, 12, 16),]$"75%", lwd=3, col="gray")
points((6:10)-0.1, qa2[c(1, 5, 9, 12, 16),]$mean, cex=0.5)

segments(x0=(6:10)-0.2, x1=(6:10)-0.2,
         y0=qa2[c(1, 5, 9, 12, 16)+1,]$"2.5%",
         y1=qa2[c(1, 5, 9, 12, 16)+1,]$"97.5%", col="gray")
segments(x0=(6:10)-0.2, x1=(6:10)-0.2,
         y0=qa2[c(1, 5, 9, 12, 16)+1,]$"25%",
         y1=qa2[c(1, 5, 9, 12, 16)+1,]$"75%", lwd=3, col="gray")
points((6:10)-0.2, qa2[c(1, 5, 9, 12, 16)+1,]$mean, cex=0.5)

segments(x0=(6:10)-0.3, x1=(6:10)-0.3,
         y0=qa2[c(1, 5, 9, 12, 16)+2,]$"2.5%",
         y1=qa2[c(1, 5, 9, 12, 16)+2,]$"97.5%", col="gray")
segments(x0=(6:10)-0.3, x1=(6:10)-0.3,
         y0=qa2[c(1, 5, 9, 12, 16)+2,]$"25%",
         y1=qa2[c(1, 5, 9, 12, 16)+2,]$"75%", lwd=3, col="gray")
points((6:10)-0.3, qa2[c(1, 5, 9, 12, 16)+2,]$mean, cex=0.5)

segments(x0=c(6,7,9,10)-0.4, x1=c(6,7,9,10)-0.4,
         y0=qa2[c(1, 5, 12, 16)+3,]$"2.5%",
         y1=qa2[c(1, 5, 12, 16)+3,]$"97.5%", col="gray")
segments(x0=c(6,7,9,10)-0.4, x1=c(6,7,9,10)-0.4,
         y0=qa2[c(1, 5, 12, 16)+3,]$"25%",
         y1=qa2[c(1, 5, 12, 16)+3,]$"75%", lwd=3, col="gray")
points(c(6,7,9,10)-0.4, qa2[c(1, 5, 12, 16)+3,]$mean, cex=0.5)

## qa3
segments(x0=(11:15)-0.1, x1=(11:15)-0.1,
         y0=qa3[c(1, 5, 9, 12, 16),]$"2.5%",
         y1=qa3[c(1, 5, 9, 12, 16),]$"97.5%", col="gray")
segments(x0=(11:15)-0.1, x1=(11:15)-0.1,
         y0=qa3[c(1, 5, 9, 12, 16),]$"25%",
         y1=qa3[c(1, 5, 9, 12, 16),]$"75%", lwd=3, col="gray")
points((11:15)-0.1, qa3[c(1, 5, 9, 12, 16),]$mean, cex=0.5)

segments(x0=(11:15)-0.2, x1=(11:15)-0.2,
         y0=qa3[c(1, 5, 9, 12, 16)+1,]$"2.5%",
         y1=qa3[c(1, 5, 9, 12, 16)+1,]$"97.5%", col="gray")
segments(x0=(11:15)-0.2, x1=(11:15)-0.2,
         y0=qa3[c(1, 5, 9, 12, 16)+1,]$"25%",
         y1=qa3[c(1, 5, 9, 12, 16)+1,]$"75%", lwd=3, col="gray")
points((11:15)-0.2, qa3[c(1, 5, 9, 12, 16)+1,]$mean, cex=0.5)

segments(x0=(11:15)-0.3, x1=(11:15)-0.3,
         y0=qa3[c(1, 5, 9, 12, 16)+2,]$"2.5%",
         y1=qa3[c(1, 5, 9, 12, 16)+2,]$"97.5%", col="gray")
segments(x0=(11:15)-0.3, x1=(11:15)-0.3,
         y0=qa3[c(1, 5, 9, 12, 16)+2,]$"25%",
         y1=qa3[c(1, 5, 9, 12, 16)+2,]$"75%", lwd=3, col="gray")
points((11:15)-0.3, qa3[c(1, 5, 9, 12, 16)+2,]$mean, cex=0.5)

segments(x0=c(11,12,14,15)-0.4, x1=c(11,12,14,15)-0.4,
         y0=qa3[c(1, 5, 12, 16)+3,]$"2.5%",
         y1=qa3[c(1, 5, 12, 16)+3,]$"97.5%", col="gray")
segments(x0=c(11,12,14,15)-0.4, x1=c(11,12,14,15)-0.4,
         y0=qa3[c(1, 5, 12, 16)+3,]$"25%",
         y1=qa3[c(1, 5, 12, 16)+3,]$"75%", lwd=3, col="gray")
points(c(11,12,14,15)-0.4, qa3[c(1, 5, 12, 16)+3,]$mean, cex=0.5)
dev.off()

## calculating bias
### BHM

BHbias <- rv(0)
for (i in 1:15)
    BHbias <- c(BHbias, abs(BHest[i]-mean(log(realC)[smpID==i])))
BHbiasM <- c(mean(BHbias[1:5], na.rm=T), mean(BHbias[6:10], na.rm=T),
             mean(BHbias[11:15], na.rm=T))
### Linear model
lmbias <- c(mean(abs(qa1rv-log(qa1real)), na.rm=T),
            mean(abs(qa2rv-log(qa2real)), na.rm=T),
            mean(abs(qa3rv-log(qa3real)), na.rm=T))

### Figure 7
tikz(file=paste(plotDIR, "biasAcet_comp.tex", sep="/"),
     height=3, width=3, standAlone=F)
##png(file="bias_comp.png", height=3*120, width=3*120)
par(mar=c(3,1,1,0.5), mgp=c(1.75, 0.15,0), las=1, tck=-0.01)
mlplot(rvmatrix(c(BHbiasM,lmbias), ncol=2), xlab="bias (\\%)",
       axes=F, xlim=log(c(1.045,1.175)), ylim=c(1,3))
text(x=rep(log(1.045), 3), y=1:3, text=1:3)
axis(1, at=log(seq(5,20,5)/100 +1),labels=seq(5,20,5)) 

dev.off()

## multivariate multiple regression for matrix effect
betas <- extract(fit2keep, pars=c("beta0","beta1"))

betas_rv <- rvsims(as.matrix(as.data.frame(betas)))
write.csv(summary(betas_rv), file=paste(dataDIR, "betas_Acet.csv", sep="/"))

## regression (not used in the paper)
bov2hum <- lm(cbind(betas$beta0[,2], betas$beta1[,2]) ~
                  betas$beta0[,1]+betas$beta1[,1])


summary(bov2hum)

png(file="post_Bov_Hum.png", height=4*120, width=4*120)
par(mfrow=c(1,2), mar=c(3,3,1,0.5), mgp=c(1.75, 0.125,0), las=1, tck=-0.01)
pairs(betas$beta0[,2]~ betas$beta1[,2] +
          betas$beta0[,1]+betas$beta1[,1],
      labels=c("Int.bov", "Slp.Bov", "Int.Hum", "Slp.Hum"))
dev.off()

## ANOVA (not used in the paper)
b0 <- data.frame(b0=as.vector(betas$beta0),
                 medium= rep(1:5, each=dim(betas$beta0)[1]))
b1 <- data.frame(b1=as.vector(betas$beta1),
                 medium= rep(1:5, each=dim(betas$beta1)[1]))

b0aov <- aov(b0~factor(medium), data=b0)
b1aov <- aov(b1~factor(medium), data=b1)
summary(b0aov)
plot(TukeyHSD(b0aov))

### Figure 9 (version 1, not used)
tikz(file=paste(plotDIR, "matrix_correct.tex", sep="/"),
     height=3, width=4.75, standAlone=F)
##png(file="matrix_correctAcet.png", height=2.25*120, width=5*125)
par(mfrow=c(1,2), mar=c(3,3,1,0.5), mgp=c(1.75, 0.125,0), las=1, tck=-0.01)
boxplot(b0~medium, data=b0, axes=F, ylab="Intercept")
axis(1, at=1:5, labels=levels(ordered(ChemData$medium)))
axis(2)
box()
boxplot(b1~medium, data=b1, axes=F, ylab="Slope")
axis(1, at=1:5, labels=levels(ordered(ChemData$medium)))
axis(2)
box()
dev.off()

beta0 <- betas$beta0
dimnames(beta0) [[2]]  <- levels(ordered(ChemData$medium))
beta0 <- rvsims(beta0)

beta1 <- betas$beta1
dimnames(beta1) [[2]]  <- levels(ordered(ChemData$medium))
beta1 <- rvsims(beta1)

## Figure 9 
tikz(file=paste(plotDIR, "matrix_correct2.tex", sep="/"),
     height=3, width=4.75, standAlone=F)
##png(file="matrix_correctAcet.png", height=2.25*120, width=5*125)
par(mfrow=c(1,2), mar=c(3,3,1,0.5), mgp=c(1.75, 0.125,0), las=1, tck=-0.01)
mlplot(beta0, xlab="intercept")
mlplot(beta1, xlab="slope", xlim=c(0.9,1.2))
dev.off()
