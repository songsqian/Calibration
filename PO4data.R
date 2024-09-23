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

packages(rv)
packages(rstan)
packages(tidyverse)
packages(tikzDevice)
packages(readxl)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

DUWC <- read_excel(paste(dataDIR, "Statistics for typical Ortho-P calibration_FINAL.xlsx", sep="/"),
                   col_names=T, skip=5)[,1:5]

names(DUWC)  <- c("Test", "Analyte", "Type", "Conc", "Abs")

DUWC <- DUWC[!is.na(DUWC$Test), ]
tmp <- substring(DUWC$Test, 1, 5)!="Track"
DUWC <- DUWC[tmp,]
DUWC$Conc[DUWC$Conc=="SAMPLE"] <- NA
DUWC$Conc <- as.numeric(DUWC$Conc)
DUWC$Abs <- as.numeric(DUWC$Abs)

Tests <- sort(unique(DUWC$Test))

pdf(file="log-log_or_not.pdf", height=5, width=3)
par(mfrow=c(2,1), mar=c(3,3,0.25,0.25), mgp=c(1.25,0.125,0), tck=-0.01)
plot(Abs~Conc, data=DUWC, ylim=c(0,1))
plot(Abs~Conc, data=DUWC, log="xy")
dev.off()

### using rv::posterior
LmMLEBCrv_cntrl <- rv(0)
for (i in 1:length(Tests)){
    temp <- DUWC[DUWC$Test==Tests[i], ]
    M <- rv::posterior(lm(Abs ~ Conc, data=temp, sub=Type=="CAL STD"), 4000)
    y0 <- temp$Abs[temp$Type=="INT CHK"]
    x0 <- (y0 - M$beta[1]-rvnorm(1,0,M$sigma))/M$beta[2]
    LmMLEBCrv_cntrl <- c(LmMLEBCrv_cntrl, x0)
}

LmMLEBCrvlog_cntrl <- rv(0)
for (i in 1:length(Tests)){
    temp <- DUWC[DUWC$Test==Tests[i], ]# & DUWC$Conc>0, ]
    M <- rv::posterior(lm(log(Abs) ~ log(Conc), data=temp, sub=Conc>0), 4000)
    y0 <- temp$Abs[temp$Type=="INT CHK"]
    x0 <- exp((log(y0) - M$beta[1]-rvnorm(1,0,M$sigma))/M$beta[2])
    LmMLEBCrvlog_cntrl <- c(LmMLEBCrvlog_cntrl, x0)
}


### MLE (Bayes fitting with and without samples)
LmMLE <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // response data
  real x[N]; // predictor data
  real y0[M]; // test data response
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
  vector[M] logx0;
}
transformed parameters{
  real mu[N];
  real mu0[M];
  vector[M] x0;

  x0 = exp(logx0);

  for (i in 1:N){
    mu[i] = beta0+beta1*(x[i]);
  }
  for (i in 1:M){
    mu0[i] = beta0+beta1*(x0[i]);
  }
}
model {
  logx0 ~ normal(0,10);
  beta0~ normal(0,5);
  beta1 ~ normal(0,5);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"

LmMLE_log <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // response data (log)
  real x[N]; // predictor data (log)
  real y0[M]; // test data response (log)
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
  vector[M] x0;
}
transformed parameters{
  real mu[N];
  real mu0[M];
//  vector[M] x0;
//  x0 = exp(logx0);


  for (i in 1:N){
    mu[i] = beta0+beta1*(x[i]);
  }
  for (i in 1:M){
    mu0[i] = beta0+beta1*(x0[i]);
  }
}
model {
  x0 ~ normal(0,10);
  beta0~normal(0,5);
  beta1~normal(0,5);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"

stan.fit1 <- stan_model(model_code=LmMLE)
stan.fit1Log <- stan_model(model_code=LmMLE_log)

stan.in1 <- function(data = DUWC, chains=nchains, stdOnly=F, Log=F){
    n <- dim(data)[1]
    temp  <- data$Type=="CAL STD"
    y <- data$Abs[temp]
    x <- data$Conc[temp]
    N <- sum(temp)
    if (stdOnly){
        M <- sum(data$Type=="INT CHK")
        y0 <- data$Abs[data$Type=="INT CHK"]
        cntrl <- NULL
    } else {
        M <- n-N
        y0 <- data$Abs[!temp]
        cntrl <- data$Type[!temp]=="INT CHK"
    }
    if (Log) {
        N <- sum(x>0)
        y0 <- log(y0)
        y <- log(y[x>0])
        x <- log(x[x>0])
    }
    stan.dat <- list(N=N, M=M, y=y, x=x-mean(x), y0=array(y0))
    inits <- list()
    for (i in 1:chains)
        if(Log) inits[[i]] <- list(beta0 = rnorm(1), beta1=runif(1),
                                   x0 = array(rnorm(M)), sigma = runif(1))
        else inits[[i]] <- list(beta0 = rnorm(1), beta1=runif(1),
                                logx0 = array(rnorm(M)), sigma = runif(1))
    parameters <- c("beta0","beta1","sigma","x0")
    return(list(para=parameters, data=stan.dat, inits=inits,
              n.chains=chains, ControlSample=cntrl, NoSmpl=M, x_bar=mean(x)))
}

## Estimation uncertainty when there are multiple analytes to estimate 
input.to.stan <- stan.in1(data=DUWC[DUWC$Test==Tests[1], ], stdOnly=T)
x_bar0 <- input.to.stan$x_bar
fit2keep <- sampling(stan.fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)
LmMLE_OneX0 <- rstan::extract(fit2keep)

input.to.stan <- stan.in1(data=DUWC[DUWC$Test==Tests[1], ], stdOnly=F)
x_bar1 <- input.to.stan$x_bar
fit2keep <- sampling(stan.fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)
LmMLE_moreX0 <- rstan::extract(fit2keep)
cntl <- input.to.stan$ControlSample
OneX0 <- rvsims(as.matrix(as.data.frame(LmMLE_OneX0)))
OneX0 <- OneX0[names(OneX0)=="x0"]

MoreX0 <- rvsims(as.matrix(as.data.frame(LmMLE_moreX0)))
MoreX0 <- MoreX0[substring(names(MoreX0),1,3)=="x0."][cntl]
X0comp <-c(OneX0+x_bar0, MoreX0+x_bar1)
names(X0comp) <- c("One", "More")
X0devia <- abs(X0comp-50)
mlplot(X0comp)
mlplot(X0devia)
## end of compare
## Estimation uncertainty increases when more unknown samples are included


LmMLE_stan_out <- list()
LmMLE_Cntrl <- list()
x_bars_mle <- numeric()

for (i in 1:length(Tests)){
    input.to.stan <- stan.in1(data=DUWC[DUWC$Test==Tests[i], ])
    LmMLE_Cntrl[[i]] <- input.to.stan$ControlSample
    x_bars_mle[i] <- input.to.stan$x_bar
}
for (i in 1:length(Tests)){
    input.to.stan <- stan.in1(data=DUWC[DUWC$Test==Tests[i], ])
    x_bars_mle[i] <- input.to.stan$x_bar

    print(paste("Test = ", i, " of ", length(Tests), sep=""))
    fit2keep <- sampling(stan.fit1, data=input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    LmMLE_Cntrl[[i]] <- input.to.stan$ControlSample
    LmMLE_stan_out[[i]] <- rstan::extract(fit2keep)
}

save(x_bars_mle, LmMLE_stan_out, file="MLEoutputLM.RData")
## load("MLEoutputLM.RData")
##pairs(fit2keep, pars=c("beta0","beta1","sigma"))

LmMLEcntrl <- rv(0)
for (i in 1:length(Tests)){
    temp <- rvsims(as.matrix(as.data.frame(LmMLE_stan_out[[i]])))
    tmp <- substring(names(temp), 1, 3)=="x0."
    LmMLEcntrl <- c(LmMLEcntrl, temp[tmp][LmMLE_Cntrl[[i]]]+x_bars_mle[i])
}


### BHM within (use centered x)
LmBHM1 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // response data
  real x[N]; // predictor data
  real y0[M]; // test data response
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
  vector[M] logx0;
  real mux0;
  real<lower=0> sigx0;
}
transformed parameters{
  real mu[N];
  real mu0[M];
  vector[M] x0;

  x0 = exp(logx0);
  for (i in 1:N){
    mu[i] = beta0+beta1*(x[i]);
  }
  for (i in 1:M){
    mu0[i] = beta0+beta1*(x0[i]);
  }
}
model {
  mux0 ~ normal(0,2.5);
  sigx0 ~ normal(0,2.5);
  logx0 ~ normal(mux0,sigx0);
  beta0~normal(0,5);
  beta1~normal(0,5);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"
stan.fit2 <- stan_model(model_code=LmBHM1)
stan.in2 <- function(data = DUWC, chains=nchains){
    n <- dim(data)[1]
    temp  <- data$Type=="CAL STD"
    y <- data$Abs[temp]
    x <- data$Conc[temp]
    N <- sum(temp)
    M <- n-sum(temp)
    y0 <- data$Abs[!temp]
    cntrl <- data$Type[!temp]=="INT CHK"
    stan.dat <- list(N=N, M=M, y=y, x=x-mean(x), y0=y0)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(beta0 = runif(1), beta1=runif(1),
                           x0 = runif(M), sigma = runif(1),
                           mux0=runif(1), sigx0=runif(1))
    parameters <- c("beta0","beta1","sigma","x0","mux0","sigx0")
    return(list(para=parameters, data=stan.dat, inits=inits,
              n.chains=chains, ControlSample=cntrl, x_bar=mean(x)))
}

LmBHM1_stan_out <- list()
x_bars_bhm1 <- numeric()
for (i in 1:length(Tests)){
    input.to.stan <- stan.in2(data=DUWC[DUWC$Test==Tests[i], ])
    x_bars_bhm1[i] <- input.to.stan$x_bar
    print(paste("Test = ", i, " of ", length(Tests), sep=""))

    fit2keep <- sampling(stan.fit2, data=input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    LmBHM1_stan_out[[i]] <- rstan::extract(fit2keep)
}

save(x_bars_bhm1, LmBHM1_stan_out, file="BHM1outputLM.RData")
## load("BHM1outputLM.RData")
#pairs(fit2keep, pars=c("beta0","beta1","sigma", "mux0","sigx0"))
LmBHM1cntrl <- rv(0)
for (i in 1:length(Tests)){
    temp <- rvsims(as.matrix(as.data.frame(LmBHM1_stan_out[[i]])))
    tmp <- substring(names(temp), 1, 3)=="x0."
    LmBHM1cntrl <- c(LmBHM1cntrl,
                     temp[tmp][ LmMLE_Cntrl[[i]]]+x_bars_bhm1[i])
}

### BHM within and across (use centered x)

LmBHM2 <- "
data {
  int<lower=0> N; // sample size
  int<lower=0> K; // number of test kits
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // std solution response data
  real x[N]; // std solution conc
  real y0[M]; // test data response
  int<lower=1,upper=K> test1[N];
  int<lower=1,upper=K> test2[M];
}
parameters {
  vector[K] zbeta0;
  vector[K] zbeta1;
  vector<lower=0>[K] sigma;
  vector[M] logx0;
  vector<lower=0>[K] sigx0;
  vector[K] mux0;
  real mu_beta0;
  real<lower=0> sig_beta0;
  real mu_beta1;
  real<lower=0> sig_beta1;
}
transformed parameters{
  vector[K] beta0;
  vector[K] beta1;
  vector[N] mu;
  vector[M] mu0;
  vector<lower=0>[N] sig;
  vector<lower=0>[M] sig0;
  vector<lower=0>[M] x0;

  x0 = exp(logx0);
  beta0 = mu_beta0+sig_beta0*zbeta0;
  beta1 = mu_beta1+sig_beta1*zbeta1;
  for (i in 1:N){
    mu[i] = beta0[test1[i]] + beta1[test1[i]]*x[i];
    sig[i] = sigma[test1[i]];
  }
  for (i in 1:M){
    mu0[i] = beta0[test2[i]] + beta1[test2[i]]*x0[i];
    sig0[i] = sigma[test2[i]];
  }
}
model {
  mux0 ~ normal(0, 2.5);
  sigx0 ~ normal(0, 2.5);
  mu_beta0 ~ normal(0, 2.5);
  sig_beta0 ~ normal(0, 2.5);
  mu_beta1 ~ normal(0, 2.5);
  sig_beta1 ~ normal(0, 2.5);

  for (j in 1:M)
    logx0[j] ~ normal(mux0[test2[j]], sigx0[test2[j]]);
  zbeta0~std_normal();
  zbeta1~std_normal();
  target += normal_lpdf(y| mu, sig);
  target += normal_lpdf(y0| mu0, sig0);
}
"

stan.fit3 <- stan_model(model_code=LmBHM2)

stan.in3 <- function(data = DUWC, chains=nchains){ ## not centering x
    n <- dim(data)[1]
    temp  <- data$Type=="CAL STD"
    y <- data$Abs[temp]
    x <- data$Conc[temp]
    N <- sum(temp)
    M <- n-sum(temp)
    y0 <- data$Abs[!temp]
    test1 <- as.numeric(ordered(data$Test[temp]))
    test2 <- as.numeric(ordered(data$Test[!temp]))
    K <- max(test1)
    cntrl <- data$Type[!temp]=="INT CHK"
    stan.dat <- list(N=N, M=M, K=K, y=y, x=x-mean(x), y0=y0,
                     test1=test1, test2=test2)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(zbeta0 = rnorm(K), zbeta1=rnorm(K),
                           logx0 = rnorm(M), sigma = runif(K),
                           mux0=rnorm(K), sigx0=runif(K),
                           mu_beta0=runif(1), sig_beta0=runif(1),
                           mu_beta1=runif(1), sig_beta1=runif(1))
    parameters <- c("beta0","beta1","sigma","x0","mux0","sigx0",
                    "mu_beta0", "sig_beta0", "mu_beta1", "sig_beta1")
    return(list(para=parameters, data=stan.dat, inits=inits,
              n.chains=chains, ControlSample=cntrl, x_bar=mean(x)))
}

input.to.stan <- stan.in3(data=DUWC)
x_bar_bhm2 <- input.to.stan$x_bar
fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)

test_stanBHM2 <- rstan::extract(fit2keep)

save(x_bar_bhm2, test_stanBHM2, file="BHM2outputLM.RData")
## load("BHM2outputLM.RData")

temp <- rvsims(as.matrix(as.data.frame(test_stanBHM2)))
LmBHM2cntrl <- temp[substring(names(temp), 1, 3)=="x0."][ input.to.stan $ControlSample] + x_bar_bhm2

nn <- dim(ChkCIs)[1]
par(mfrow=c(2,2))
mlplot((LmBHM2cntrl), xlim=c(0,80), xlab="BHM2")
abline(v=50)
mlplot(LmBHM1cntrl, xlim=c(0, 80), xlab="BHM1")
abline(v=50)
mlplot(rev(LmMLEcntrl), xlim=c(0, 80), xlab="Bayes")
abline(v=50)
plot(c(0,1), c(0,1), xlim=c(0,80), ylim=c(1,nn), xlab="Invs", type="n")
segments(x0=ChkCIs[,4],x1=ChkCIs[,5], y0=1:nn, y1=1:nn)
segments(x0=ChkCIs[,2],x1=ChkCIs[,3], y0=1:nn, y1=1:nn, lwd=3)
points(x=ChkCIs[,1], 1:nn)
abline(v=50)

abserr <- c(mean(abs(LmBHM2cntrl-50)),
            mean(abs(LmBHM1cntrl-50)),
            mean(abs(LmMLEcntrl-50)),
            mean(abs(LmMLEBCrv_cntrl-50)))
names(abserr) <- c("$BHM_2$","$BHM_1$","Bayes","Invs")

tikz(file=paste(plotDIR, "bias_PO4.tex", sep="/"),
     height=2.5, width=2.5, standAlone=F)
par(mar=c(3, 3, 2, 1.25), mgp=c(1.25,0.25,0), tck=-0.01)
mlplot(abserr[-4], xlab="mean absolute error",
       xlim=c(1, 20), cex=0.75, cex.axis=0.75,
       mgp=c(1.25,0.25,0), tck=-0.01)
dev.off()
save(abserr, file="PO4abserr.RData")
