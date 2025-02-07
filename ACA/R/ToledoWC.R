## Toledo Water Crisis Example

### Frontmatters
packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos, ...)
    require(x,character.only=TRUE)
  }
}

### From Qian (2016) -- Monte Carlo simulation of nonlinear models
sim.nls <- function (object, n.sims=100){
# sim.nls:  get posterior simulations of sigma and beta from a nls object
#
# Arguments:
#
#     object:  the output of a call to "nls"
#              with n data points and k predictors
#     n.sims:  number of independent simulation draws to create
#
# Output is a list (sigma.sim, beta.sim):
#
#     sigma.sim:  vector of n.sims random draws of sigma
#       (for glm's, this just returns a vector of 1's or else of the
#       square root of the overdispersion parameter if that is in the model)
#     beta.sim:  matrix (dimensions n.sims x k) of n.sims random draws of beta
#

  object.class <- class(object)[[1]]
  if (object.class!="nls") stop("not a nls object")

    summ <- summary (object)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    sigma.hat <- summ$sigma
    beta.hat <- coef[,1]
    V.beta <- summ$cov.unscaled
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    sigma <- rep (NA, n.sims)
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, names(beta.hat))
    for (s in 1:n.sims){
      sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
      beta[s,] <- mvrnorm (1, beta.hat, V.beta*sigma[s]^2)
    }
    return (list (beta=beta, sigma=sigma))
  }

## self-start function of a 4-parameter logistic function:
## y = al4+(al1-al4)/(1+(x/al3)^al2)
## the function SSfpl is for a different form of fpl:
## A+(B-A)/(1+exp((xmid-input)/scal))
## these two functions are the same (input = log(x))

### Self-starter function for the four parameter logistic function
## the mean function
fplModel <- function(input, al1, al2, al3, al4){
    .x <- input+0.0001
    .expr1 <- (.x/al3)^al2
    .expr2 <- al1-al4
    .expr3 <- 1 + .expr1
    .expr4 <- .x/al3
    .value <- al4 + .expr2/.expr3
    .grad <- array(0, c(length(.value), 4L),
                   list(NULL, c("al1","al2","al3","al4")))
    .grad[,"al1"] <- 1/.expr3
    .grad[,"al2"] <- -.expr2*.expr1*log(.expr4)/.expr3^2
    .grad[,"al3"] <- .expr1*.expr2*(al2/al3)/.expr3^2
    .grad[,"al4"] <- .expr1/(1+.expr1)
    attr(.value, "gradient") <- .grad
    .value
}

## initial values
fplModelInit <- function(mCall, LHS, data, ...){
    xy <- sortedXyData(mCall[["input"]], LHS, data)
   if (nrow(xy) < 5) {
        stop("too few distinct input values to fit a four-parameter logistic")
    }
    rng <- range(xy$y)
    drng <- diff(rng)
    xy$prop <- (xy$y-rng[1]+0.05*drng)/(1.1*drng)
    xy$logx <- log(xy$x+0.0001)
    ir <- as.vector(coef(lm(I(log(prop/(1-prop))) ~ logx, data=xy)))
    pars <- as.vector(coef(nls(y~cbind(1, 1/(1+(x/exp(lal3))^al2)),
                               data=xy,
                               start=list(al2=-ir[2], lal3=-ir[1]/ir[2]),
                               algorithm="plinear")))
    value <- c(pars[4]+pars[3], pars[1], exp(pars[2]), pars[3])
    names(value) <- mCall[c("al1","al2","al3","al4")]
    value
}

SSfpl2 <- selfStart(fplModel, fplModelInit, c("al1","al2","al3","al4"))

base <- getwd()
dataDIR <- paste(base, "Data", sep="/")
plotDIR <- paste(base, "figures", sep="/")

packages(rv)
packages(rstan)
packages(tidyverse)
packages(tikzDevice)
packages(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Importing data
Toledo <- read.csv(paste(dataDIR, "ToledoCrisis.csv", sep="/"),
                   header=T)
names(Toledo)[1] <- "SampleID"

head(Toledo)

## plotting the raw data (not shown in the paper)
tikz(file=paste(plotDIR, "TolData.tex", sep="/"),
     height=3, width=4, standAlone=T)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2, 0), las=1, tck=-0.01)
plot(Absorbance ~ jitter(Concentration), data=Toledo, col=Toledo$Test,
     xlab="MC concentration ($\\mu$g/L)", ylab="Abs")
##ggplot(data=Toledo, aes(y=Absorbance, x=jitter(Concentration), color=Test))+
##    geom_point()
dev.off()

##########Part 1: Stan models and rstan implementation ###########
### Default ELISA model -- simulating nonlinear regression
###  with standard solutions only
### Stan Model (single test kit)
### The model: y=th4+(th1-th4)/[1+(x/th3)^(-th2)]
ELISA1 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  real y[N]; // response data
  real x[N]; // predictor data
}
parameters {
  real<lower=0> delta;
  real<lower=0> th1;
  real<lower=0> th2;
  real<lower=0> th3;
  real<lower=0> sigma;
}
transformed parameters{
  real mu[N];
  real<lower=0> th4;
  th4 = th1+delta;
  for (i in 1:N){
    mu[i] = th4 - delta/(1+(x[i]/th3)^(-th2));
  }
}
model {
  delta~normal(1,5);
  th1~normal(0,5);
  th2~normal(0,5);
  th3~normal(0,5);
  target += normal_lpdf(y| mu, sigma);
}
"
stan.fit1 <- stan_model(model_code=ELISA1)

stan.in1 <- function(data = Toledo, chains=nchains){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    y <- data$Absorbance[!temp]
    x <- data$Concentration[!temp]
    N <- n-sum(temp)
    stan.dat <- list(N=N, y=y, x=x)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(delta = runif(1,0.5,1), th2=runif(1),
                           th3 = runif(1), th1=runif(1,0,0.2),
                           sigma = runif(1))
    parameters <- c("th1", "th2","th3","th4","sigma")
    return(list(para=parameters, data=stan.dat, inits=inits,
                n.chains=chains))
}

test_stan1 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in1(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit1, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    test_stan1[[i]] <- rstan::extract(fit2keep)
}
save(test_stan1, file="ELISA_stan1.RData")
##load("ELISA_stan1.RData")

plot(test_stan1[[1]]$th1, log(test_stan1[[1]]$sigma))
## Neal's funnel -- explaining the reason of high levels of
##  uncertainty when using small sample sizes
##  Add a case where only five data points were used
## use the ELISA test data with all known concentrations to
##  show if the Neal's funnel disappears

## five data point case

Toledo5 <- Toledo[!is.na(Toledo$Concentration), ]
Toledo5mean <- data.frame(
    Absorbance=tapply(Toledo5$Absorbance,
                      paste(Toledo5$SampleID, Toledo5$Test),
                      mean),
    Concentrationc=tapply(Toledo5$Concentration,
                      paste(Toledo5$SampleID, Toledo5$Test),
                      mean),
    Test=tapply(Toledo5$Test,
                      paste(Toledo5$SampleID, Toledo5$Test),
                      mean))

stan.in2 <- function(data = Toledo5mean, chains=nchains){
    Abs0 <- data$Absorbance[data$Concentration==0]
    data <- data[data$Concentration!=0,]
    n <- dim(data)[1]
    y <- data$Absorbance/Abs0
    x <- data$Concentration
    stan.dat <- list(N=n, y=y, x=x)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(delta = runif(1,0.5,1), th2=runif(1),
                           th3 = runif(1), th1=runif(1,0,0.2),
                           sigma = runif(1))
    parameters <- c("th1", "th2","th3","th4","sigma")
    return(list(para=parameters, data=stan.dat, inits=inits,
                n.chains=chains))
}

test_stan2 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <-
        stan.in2(data=Toledo5mean[Toledo5mean$Test==i,])
    fit2keep <- sampling(stan.fit1, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    test_stan2[[i]] <- rstan::extract(fit2keep)
}
save(test_stan2, file="ELISA_stan2.RData")
## load("ELISA_stan2.RData")


### Default ELISA model with water samples
###    unknown concentrations estimated using MLE
### Stan Model (single test kit)
### The model: y=th4+(th1-th4)/[1+(x/th3)^(-th2)]
ELISA2 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // response data
  real x[N]; // predictor data
  real<lower=1> dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples
}
parameters {
  real<lower=0> delta;
  real<lower=0> th1;
  real<lower=0> th2;
  real<lower=0> th3;
  real<lower=0> sigma;
  vector[MM] logx0;
}
transformed parameters{
  real mu[N];
  real mu0[M];
  vector<lower=0>[MM] x0;
  real<lower=0> th4;
  th4 = th1+delta;
  x0 = exp(logx0);
  for (i in 1:N){
    mu[i] = th4 - delta/(1+(x[i]/th3)^(-th2));
  }
  for (i in 1:M){
    mu0[i] = th4 - delta/(1+((x0[ws[i]]/dlf[i])/th3)^(-th2));
  }
}
model {
  logx0 ~ normal(0,2.5);
  delta~normal(1,1);
  th1~normal(0,1);
  th2~normal(0,5);
  th3~normal(0,5);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"
stan.fit2 <- stan_model(model_code=ELISA2)

stan.in3 <- function(data = Toledo, chains=nchains){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    y <- data$Absorbance[!temp]
    x <- data$Concentration[!temp]
    N <- n-sum(temp)
    M <- sum(temp)
    y0 <- data$Absorbance[temp]
    dlf <- data$Dilution[temp]
    wsC<-ordered(data$SampleID[temp])
    ws <- as.numeric(wsC)
    stan.dat <- list(N=N, M=M, MM=max(ws), y=y, x=x,
                     y0=y0, ws=ws, dlf=dlf)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(delta = runif(1,0.5,1), th2=runif(1),
                           th3 = runif(1), th1=runif(1,0,0.1),
                           logx0 = rnorm(max(ws)),
                           sigma = runif(1))
    parameters <- c("th1","th2","th3","th4","sigma","x0")
    return(list(para=parameters, data=stan.dat, inits=inits,
                n.chains=chains, ControlSample=(1:max(ws))[levels(wsC)=="NConl"],
                TheSample=(1:max(ws))[levels(wsC)=="TUL8_1"]))
}

### extracting information from the input data sets for plotting purposes
test_stanMLE1 <- list()
controlS <- numeric()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in3(data=Toledo[Toledo$Test==i,])
    controlS[i] <- input.to.stan$ControlSample
}
thesample <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in3(data=Toledo[Toledo$Test==i,])
    thesample[[i]] <- input.to.stan$TheSample
}

## Running the stan model separately for 6 tests
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in3(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    test_stanMLE1[[i]] <- rstan::extract(fit2keep)
    controlS[i] <- input.to.stan$ControlSample
}
save(test_stanMLE1, file="ELISA_stan3.RData")
## load("ELISA_stan3.RData")

### rmss -- compare the fitted and predicted rmss (wrt y, not used in paper)
rmssF <- rmssP <- rv(0)
for (i in 1:6){
    input.to.stan <- stan.in3(data=Toledo[Toledo$Test==i,])
    tmp <- rvsims(as.matrix(as.data.frame(test_stanMLE1[[i]])))
    yhat <- tmp[4] -
        (tmp[4]-tmp[1])/(1+(input.to.stan$data$x/tmp[3])^(-tmp[2]))
    rmssF <- c(rmssF, mean((input.to.stan$data$y - yhat)^2))
    x0 <- tmp[substring(names(tmp), 1,2)=="x0"]
    x0hat <- x0[input.to.stan$data$ws]/input.to.stan$data$dlf
##print(c(summary(x0[input.to.stan$TheSample])))
    y0hat <- tmp[4] - (tmp[4]-tmp[1])/(1+(x0hat/tmp[3])^(-tmp[2]))
    rmssP <- c(rmssP, mean((input.to.stan$data$y0-y0hat)^2))
}

rmss <- sqrt(rvmatrix(c(rmssF, rmssP), ncol=2,
                      dimnames=list(paste("Test", 1:6),
                                     c("Fitted","Predicted"))))
mlplot(rmss, xlab="Root mean sum of squares")
## load("ELISA_stan3.RData")

### Default ELISA model with water samples
###    unknown concentrations estimaterd using hierarchical model
### Stan Model (single test kit)
### The model: y=th4+(th1-th4)/[1+(x/th3)^(-th2)]

ELISA3 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // response data
  real x[N]; // predictor data
  real<lower=1> dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples
}
parameters {
  real<lower=0> delta;
  real<lower=0> th1;
  real<lower=0> th2;
  real<lower=0> th3;
  real<lower=0> sigma;
  vector[MM] logx0;
  real<lower=0> x0sig;
  real x0mu;
}
transformed parameters{
  vector[N] mu;
  vector[M] mu0;
  vector<lower=0>[MM] x0;
  real<lower=0> th4;
  x0 = exp(logx0);
  th4 = th1+delta;
  for (i in 1:N){
    mu[i] = th4 - delta/(1+(x[i]/th3)^(-th2));
  }
  for (i in 1:M){
    mu0[i] = th4 - delta/(1+((x0[ws[i]]/dlf[i])/th3)^(-th2));
  }
}
model {
  x0mu ~ normal(0, 2.5);
  x0sig ~ normal(0, 2.5);
  delta~normal(1,1);
  logx0 ~ normal(x0mu, x0sig);
  th1~normal(0,1);
  th2~normal(0,5);
  th3~normal(0,5);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"
stan.fit3 <- stan_model(model_code=ELISA3)

stan.in4 <- function(data = Toledo, chains=nchains){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    y <- data$Absorbance[!temp]
    x <- data$Concentration[!temp]
    N <- n-sum(temp)
    M <- sum(temp)
    y0 <- data$Absorbance[temp]
    dlf <- data$Dilution[temp]
    wsC <- ordered(data$SampleID[temp])
    ws <- as.numeric(wsC)
    stan.dat <- list(N=N, M=M, MM=max(ws), y=y, x=x,
                     y0=y0, ws=ws, dlf=dlf)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(delta = runif(1, 0.5, 1), th2=runif(1, 0.4, 0.6),
                           th3 = runif(1, 0.3, 0.7), th1=runif(1, 0, 0.1),
                           logx0 = rnorm(max(ws)),
                           sigma = runif(1), x0mu=rnorm(1), x0sig=runif(1))
    parameters <- c("th1", "th2","th3","th4","sigma","x0", "x0mu","x0sig")
    return(list(para=parameters, data=stan.dat, inits=inits, n.chains=chains,
                ControlSample=(1:max(ws))[levels(wsC)=="NConl"],
                TheSample=(1:max(ws))[levels(wsC)=="TUL8_1"]))
}

test_stanBHM1 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in4(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    test_stanBHM1[[i]] <- rstan::extract(fit2keep)
}
save(test_stanBHM1, file="ELISA_stan4.RData")
## load("ELISA_stan4.RData")
##pairs(fit2keep, pars=c("th1","th2","th3","th4","sigma"))

### Default ELISA model with water samples
###    unknown concentrations estimaterd using hierarchical model
### Stan Model (single test kit)
### The model: y= A-D/(1+exp((xmid-x)/scal))
### Comparing fitting the four-parameter logistic function in the log-concentration scale. (Not used)

ELISA3.5 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  int<lower=1> N0;
  real zeroY[N0]; // observed 0 conc OD
  real y[N]; // response data
  real x[N]; // predictor data (log)
  real<lower=1> dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples
}
parameters {
  real<lower=0> delta;
  real<lower=0> A;
  real<lower=0> scal;
  real xmid;
  real<lower=0> sigma;
  vector[MM] x0;
  real<lower=0> x0sig;
  real x0mu;
}
transformed parameters{
  vector[N] mu;
  vector[M] mu0;
  real<lower=0> B;

  B = A+delta;
  for (i in 1:N){
    mu[i] = A - delta/(1+exp((xmid - x[i])/scal));
  }
  for (i in 1:M){
    mu0[i] = A - delta/(1+exp((xmid-(x0[ws[i]]-log(dlf[i])))/scal));
  }
}
model {
  x0mu ~ normal(0, 2.5);
  x0sig ~ normal(0, 2.5);
  delta~normal(1,1);
  x0 ~ normal(x0mu, x0sig);
  A~normal(0,1);
  xmid~normal(0,2.5);
  scal~normal(0,2.5);
  target += normal_lpdf(zeroY | A, sigma);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"
stan.fit3.5 <- stan_model(model_code=ELISA3.5)

stan.in4.5 <- function(data = Toledo, chains=nchains){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    y <- data$Absorbance[!temp]
    x <- data$Concentration[!temp]
    N0 <- sum(x==0)
    zeroy <- y[x==0]
    y <- y[x!=0]
    x <- log(x[x!=0])
    N <- length(x)
    M <- sum(temp)
    y0 <- data$Absorbance[temp]
    dlf <- data$Dilution[temp]
    wsC <- ordered(data$SampleID[temp])
    ws <- as.numeric(wsC)
    stan.dat <- list(N=N, M=M, MM=max(ws), N0=N0,
                     y=y, x=x, zeroY=zeroy,
                     y0=y0, ws=ws, dlf=dlf)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(delta = runif(1, 0.5, 1), xmid=rnorm(1),
                           scal = runif(1, 0.3, 0.7), A=runif(1, 0, 0.1),
                           x0 = runif(max(ws)),
                           sigma = runif(1), x0mu=runif(1), x0sig=runif(1))
    parameters <- c("A", "B","xmid","scal","sigma","x0", "x0mu","x0sig")
    return(list(para=parameters, data=stan.dat, inits=inits, n.chains=chains,
                ControlSample=(1:max(ws))[levels(wsC)=="NConl"]))
}

test_stanBHM3.5 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in4.5(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit3.5, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    test_stanBHM3.5[[i]] <- rstan::extract(fit2keep)
}

save(test_stanBHM3.5, file="ELISA_stan4.5.RData")
## load("ELISA_stan4.5.RData")
##pairs(fit2keep, pars=c("th1","th2","th3","th4","sigma"))


### Hierarchical ELISA model with water samples
###    unknown concentrations estimaterd using hierarchical model
### Stan Model (multiple test kits)
### The model: y=th4+(th1-th4)/[1+(x/th3)^(-th2)]

ELISA4 <- "
data {
  int<lower=0> N; // sample size
  int<lower=0> K; // number of test kits
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // std solution response data
  real x[N]; // std solution conc
  real<lower=1> dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples
  int<lower=1,upper=K> test1[N];
  int<lower=1,upper=K> test2[M];
  int<lower=1,upper=K> test3[MM];
}
parameters {
  vector[K] zdelta;
  vector[K] zth1;
  vector[K] zth2;
  vector[K] zth3;
  vector<lower=0>[K] sigma;
  vector[MM] logx0;
  vector<lower=0>[K] x0sig;
  vector[K] x0mu;
  real mu_dlt;
  real<lower=0> sig_dlt;
  real mu_th1;
  real<lower=0> sig_th1;
  real mu_th2;
  real<lower=0> sig_th2;
  real mu_th3;
  real<lower=0> sig_th3;
}
transformed parameters{
  vector<lower=0>[K] delta;
  vector<lower=0>[K] th1;
  vector<lower=0>[K] th2;
  vector<lower=0>[K] th3;
  vector[N] mu;
  vector[M] mu0;
  vector<lower=0>[N] sig;
  vector<lower=0>[M] sig0;
  vector<lower=0>[K] th4;
  vector<lower=0>[MM] x0;

  x0 = exp(logx0);
  delta = mu_dlt+sig_dlt*zdelta;
  th1 = mu_th1+sig_th1*zth1;
  th2 = mu_th2+sig_th2*zth2;
  th3 = mu_th3+sig_th3*zth3;
  th4 = th1+delta;
  for (i in 1:N){
    mu[i] = th4[test1[i]] -
              delta[test1[i]]/(1+(x[i]/th3[test1[i]])^(-th2[test1[i]]));
    sig[i] = sigma[test1[i]];
  }
  for (i in 1:M){
    mu0[i] = th4[test2[i]] -
       delta[test2[i]]/(1+((x0[ws[i]]/dlf[i])/th3[test2[i]])^(-th2[test2[i]]));
    sig0[i] = sigma[test2[i]];
  }
}
model {
  x0mu ~ normal(0, 2.5);
  x0sig ~ normal(0, 2.5);
  mu_dlt ~ normal(0, 2.5);
  sig_dlt ~ normal(0, 2.5);
  mu_th1 ~ normal(0, 2.5);
  sig_th1 ~ normal(0, 2.5);
  mu_th2 ~ normal(0, 2.5);
  sig_th2 ~ normal(0, 2.5);
  mu_th3 ~ normal(0, 2.5);
  sig_th3 ~ normal(0, 2.5);

  for (j in 1:MM)
    logx0[j] ~ normal(x0mu[test3[j]], x0sig[test3[j]]);
  zth1~std_normal();
  zth2~std_normal();
  zth3~std_normal();
  zdelta~std_normal();
  target += normal_lpdf(y| mu, sig);
  target += normal_lpdf(y0| mu0, sig0);
}
"

stan.fit4 <- stan_model(model_code=ELISA4)

stan.in5 <- function(data = Toledo, chains=nchains){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    tmp1 <- data[!temp,]
    tmp2 <- data[temp,]
    test1 <- as.numeric(ordered(tmp1$Test))
    y <- tmp1$Absorbance
    x <- tmp1$Concentration
    N <- dim(tmp1)[1]
    M <- dim(tmp2)[1]
    K <- max(test1)
    y0 <- tmp2$Absorbance
    dlf <- tmp2$Dilution
    wsC <- ordered(paste(tmp2$Test,tmp2$SampleID))
    wsLevels <- levels(wsC)
    tmpp <- substring(wsLevels, 3, 8)=="NConl"
    tmpp2 <- substring(wsLevels, 3, 8) == "TUL8_1"
    ws <- as.numeric(wsC)
    test2 <- as.numeric(ordered(tmp2$Test))
    test3 <- as.numeric(substring(wsLevels, 1, 1))
    stan.dat <- list(N=N, M=M, MM=max(ws), K=K, y=y, x=x,
                     y0=y0, ws=ws, dlf=dlf, test1=test1,
                     test2=test2, test3=test3)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(zdelta=runif(K,0.5,1), zth1=runif(K,0,0.1),
                           zth2=runif(K,0.4,0.6), zth3=runif(K,0.3,0.7),
                           sigma=runif(K), logx0=rnorm(max(ws)),
                           x0sig=runif(K), x0mu=rnorm(K),
                           mu_dlt=runif(1), sig_dlt=runif(1),
                           mu_th1=runif(1), sig_th1=runif(1),
                           mu_th2=runif(1), sig_th2=runif(1),
                           mu_th3=runif(1), sig_th3=runif(1))

    parameters <- c("th1", "th2","th3","th4","sigma","x0",
                    "x0mu","x0sig", "mu_dlt", "sig_dlt", "mu_th1",
                    "sig_th1", "mu_th2", "sig_th2", "mu_th3", "sig_th3")
    return(list(para=parameters, data=stan.dat, inits=inits, n.chains=chains,
                ControlSample=(1:max(ws))[tmpp],
                TheSample=(1:max(ws))[tmpp2]))
}

input.to.stan <- stan.in5(data=Toledo)
BHM2ControlSamples <- input.to.stan$ControlSample
TheSample <- input.to.stan$TheSample
fit2keep <- sampling(stan.fit4, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)
test_stanBHM2 <- rstan::extract(fit2keep)
save(test_stanBHM2, file="ELISA_stan5.RData")
## load("ELISA_stan5.RData")

### Sequential updating ELISA model with water samples
###    unknown concentrations estimaterd using hierarchical model
### Stan Model (single multiple test kit)
### The model: y=th4+(th1-th4)/[1+(x/th3)^(-th2)]
### Not used in the paper

ELISA5 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; //observed ODs (for estimating concentration)
  real y[N]; // response data
  real x[N]; // predictor data
  real dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples

  real deltamu0; //prior means
  real th1mu0;
  real th2mu0;
  real th3mu0;
  real<lower=0> deltan0; // prior lambda
  real<lower=0> th1n0;
  real<lower=0> th2n0;
  real<lower=0> th3n0;
  real<lower=0> deltalpha; // prior alpha
  real<lower=0> th1alpha;
  real<lower=0> th2alpha;
  real<lower=0> th3alpha;
  real<lower=0> deltabeta; // prior beta
  real<lower=0> th1beta;
  real<lower=0> th2beta;
  real<lower=0> th3beta;
  real<lower=0> aa; // prior a (sigma_y)
  real<lower=0> bb; // prior b
}
parameters {
  real<lower=0> delta;
  real<lower=0> th1;
  real<lower=0> th2;
  real<lower=0> th3;
  real<lower=0> sig2;
  vector[MM] logx0;
  real<lower=0> x0sig;
  real x0mu;
  real mu_dlt;
  real<lower=0> sig2_dlt;
  real mu_th1;
  real<lower=0> sig2_th1;
  real mu_th2;
  real<lower=0> sig2_th2;
  real mu_th3;
  real<lower=0> sig2_th3;
}
transformed parameters{
  vector[N] mu;
  vector[M] mu0;
  real<lower=0> th4;
  vector<lower=0>[MM] x0;

  x0 = exp(logx0);
  th4 = th1+delta;
  for (i in 1:N){
    mu[i] = th4 - delta/(1+(x[i]/th3)^(-th2));
  }
  for (i in 1:M){
    mu0[i] = th4 - delta/(1+((x0[ws[i]]/dlf[i])/th3)^(-th2));
  }
}
model {
  x0mu ~ normal(0, 2.5);
  x0sig ~ normal(0, 2.5);
  logx0 ~ normal(x0mu, x0sig);

  sig2 ~ inv_gamma(aa, bb);
  sig2_dlt ~ inv_gamma(deltalpha, deltabeta);
  mu_dlt ~ normal(deltamu0, sqrt(sig2_dlt/deltan0));
  sig2_th1 ~ inv_gamma(th1alpha, th1beta);
  mu_th1 ~ normal(th1mu0, sqrt(sig2_th1/th1n0));
  sig2_th2 ~ inv_gamma(th2alpha, th2beta);
  mu_th2 ~ normal(th2mu0, sqrt(sig2_th2/th2n0));
  sig2_th3 ~ inv_gamma(th3alpha, th3beta);
  mu_th3 ~ normal(th3mu0, sqrt(sig2_th3/th3n0));

  delta~normal(mu_dlt,sqrt(sig2_dlt));
  th1~normal(mu_th1, sqrt(sig2_th1));
  th2~normal(mu_th2, sqrt(sig2_th2));
  th3~normal(mu_th3, sqrt(sig2_th3));

  target += normal_lpdf(y| mu, sqrt(sig2));
  target += normal_lpdf(y0| mu0, sqrt(sig2));
}
"

stan.fit5 <- stan_model(model_code=ELISA5)

stan.in6 <- function(data = Toledo[Toledo$Test==1,], chains=nchains,
                     Dstan, th1stan, th2stan, th3stan, sig2stan){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    tmp1 <- data[!temp,]
    tmp2 <- data[temp,]
    test1 <- as.numeric(ordered(tmp1$Test))
    y <- tmp1$Absorbance
    x <- tmp1$Concentration
    N <- dim(tmp1)[1]
    M <- dim(tmp2)[1]
    y0 <- tmp2$Absorbance
    dlf <- tmp2$Dilution
    wsC <- ordered(paste(tmp2$Test,tmp2$SampleID))
    wsLevels <- levels(wsC)
    tmpp <- substring(wsLevels, 3, 8)=="NConl"
    tmpp2 <- substring(wsLevels, 3, 8) == "TUL8_1"
    ws <- as.numeric(wsC)
    stan.dat <- list(N=N, M=M, MM=max(ws), y=y, x=x,
                     y0=y0, ws=ws, dlf=dlf,
                     deltamu0=Dstan$mu0, th1mu0=th1stan$mu0,
                     th2mu0=th2stan$mu0, th3mu0=th3stan$mu0,
                     deltan0=Dstan$lambda, th1n0=th1stan$lambda,
                     th2n0=th2stan$lambda, th3n0=th3stan$lambda,
                     deltalpha=Dstan$alpha, th1alpha=th1stan$alpha,
                     th2alpha=th2stan$alpha, th3alpha=th3stan$alpha,
                     deltabeta=Dstan$beta, th1beta=th1stan$beta,
                     th2beta=th2stan$beta, th3beta=th3stan$beta,
                     aa=sig2stan$alpha, bb=sig2stan$beta)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(zdelta=runif(1,0.5,1), zth1=runif(1,0,0.1),
                           zth2=runif(1,0.4,0.6), zth3=runif(1,0.3,0.7),
                           sig2=runif(1), logx0=rnorm(max(ws)),
                           x0sig=runif(1), x0mu=runif(1),
                           mu_dlt=runif(1), sig2_dlt=runif(1),
                           mu_th1=runif(1), sig2_th1=runif(1),
                           mu_th2=runif(1), sig2_th2=runif(1),
                           mu_th3=runif(1), sig2_th3=runif(1))

    parameters <- c("th1", "th2","th3","th4","sig2","x0",
                    "x0mu","x0sig", "mu_dlt", "sig2_dlt", "mu_th1",
                    "sig2_th1", "mu_th2", "sig2_th2", "mu_th3", "sig2_th3")
    return(list(para=parameters, data=stan.dat, inits=inits, n.chains=chains,
                ControlSample=(1:max(ws))[tmpp],
                TheSample=(1:max(ws))[tmpp2]))
}

prior_pars_NIG <- function(mus, sigs2){
  Ex <- mean(mus)
  Vx <- sd(mus)^2
  Esig2 <- mean(sigs2)
  Vsig2 <- sd(sigs2)^2
  return(list(mu0=Ex, beta=Esig2*(1+Esig2^2/Vsig2),
              alpha=2+Esig2^2/Vsig2, lambda=Esig2/Vx))
}

prior_pars_IG <- function(sig2){
    Esig2 <- mean(sig2)
    Vsig2 <- sd(sig2)^2
    return(list(alpha=2+Esig2^2/Vsig2, beta=Esig2*(1+Esig2^2/Vsig2)))
}

sig2_stan <- prior_pars_IG(test_stanBHM2$sigma^2)
D_stan <- prior_pars_NIG(test_stanBHM2$mu_dlt, test_stanBHM2$sig_dlt^2)
th1_stan <- prior_pars_NIG(test_stanBHM2$mu_th1, test_stanBHM2$sig_th1^2)
th2_stan <- prior_pars_NIG(test_stanBHM2$mu_th2, test_stanBHM2$sig_th2^2)
th3_stan <- prior_pars_NIG(test_stanBHM2$mu_th3, test_stanBHM2$sig_th3^2)

seqBHM <- list()
for (i in 1:6){
    print(paste("Test: ", i))
          input.to.stan <- stan.in6(data=Toledo[Toledo$Test==i,],
                              Dstan=D_stan, th1stan=th1_stan,
                              th2stan=th2_stan, th3stan=th3_stan,
                              sig2stan=sig2_stan)

    fit2keep <- sampling(stan.fit5, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    seqBHM[[i]]<- rstan::extract(fit2keep)
    sig2_stan <- prior_pars_IG(seqBHM[[i]]$sig2)
    D_stan <- prior_pars_NIG(seqBHM[[i]]$mu_dlt, seqBHM[[i]]$sig2_dlt)
    th1_stan <- prior_pars_NIG(seqBHM[[i]]$mu_th1, seqBHM[[i]]$sig2_th1)
    th2_stan <- prior_pars_NIG(seqBHM[[i]]$mu_th2, seqBHM[[i]]$sig2_th2)
    th3_stan <- prior_pars_NIG(seqBHM[[i]]$mu_th3, seqBHM[[i]]$sig2_th3)
}

save(seqBHM, file="ELISA_stan6.RData")
## load("ELISA_stan6.RData")


######################## Part 2: Figures ########################

### 1. Inverse-function method
#### Sample size comparisons (Figure 2)
smpl <- sample(1:2504, size=200)
tikz(file=paste(plotDIR, "samplesize_eff.tex", sep="/"),
     height=2.5, width=5, standAlone=F)
par(mfrow=c(2,4), mar=c(3,3,0,0.5), mgp=c(1.25,0.125,0), las=1, tck=-0.01)
plot(test_stan2[[6]]$th1[smpl], log(test_stan2[[6]]$sigma[smpl]),
     xlim=c(0, 0.6), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_1$", ylab="$\\log(\\sigma)$")
text(0.4, -3.5, "$n=5$")
plot(test_stan2[[6]]$th2[smpl], log(test_stan2[[6]]$sigma[smpl]),
     xlim=c(0, 2), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_2$", ylab="$\\log(\\sigma)$")
plot(test_stan2[[6]]$th3[smpl], log(test_stan2[[6]]$sigma[smpl]),
     xlim=c(0, 2), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_3$", ylab="$\\log(\\sigma)$")
plot(test_stan2[[6]]$th4[smpl], log(test_stan2[[6]]$sigma[smpl]),
     xlim=c(0, 2), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_4$", ylab="$\\log(\\sigma)$")
plot(test_stan1[[6]]$th1[smpl], log(test_stan1[[6]]$sigma[smpl]),
     xlim=c(0, 0.6), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_1$", ylab="$\\log(\\sigma)$")
text(0.4, -3.5, "$n=12$")
plot(test_stan1[[6]]$th2[smpl], log(test_stan1[[6]]$sigma[smpl]),
     xlim=c(0, 2), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_2$", ylab="$\\log(\\sigma)$")
plot(test_stan1[[6]]$th3[smpl], log(test_stan1[[6]]$sigma[smpl]),
     xlim=c(0, 2), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_3$", ylab="$\\log(\\sigma)$")
plot(test_stan1[[6]]$th4[smpl], log(test_stan1[[6]]$sigma[smpl]),
     xlim=c(0, 2), ylim=c(-5, 0), cex=0.75,
     xlab="$\\theta_4$", ylab="$\\log(\\sigma)$")
dev.off()

## Predictive uncertainty using inverse function
##  The inverse function log(x) = log(th3) + log((th1-y)/(y-th4))/th2
## through nls
## residual variance corrected based on Mandel (1958), Lieberman (1961),
## Lieberman and Miller (1963), Lieberman, Miller, and Hamilton (1967)
## Simultaneous discreminant intervals

predToledo <- function(data=Toledo[Toledo$Test==2,]){##  12 data points
    data_std <- data[!is.na(data$Concentration),]
    M <- nls(Absorbance ~ SSfpl2(Concentration, al1, al2, al3, al4),
             data=data_std)
    SXX <- sum((data_std$Concentration-mean(data_std$Concentration))^2)
    sigma2_x <- sd(data_std$Concentration)^2
    Msims <- sim.nls(M, 4000)
    betas <- rvsims(Msims$beta)
    sigma <- rvsims(Msims$sigma)
    data_smp <- data[is.na(data$Concentration),]
    m <- dim(data_smp)[1]
    n <- length(data_smp$Absorbance)
##  y <- data_smp$Absorbance+rvnorm(n, mean=0, sd=sigma)#sqrt((8+m)/8)*sigma)
    y <- data_smp$Absorbance+
        rvnorm(n, mean=0,sd=sqrt((1+1/n+sigma2_x*rvchisq(1,m)/SXX))*sigma)
    ## based on linear regression predictive variance formula
    temp1 <- sims(betas[1]-y)
    numm <- rvsims(t(apply(temp1, 1, function(x)ifelse(x<0, 0.0001, x))))
    temp2 <- sims(y-betas[4])
    denm <- rvsims(t(apply(temp2, 1, function(x)ifelse(x<0, 0.0001, x))))
    conc_pred <- exp(log(betas[3])+log(numm/denm)/betas[2])*data_smp$Dilution
    temp <- summary(conc_pred)
    temp$SampleID <- data_smp$SampleID
    temp$Dilution <- data_smp$Dilution
    return(list(summ=temp, rvs=conc_pred))
}

predToledo5 <- function(data=Toledo[Toledo$Test==1,]){##  5 data points
    data$ID <- paste(data$SampleID, data$Dilution)
    data5 <- data.frame(Concentration=tapply(data$Concentration, data$ID, mean),
                        Absorbance=tapply(data$Absorbance, data$ID, mean),
                        Test=tapply(data$Test, data$ID, unique),
                        SampleID=tapply(data$SampleID, data$ID,unique))
    data_std <- data5[!is.na(data5$Concentration),]
    std0mean <- data_std$Absorbance[data_std$Concentration==0]
    data_std$Absorbance <- data_std$Absorbance/std0mean
    M <- nls(Absorbance ~ SSfpl2(Concentration, al1, al2, al3, al4),
             data=data_std)
    SXX <- sum((data_std$Concentration-mean(data_std$Concentration))^2)
    sigma2_x <- sd(data_std$Concentration)^2
    Msims <- sim.nls(M, 4000)
    betas <- rvsims(Msims$beta)
    sigma <- rvsims(Msims$sigma)
    data_smp <- data5[is.na(data5$Concentration),]
    m <- dim(data_smp)[1]/2
    n <- length(data_smp$Absorbance)
    data_smp$Absorbance <- data_smp$Absorbance/std0mean
##    y <- data_smp$Absorbance+rvnorm(mean=0, sd=sigma)# sqrt(1+m)*sigma)
    y <- data_smp$Absorbance+
        rvnorm(n,mean=0,sd=sqrt((1+1/n+sigma2_x*rvchisq(1,m)/SXX))*sigma)
    ## based on linear regression predictive variance formula
    temp1 <- sims(betas[1]-y)
    numm <- rvsims(t(apply(temp1, 1, function(x)ifelse(x<0, 0.0001, x))))
    temp2 <- sims(y-betas[4])
    denm <- rvsims(t(apply(temp2, 1, function(x)ifelse(x<0, 0.0001, x))))
    conc_pred <- exp(log(betas[3])+log(numm/denm)/betas[2])
    temp <- summary(conc_pred)
    temp$SampleID <- data_smp$SampleID
    temp$Dilution <- data_smp$Dilution
    return(list(summ=temp, rvs=conc_pred))
}


pred_summ12 <- list()
pred_rvs12 <- list()
pred_summ5 <- list()
pred_rvs5 <- list()
for (i in 1:6){
    print(i)
    temp <- predToledo(data=Toledo[Toledo$Test==i,])
    pred_summ12[[i]] <- temp$summ
    pred_rvs12[[i]] <- temp$rvs

    temp <- predToledo5(data=Toledo[Toledo$Test==i,])
    pred_summ5[[i]] <- temp$summ
    pred_rvs5[[i]] <- temp$rvs
}

### Figure 3
tikz(file=paste(plotDIR, "mlelogstdv.tex", sep="/"),
     height=2.25, width=2, standAlone=F)
par(mar=c(3,3,1,0), mgp=c(1.25,0.125,0), las=1, tck=-0.01)
boxplot(log(pred_summ12[[3]]$sd), log(pred_summ5[[3]]$sd),
        names=c("$n=12$", "$n=5$"), ylab="log standard deviation",
        cex.lab=0.75)
dev.off()
####

### Comparing estimated values between BHM and inverse function method
### medians and middle 50% CIs (Not used in the paper)
par(mfrow=c(2,3), oma=c(2, 2, 0, 0), mar=c(2,2,1,1), mgp=c(1.25,0.125,0),
    tck=-0.01, las=0)
for (i in 1:6){
    smples <- stan.in4(data=Toledo[Toledo$Test==i,])$data$ws
    tmpX <- pred_summ12[[i]]
    tmpY <- summary(rvsims(test_stanMLE1[[i]]$x0))[smples, ]
    plot(x=tmpX[,6], y=tmpY[,6], pch=16, cex=1.25,
         xlab="",ylab="", log="xy",
         xlim=range(tmpX[, c(5,7)]), ylim=range(tmpY[, c(5,7)]))
    segments(x0=tmpX[,5], x1=tmpX[,7],
             y0=tmpY[,6], y1=tmpY[,6],
             col="grey")
    segments(y0=tmpY[,5], y1=tmpY[,7],
             x0=tmpX[,6], x1=tmpX[,6],
             col="grey")
}
mtext("inverse",side=1, outer=T)
mtext("$MLE_{12}$", side=2, outer=T)

par(mfrow=c(2,3), oma=c(2, 2, 0, 0), mar=c(2,2,1,1), mgp=c(1.25,0.125,0),
    tck=-0.01, las=0)
for (i in 1:6){
    smples <- stan.in4(data=Toledo[Toledo$Test==i,])$data$ws
    tmpX <- pred_summ12[[i]]
    tmpY <- summary(rvsims(test_stanBHM1[[i]]$x0))[smples, ]
    plot(x=tmpX$"50%", y=tmpY$"50%",
         xlab="", ylab="", log="xy",
         xlim=range(tmpX[, c(5,7)]), ylim=range(tmpY[, c(5,7)]))
    segments(x0=tmpX[,5], x1=tmpX[,7],
             y0=tmpY[,6], y1=tmpY[,6],
             col="grey")
    segments(y0=tmpY$"25%", y1=tmpY$"75%",
             x0=tmpX$"50%", x1=tmpX$"50%",
             col="grey")
}
mtext("inverse",side=1, outer=T)
mtext("$BHM_1$", side=2, outer=T)

par(mfrow=c(2,3), oma=c(2, 2, 0, 0), mar=c(2,2,1,1), mgp=c(1.25,0.125,0),
    tck=-0.01)
for (i in 1:6){
    smples <- stan.in4(data=Toledo[Toledo$Test==i,])$data$ws
    tmpX <- pred_summ12[[i]]
    tmpY <- summary(rvsims(seqBHM[[i]]$x0))[smples, ]
    plot(x=tmpX$"50%", y=tmpY$"50%",
         xlab="", ylab="", log="xy",
         xlim=range(tmpX[, c(5,7)]), ylim=range(tmpY[, c(5,7)]))
    segments(x0=tmpX[,5], x1=tmpX[,7],
             y0=tmpY[,6], y1=tmpY[,6],
             col="grey")
    segments(y0=tmpY$"25%", y1=tmpY$"75%",
             x0=tmpX$"50%", x1=tmpX$"50%",
             col="grey")
}
mtext("inverse",side=1, outer=T)
mtext("$sBHM$", side=2, outer=T)

par(mfrow=c(2,3), oma=c(2, 2, 0, 0), mar=c(2,2,1,1), mgp=c(1.25,0.125,0),
    tck=-0.01, las=0)
par(mfrow=c(2,3), mar=c(3,3,1,1), mgp=c(1.25,0.125,0))
for (i in 1:6){
    smples <- stan.in4(data=Toledo[Toledo$Test==i,])$data$ws
    tmpX <- summary(rvsims(test_stanBHM1[[i]]$x0))[smples, ]
    tmpY <- summary(rvsims(seqBHM[[i]]$x0))[smples, ]
    plot(x=tmpX$"50%", y=tmpY$"50%",
         xlab="",ylab="", log="xy",
         xlim=range(tmpX[, c(5,7)]), ylim=range(tmpY[, c(5,7)]))
    segments(x0=tmpX[,5], x1=tmpX[,7],
             y0=tmpY[,6], y1=tmpY[,6],
             col="grey")
    segments(y0=tmpY$"25%", y1=tmpY$"75%",
             x0=tmpX$"50%", x1=tmpX$"50%",
             col="grey")
}
    mtext("$BHM_1$",side=1, outer=T)
mtext("$BHM_2$", side=2, outer=T)


## comparing control sample fits (Figure 4)

## MLE (back-calculate)

cntrlID <- (1:length(pred_summ12[[1]]))[substring(pred_summ12[[1]]$SampleID, 1, 5)=="NConl"]
CntrlMLE12 <- pred_summ12[[1]][cntrlID,]
CntrlMLErvs12 <- pred_rvs12[[1]][cntrlID]

for (i in 2:6){
    cntrlID <- (1:length(pred_summ12[[i]]))[substring(pred_summ12[[i]]$SampleID, 1, 5)=="NConl"]
    CntrlMLE12 <- rbind(CntrlMLE12[,1:9], pred_summ12[[i]][cntrlID,1:9])
    CntrlMLErvs12 <- c(CntrlMLErvs12, pred_rvs12[[i]][cntrlID])
}

cntrlID <- (1:dim(pred_summ5[[1]])[1])[substring(pred_summ5[[1]]$SampleID, 1, 5)=="NConl"]
CntrlMLE5 <- pred_summ5[[1]][cntrlID,]
CntrlMLErvs5 <- pred_rvs5[[1]][cntrlID]
for (i in 2:6){
    cntrlID <- (1:dim(pred_summ5[[i]])[1])[substring(pred_summ5[[i]]$SampleID, 1, 5)=="NConl"]
    CntrlMLE5 <- rbind(CntrlMLE5[,1:9], pred_summ5[[i]][cntrlID,1:9])
    CntrlMLErvs5 <- c(CntrlMLErvs5, pred_rvs5[[i]][cntrlID])
}

##Bayesian models
CntrlStan1 <- rvsims(test_stanMLE1[[1]]$x0)[controlS[1]]
for (i in 2:6)
    CntrlStan1 <- c(CntrlStan1, rvsims(test_stanMLE1[[i]]$x0)[controlS[i]])

CntrlStan3 <- rvsims(test_stanBHM1[[1]]$x0)[controlS[1]]
for (i in 2:6)
    CntrlStan3 <- c(CntrlStan3, rvsims(test_stanBHM1[[i]]$x0)[controlS[i]])

##CntrlStan3.5 <- rvsims(test_stanBHM3.5[[1]]$x0)[controlS[1]]
##or (i in 2:6)
##CntrlStan3.5 <- c(CntrlStan3.5, rvsims(test_stanBHM3.5[[i]]$x0)[controlS[i]])

CntrlStanSeq <- rvsims(seqBHM[[1]]$x0)[controlS[1]]
for (i in 2:6)
    CntrlStanSeq <- c(CntrlStanSeq, rvsims(seqBHM[[i]]$x0)[controlS[i]])

CntrlStan4 <- rvsims(test_stanBHM2$x0[, BHM2ControlSamples])

Bayesian1.2 <- rvmatrix(c(CntrlStan1, CntrlStan3, CntrlStanSeq, CntrlStan4),
                        nrow=6)
mlplot(t(Bayesian1.2))
abline(v=0.75)

CntrlMLE5
CntrlMLE12
CntrlBay1 <- summary(CntrlStan1)
CntrlBayH <- summary(CntrlStan3)
CntrlBayH1 <- summary(CntrlStanSeq)
CntrlBayH2 <- summary(CntrlStan4)

mlplot(CntrlStan1)
mlplot(CntrlStan3)
mlplot(CntrlStan4)

tikz(file=paste(plotDIR, "ToledoComp.tex", sep="/"),
     width=3, height=3.75, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25, 0.25, 0), las=1, tck=-0.01)
plot(c(0,1),c(0,1), ylim=c(0.75,6.25), xlim=c(0.,4), type="n",
     ylab="Test number", xlab="Estimated MC concentrations ($\\mu$g/L)", cex.lab=0.85)

segments(y0=(1:6)-0.25/1.25, y1=(1:6)-0.25/1.25,
         x0=CntrlMLE12[seq(1,12,2), 4],
         x1=CntrlMLE12[seq(1,12,2), 8])
segments(y0=(1:6)-0.25/1.25, y1=(1:6)-0.25/1.25,
         x0=CntrlMLE12[seq(1,12,2), 5],
         x1=CntrlMLE12[seq(1,12,2), 7], lwd=3)
points(y=(1:6)-0.25/1.25, x=CntrlMLE12[seq(1,12,2), 6], cex=0.75)

segments(y0=(1:6)-0.125/1.25, y1=(1:6)-0.125/1.25,
         x0=CntrlMLE12[seq(2,12,2), 4],
         x1=CntrlMLE12[seq(2,12,2), 8])
segments(y0=(1:6)-0.125/1.25, y1=(1:6)-0.125/1.25,
         x0=CntrlMLE12[seq(2,12,2), 5],
         x1=CntrlMLE12[seq(2,12,2), 7], lwd=3)
points(y=(1:6)-0.125/1.25, x=CntrlMLE12[seq(2,12,2), 6], cex=0.75)

##segments(y0=(1:6), y1=(1:6),
##         x0=CntrlBay1[, 4],
##         x1=CntrlBay1[, 8], col="purple")
##segments(y0=(1:6), y1=(1:6),
##         x0=CntrlBay1[, 5],
##         x1=CntrlBay1[, 7], lwd=3, col="purple")
##points(y=(1:6), x=CntrlBay1[, 6], col="purple", cex=0.75)

segments(y0=(1:6)+0.125/1.25, y1=(1:6)+0.125/1.25,
         x0=CntrlBayH[, 4],
         x1=CntrlBayH[, 8], col="blue")
segments(y0=(1:6)+0.125/1.25, y1=(1:6)+0.125/1.25,
         x0=CntrlBayH[, 5],
         x1=CntrlBayH[, 7], lwd=3, col="blue")
points(y=(1:6)+0.125/1.25, x=CntrlBayH[, 6], col="blue", cex=0.75)

segments(y0=(1:6)+0.25/1.25, y1=(1:6)+0.25/1.25,
         x0=CntrlBayH2[, 4],
         x1=CntrlBayH2[, 8], col="red")
segments(y0=(1:6)+0.25/1.25, y1=(1:6)+0.25/1.25,
         x0=CntrlBayH2[, 5],
         x1=CntrlBayH2[, 7], lwd=3, col="red")
points(y=(1:6)+0.25/1.25, x=CntrlBayH2[, 6], col="red", cex=0.75)

abline(v=0.75, col="gray", lwd=2)
abline(h=seq(1.5,5.5, 1), col="gray", lwd=2)
dev.off()

## Evaluating biasness (Figure 5)

bias <- c(mean(abs(CntrlStan4-0.75)),
          mean(abs(CntrlStanSeq-0.75)),
          mean(abs(CntrlStan3-0.75)),
          mean(abs(CntrlStan1-0.75)),
          mean(abs(CntrlMLErvs12-0.75),na.rm=T),
          mean(abs(CntrlMLErvs5-0.75),na.rm=T))
names(bias) <- c("$BHM_2$", "$BHM_{seq}$", "$BHM_1$",
                 "Bayes","$IFE_{12}$", "$IFE_5$")

tikz(file=paste(plotDIR, "bias_comp_elisa.tex", sep="/"),
     height=4, width=4, standAlone=T)
par(mar=c(3, 3, 2, 0.125), mgp=c(1.25,0.25,0), tck=-0.01)
mlplot(log(bias[-c(2)]), xlim=c(-3.5,6.5),
       xlab="accuracy ($\\mu$g/L)", cex=0.75, cex.axis=0.75, axes=F,
       adj=0.3, mgp=c(1.25,0.25,0), tck=-0.01)
text(x=rep(log(0.039),5), y=1:5, as.character(names(bias)[-2]),
     cex=0.7)
##mlplot(log(bias), xlab="log bias", xlim=c(-3,3))
axis(1, at=log(c(0.05,0.1, 0.2, 0.5, 1, 2,5, 10)),
     labels=c(0.05,0.1,0.2,0.5,1,2, 5, 10), mgp=c(1.25,0.25,0), tck=-0.01)
##axis(2, at=1:5, labels=as.character(names(bias)[-2]), las=1)
##box()
dev.off()

save(bias, file="elisaabserr.RData")
## load("PO4abserr.RData")

tikz(file=paste(plotDIR, "bias_comp.tex", sep="/"),
     height=3.5, width=3.5, standAlone=T)
par(mar=c(3, 2, 2, 0.125), mgp=c(1.25,0.25,0), tck=-0.01)
mlplot(log(abserr[1:3]), xlab="accuracy ($\\mu$g/L)", axes=F)#, xlim=c(1, 5))
##mtext("accuracy ($\\mu$g/L)", side=1, outer=T, line=-1)
axis(1, at=log(c(2.5,5, 10, 15)), labels=c(2.5,5, 10, 15))
##axis(3)
text(x=rep(log(1.75),3), y=1:3, names(abserr)[1:3],
     cex.axis=0.75, las=1)
dev.off()


### Biased versus unbiased (Figure 1)
draw_circle <- function(x0=0, y0=0, r=1, ...){
    x <- x0 + seq(-r, r, length=100)
   y1 <- y0+sqrt(r^2-(x-x0)^2)
   y2 <- y0-sqrt(r^2-(x-x0)^2)
   lines(x, y1, ...)
    lines(x, y2, ...)
    invisible()
}

tikz(file=paste(plotDIR, "biasVunbias.tex", sep="/"), height=3, width=5, standAlone=F)
par(mfrow=c(1,2), pty="s", mar=c(0,0,0,0))
plot(c(-1, 1), c(-1, 1), type="n", ann=F, axes=F)
draw_circle()
draw_circle(r=0.75)
draw_circle(r=0.5)
draw_circle(r=0.125)
points(x=runif(25, 0.1, 0.35), y=runif(25, 0.1, 0.35), col="red")
points(0,0,pch=16)
title(main="biased")

plot(c(-1, 1), c(-1, 1), type="n", ann=F, axes=F)
draw_circle()
draw_circle(r=0.75)
draw_circle(r=0.5)
draw_circle(r=0.125)
points(0,0,pch=16)
points(x=runif(25, -0.9, 0.9), y=runif(25, -0.9, 0.9), col="red")
title(main="unbiased")
dev.off()
