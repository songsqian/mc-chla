## creating a subset of data for a given year
## Grouping sampling events to at least n=n_min
erie_sub <- function(file=eriedata, Yr=2018, n_min = 8,
	       x="Chla", y="part_microcystin"){
  subdata <- file[file$Year==Yr & !is.na(file[,y]),]
  cs <- cumsum(sz <- table(subdata$Rdate))
  stepcs <- 0
  k <- 1
  gr <- rep(1, length(cs))

  for (i in 1:length(cs)){
    if (stepcs >= n_min){
	    k <- k+1
	    stepcs <- 0
	  }
	  stepcs <- stepcs + sz[i]
	  gr[i] <- k
  }
  subdata$gr <- gr[as.numeric(ordered(subdata$Rdate))]
  return(subdata)
}

prior <- function(fitrv, b0="B0", de="De", ph="Ph",
                  s0="sigma0", sD="sigmaD",sP="sigmaP",
		              n0=20, setn0=F){
	## fitrv: an rv object of stan fitted model
	## setn0: whether to use non-informative prior

	tmp <- summary(fitrv[names(fitrv)==b0])
	Eb0 <- tmp$mean
	Vb0 <- tmp$sd^2

	tmp <- summary(fitrv[names(fitrv)==de])
	EDe <- tmp$mean
	VDe <- tmp$sd^2

	tmp <- summary(fitrv[names(fitrv)==ph])
	EPh <- tmp$mean
	VPh <- tmp$sd^2

	tmp <- summary(fitrv[names(fitrv)==s0]^2)
	Esigma0 <- tmp$mean
	Vsigma0 <- tmp$sd^2

	tmp <- summary(fitrv[names(fitrv)==sD]^2)
	EsigmaD <- tmp$mean
	VsigmaD <- tmp$sd^2

	tmp <- summary(fitrv[names(fitrv)==sP]^2)
	EsigmaP <- tmp$mean
	VsigmaP <- tmp$sd^2

	if (setn0) {
	    alpha0 <- n0+1
	    alphaD <- n0+1
	    alphaP <- n0+1
	} else {
	    alpha0 <- 2+Esigma0^2/Vsigma0
	    alphaD <- 2+Esigma0^2/VsigmaD
	    alphaP <- 2+Esigma0^2/VsigmaP
	}
	beta0 <- Esigma0*(alpha0-1)
	betaD <- EsigmaD*(alphaD-1)
	betaP <- EsigmaP*(alphaP-1)
	lambda0 <- Esigma0/Vb0
	lambdaD <- EsigmaD/VDe
	lambdaP <- EsigmaP/VPh
	## limiting alpha+beta < 1000
	while (alpha0+beta0 > 1000){
	  alpha0 <- alpha0/10
	  beta0 <- beta0/10
	}
	while (alphaD+betaD > 1000){
	  alphaD <- alphaD/10
	  betaD <- betaD/10
	}
	while (alphaP+betaP > 1000){
	  alphaP <- alphaP/10
	  betaP <- betaP/10
	}

	return(list(m0=Eb0,  mD=EDe, mP=EPh,
		    lmbd0=lambda0, lmbdD=lambdaD, lmbdP=lambdaP,
		    al0=alpha0, alP=alphaP, alD=alphaD,
		    bt0=beta0, btP=betaP, btD=betaD))
}


## Stan model (reporting $\sigma$)
stan_model3 <- "
	      data{
	      int N; //the number of observations
	      vector[N] y; //the response
	      vector[N] x;

	      real theta;
	      //real beta1;

	      real m0;
	      real mD;
	      real mP;

	      real lmbd0;
	      real lmbdD;
	      real lmbdP;

	      real al0;
	      real alP;
	      real alD;

	      real bt0;
	      real btP;
	      real btD;

	    }
	    parameters {
	      real beta0; //the regression parameters
	      real<lower=0> delta;
	      real phi; //change point

	      real<lower=0> sigma;

	      real mu0;
	      real muD;
	      real muP;

	      real<lower=0> sigma0sq;
	      real<lower=0> sigmaDsq;
	      real<lower=0> sigmaPsq;
	    }
	    transformed parameters {
	      real<lower=0> sigma0;
	      real<lower=0> sigmaD;
	      real<lower=0> sigmaP;
	      vector[N] mu;

	      sigma0 = sqrt(sigma0sq);
	      sigmaD = sqrt(sigmaDsq);
	      sigmaP = sqrt(sigmaPsq);
	      for (i in 1:N)
		mu[i] = beta0 + //beta1 * (x[i]-phi) +
                        delta * theta * log1p(exp((x[i]-phi)/theta));
	    }
	    model {
	      sigma ~ cauchy(0, 1);
	      sigma0sq ~ inv_gamma(al0, bt0);
	      sigmaDsq ~ inv_gamma(alD, btD);
	      sigmaPsq ~ inv_gamma(alP, btP);

	      mu0 ~ normal(m0, sqrt(sigma0sq/lmbd0));
	      muD ~ normal(mD, sqrt(sigmaDsq/lmbdD));
	      muP ~ normal(mP, sqrt(sigmaPsq/lmbdP));

	      phi ~ normal(muP, sigmaP);
	      beta0 ~ normal(mu0, sigma0);
	      delta ~ normal(muD, sigmaD);

	      y ~ normal(mu, sigma);
	    }
"

stan.in <- function(infile, x="Chla", y="part_microcystin",
              			n.chains=nchains, grp=NULL,
			              stdz=T, info=T, prrs = NULL){
	if (info & is.null(prrs)) stop("Need informative priors")
	if (!is.null(grp)) infile=infile[grp,]
	keep <-  (infile[,x] > 0) & (infile[,y] >0)
	infile <- infile[keep & !is.na(keep),]
	x <- log(infile[,x])
	xmu <- mean(x)
	xsd <- sd(x)
	if (stdz) x <- (x - xmu)/xsd
	y <- log(infile[,y])
	n <- dim(infile)[1]
	if (info){
	    m0 = prrs$m0
	    mD = prrs$mD
	    mP = prrs$mP
	    lmbd0=prrs$lmbd0
	    lmbdD=prrs$lmbdD
	    lmbdP=prrs$lmbdP
	    al0=prrs$al0
	    alP=prrs$alP
	    alD=prrs$alD
	    bt0=prrs$bt0
	    btP=prrs$btP
	    btD=prrs$btD
	}else{
	    m0 = 0
	    mD = 0
	    mP = 0
	    lmbd0=1
	    lmbdD=1
	    lmbdP=1
	    al0=2
	    alP=2
	    alD=2
	    bt0=2
	    btP=2
	    btD=2
	}

	s0 <- sqrt(bt0/(al0-1))
	sD <- sqrt(btD/(alD-1))
	sP <- sqrt(btP/(alP-1))

	inits <- list()
	if (stdz) theta <- 0.04
	else theta <- 0.01*diff(range(x))
	bugs.data <- list(N=n, y=y, x=x,
			  theta=theta, #beta1=0,
			  m0 = m0,  mD = mD, mP = mP,
			  lmbd0=lmbd0, lmbdD=lmbdD, lmbdP=lmbdP,
			  al0=al0, alP=alP, alD=alD,
			  bt0=bt0, btP=btP, btD=btD )
	for (i in 1:n.chains)
	    inits[[i]] <- list(beta0=rnorm(1, m0, s0),
			       delta=abs(rnorm(1,mD,sD)),
			       phi=runif(1, range(x)[1], range(x)[2]),
			       sigma=runif(1), sigmaPsq=runif(1),
			       sigmaDsq=runif(1),
			       sigma0sq=runif(1),
			       mu0=rnorm(1, m0,s0),
			       muD=abs(rnorm(1, mD,sD)),
			       muP=rnorm(1, mP,sP))
	para <- c("beta0", "delta", "phi","sigma",
		        "mu0", "muD","muP", "sigmaP","sigma0", "sigmaD")
	return(list(para=para, data=bugs.data,
		          inits=inits,n.chains=n.chains,
		          mux=xmu, sdx=xsd, theta=theta))
}

