## A hockey stick model

stan_model2 <- "
          data{
            int N; //the number of observations
            vector[N] y; //the response
            vector[N] x;
            int R;
            int region[N];

            real theta;
            real beta1;
          }
          parameters {
            real beta0; //the regression parameters
            real<lower=0> delta;
            real phi; //change point

            vector[R] re0;
            vector[R] reP;
            vector[R] reD;
            real<lower=0> sigma;
            real<lower=0> sigma0;
            real<lower=0> sigmaD;
            real<lower=0> sigmaP;
          }
          transformed parameters {
            vector[N] mu;
            for (i in 1:N)
              mu[i] = (beta0+re0[region[i]]) +
                      (beta1) * (x[i]-(phi+reP[region[i]])) +
                      (delta+reD[region[i]]) * theta *
              log1p(exp((x[i]-(phi+reP[region[i]]))/theta));
          }
          model {
            phi ~ cauchy(0,1);
            beta0 ~ normal(0,5);
            delta ~ normal(0,5);
            sigma ~ cauchy(0, 1);
            sigma0 ~ cauchy(0, 1);
            sigmaD ~ cauchy(0, 1);
            sigmaP ~ cauchy(0, 1);

            re0 ~ normal(0, sigma0);
            reD ~ normal(0, sigmaD);
            reP ~ normal(0, sigmaP);

            y ~ normal(mu, sigma);
          }
          generated quantities {
            real Ph;
            vector[R] delPh;
            real B0;
            vector[R] delB0;
            real De;
            vector[R] delD;
            Ph = phi + mean(reP[]);
            B0 = beta0 + mean(re0[]);
            De = delta + mean(reD[]);
            for (i in 1:R){
              delPh[i] = reP[i] - mean(reP[]);
              delB0[i] = re0[i] - mean(re0[]);
              delD[i]  = reD[i] - mean(reD[]);
            }
          }
"

stan.in <- function(infile=eriedata, x = "POC", grp="Year",
                    n.chains=nchains){
  infile$part_microcystin[infile$part_microcystin<0.01]<-0.01
  x <- log(infile[,x])
  y <- log(infile$part_microcystin)
  gr <- as.numeric(ordered(infile[, grp]))
  n <- dim(infile)[1]
  R <- max(gr)
  xlimit <- range(x, na.rm=T)

  inits <- list()
  bugs.data <- list(N=n,  R=R, y=y, x=x, region=gr,
                    theta=0.01*diff(xlimit),
                    phi_low=xlimit[1], phi_up=xlimit[2],
                    beta1=0)
  for (i in 1:n.chains)
    inits[[i]] <- list(beta0=rnorm(1), delta=runif(1),
                       phi=runif(1, xlimit[1], xlimit[2]),
                       re0=rep(0, R), reD=rep(0,R),
                       reP=rep(0,R),
                       sigma=runif(1),
                       sigma0=runif(1),sigmaD=runif(1),
                       sigmaP=runif(1))

  para <- c("B0", "De", "Ph", "delB0", "delD", "delPh",
            "sigma", "sigma0","sigmaD","sigmaP")
  return(list(para=para, data=bugs.data,
              inits=inits,n.chains=n.chains,
              theta=0.01*diff(xlimit)))
}

running_fun <- function(Data, stanfit=fit){
  input.to.bugs <- stan.in(infile=Data, n.chains=nchains,
                           x="Chla")
  fit2keep <- sampling(stanfit, data = input.to.bugs$data,
                       init=input.to.bugs$inits,
                       pars = input.to.bugs$para,
                       iter=niters, thin=nthin,
                       chains=input.to.bugs$n.chains,
                       control = list(adapt_delta = 0.99,
                                      max_treedepth=15))
  return(list(fit=fit2keep, input=input.to.bugs))
}
