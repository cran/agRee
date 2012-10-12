#obtain the unbiased estimator of parameters from a MVN dist.
# input
#      X: samples, one row per observation
#note that the MLE of \Sigma is "(n-1)/n * hatSigma"
getUEMVN <- function(X){
    d <- dim(X)
    n <- d[1]
    p <- d[2]
    barX <- colMeans(X)
    #might be a better way to get S, say using sth provided in R
    S <- rowSums(apply(X, 1, function(x) (x - barX) %*% t(x - barX)))
    S <- matrix(S, nrow=p, ncol=p)
    list(hatMu=barX, hatSigma=S/(n-1))
}
 

#get unbiased estimator of the CCC based on CarrascoJover Biometrics 2003.
#the formula is at the end of page 851

getUECCC <- function(barY, S, n){
    if(!is.matrix(S) || nrow(S)!=ncol(S))
        stop("S has to be a square matrix.")
    if(nrow(S) != length(barY))
        stop("The leading dimenstion of 'S' has to be equal to the length of 'barY'.")
    if(! all(S[upper.tri(S)]==S[lower.tri(S)]))
        stop("'S' has to be a symmetric matrix.")
    if(! all(eigen(S)$values > 0))
        stop("'S' has to be positive definite.")

    k <- nrow(S)
    #numerator
    nu <- 2 * sum(S[upper.tri(S)])
    #denominator
    d1 <- (k-1) * sum(diag(S))
    d2 <- sum(outer(barY, barY, "-")^2) / 2
    d3 <- (k*(k-1)/n) * (d1 - nu)
    
    nu / (d1+ d2 - d3)
        
}

#Get biased estimator of the CCC. Since \sum\sum(\bar{Y_i}-\bar{Y_j})^2 is a
#biased estimator of \sum\sum(\bar{\mu_i}-\bar{\mu_j})^2
#based on CarrascoJover Biometrics 2003

#NOTE: the formula is different than that based on CarrascoJover
#      Biometrics 2003.at the end of page 851.
getBECCC <- function(barY, S, n){
    if(!is.matrix(S) || nrow(S)!=ncol(S))
        stop("S has to be a square matrix.")
    if(nrow(S) != length(barY))
        stop("The leading dimenstion of 'S' has to be equal to the length of 'barY'.")
    if(! all(S[upper.tri(S)]==S[lower.tri(S)]))
        stop("'S' has to be a symmetric matrix.")
    if(! all(eigen(S)$values > 0))
        stop("'S' has to be positive definite.")

    k <- nrow(S)
    #numerator
    nu <- 2 * sum(S[upper.tri(S)])
    #denominator
    d1 <- (k-1) * sum(diag(S))
    d2 <- sum(outer(barY, barY, "-")^2) / 2
    
    nu / (d1+ d2)
        
}

#leave one out jackknife based on biased estimator obtained using 'getBECCC()'
ccc.nonpara.jackknife <- function(X, alpha){
    d <- dim(X)
    n <- d[1]
    k <- d[2]
    ue <- getUEMVN(X)
    hatTheta <- getBECCC(ue[[1]], ue[[2]], n)
    hatTheta.Z <- 1/2 * log((1+hatTheta)/(1-hatTheta))
    
    hatThetaPartial <- hatThetaPseudo <- rep(0, n)
    hatThetaPartial.Z <- hatThetaPseudo.Z <- rep(0, n)
    for(i in 1:n){
        ue <- getUEMVN(X[-i,])
        hatThetaPartial[i] <- getBECCC(ue[[1]], ue[[2]], n-1)
        hatThetaPseudo[i] <- n*hatTheta - (n-1)*hatThetaPartial[i]

        hatThetaPartial.Z[i] <- 1/2 *
          log((1+hatThetaPartial[i])/(1-hatThetaPartial[i]))
        hatThetaPseudo.Z[i] <- n*hatTheta.Z - (n-1)*hatThetaPartial.Z[i]
    }
    point <- mean(hatThetaPseudo)
    point.Z <- mean(hatThetaPseudo.Z)
    
    ST2 <- mean((hatThetaPseudo - point)^2) / (n-1)
    ST <- sqrt(ST2)
    #z.alpha2 <- qnorm(1 - alpha/2)
    z.alpha2 <- qt(1 - alpha/2, df=n-1)
    CI <- c(point - z.alpha2*ST, point + z.alpha2*ST)
    
    ST2.Z <- mean((hatThetaPseudo.Z - point.Z)^2) / (n-1)
    ST.Z <- sqrt(ST2.Z)
    CI.Z <- c(point.Z - z.alpha2*ST.Z, point.Z + z.alpha2*ST.Z)
    CI.Z <- c((exp(2*CI.Z[1])-1) / (exp(2*CI.Z[1])+1),
              (exp(2*CI.Z[2])-1) / (exp(2*CI.Z[2])+1))

    list(CI=CI, CI.Z=CI.Z)
    
}

#-----------------------------------------------------------------------------

#page 87 of GelmanCarlinSternRubin
simMvnConjugate <- function(X, nsim){
    
    n <- nrow(X)
    k <- ncol(X)
    
    bar.x <- colMeans(X)
    #might be a better way to get S, say using sth provided in R
    S <- rowSums(apply(X, 1, function(x) (x - bar.x) %*% t(x - bar.x)))
    S <- matrix(S, nrow=k, ncol=k)

    mu0 <- bar.x
    
    k0 <- 1
    
    #v0 <- k + 2
    v0 <- k

    #it is not a good idea to let Gamma0 <- S
    #Gamma0 <- S
    Gamma0 <- S / (n-1)


    kn <- k0 + n
    mun <- (k0/kn) * mu0  +  (n/kn) * bar.x
    vn <- v0 + n
    Gamman <- Gamma0 + S + (k0*n / kn) * (bar.x-mu0) %*% t(bar.x-mu0)

    Sigma <- array(0, dim=c(k, k, nsim))
    Mu <- matrix(0, ncol=k, nrow=nsim)

    
    for(i in 1:nsim){
        #browser()
        Sigma[,,i] <- riwish(vn, Gamman)
        Mu[i,] <- mvrnorm(1, mun, Sigma[,,i]/kn)
    }
    list(Sigma=Sigma, Mu=Mu)
}

ccc.mcmc.conjugate <- function(X, alpha, nsim){
    mcmc.mvn <- simMvnConjugate(X, nsim)
    Sigma <- mcmc.mvn$Sigma
    Mu <- mcmc.mvn$Mu
    k <- ncol(X)
    ccc.Bayes <- rep(0, nsim)
    for(iter in 1:nsim){
        Sig <- Sigma[,,iter]
        M <- Mu[iter,]
        ss <- 0
        sm <- 0
        for(i in 1:(k-1)){
            for(j in (i+1):k){
                ss <- ss + Sig[i,j]
                sm <- sm + (M[i] - M[j])^2
            }
        }
        ccc.Bayes[iter] <- 2*ss / ((k-1)*sum(diag(Sig)) + sm)
    }
    ccc.Bayes <- mcmc(ccc.Bayes)
    hpd <- HPDinterval(ccc.Bayes, prob=1-alpha)
    list(icl=hpd[1], icu=hpd[2])
}


#----------------------------------------------------------------------
simBvnJeffreys <- function(X, nsim){
    a <- 1
    b <- 0

    sim.save <- matrix(0, nrow=nsim, ncol=5)
    colnames(sim.save) <- c("mu1", "mu_2", "sigma1", "sigma2", "rho")
    n <- nrow(X)
    bar.x <- colMeans(X)
    bar.x1 <- bar.x[1] 
    bar.x2 <- bar.x[2]
    #might be a better way to get S, say using sth provided in R
    S <- rowSums(apply(X, 1, function(x) (x - bar.x) %*% t(x - bar.x)))
    S <- matrix(S, nrow=2, ncol=2)
    s11 <- S[1,1]
    s22 <- S[2,2]
    r <- S[1,2]/sqrt(s11*s22)

   
    Z1 <- rnorm(nsim)
    Z2 <- rnorm(nsim)
    Z3 <- rnorm(nsim)
    chisq.na <- rchisq(nsim, n-a)
    chisq.nb <- rchisq(nsim, n-b)
 
    mu1 <- bar.x1 + 
           Z1/sqrt(chisq.na) * sqrt(s11/n)

    mu2 <- bar.x2 + 
           Z1/sqrt(chisq.na) * r*sqrt(s22/n) + 
           ( Z2/sqrt(chisq.nb) - 
             Z3/sqrt(chisq.nb) * Z1/sqrt(chisq.na)) * 
             sqrt(s22*(1-r^2)/n)

    sigma1 <- sqrt(s11/chisq.na)

    sigma2 <- sqrt(s22*(1-r^2)) *
              sqrt(1/chisq.nb + 
                   1/chisq.na * (Z3/sqrt(chisq.nb) - r/sqrt(1-r^2))^2)
   
    Y <- -Z3/sqrt(chisq.na) + 
             sqrt(chisq.nb)/sqrt(chisq.na) * r/sqrt(1-r^2)
    rho <- Y / sqrt(1+Y^2)

    sim.save <- cbind(mu1, mu2, sigma1, sigma2, rho)
    
    sim.save
}

ccc.jeffreys.bi <- function(X, alpha, nsim){
    mcmc.bvn <- simBvnJeffreys(X, nsim)
    rho <- mcmc.bvn[,5]
    sigma1 <- mcmc.bvn[,3]
    sigma2 <- mcmc.bvn[,4]
    mu1 <- mcmc.bvn[,1]
    mu2 <- mcmc.bvn[,2]
    
    ccc.Bayes <- 2 * rho * sigma1  * sigma2 / 
      (sigma1^2 + sigma2^2 + (mu1-mu2)^2)
    ccc.Bayes <- mcmc(ccc.Bayes)
    hpd <- HPDinterval(ccc.Bayes, prob=1-alpha)
    list(icl=hpd[1], icu=hpd[2])
}

simMvnJeffreys <- function(X, nsim){
    n <- nrow(X)
    k <- ncol(X)
    bar.x <- colMeans(X)
    #might be a better way to get S, say using sth provided in R
    S <- rowSums(apply(X, 1, function(x) (x - bar.x) %*% t(x - bar.x)))
    S <- matrix(S, nrow=k, ncol=k)

    Sigma <- array(0, dim=c(k, k, nsim))
    Mu <- matrix(0, ncol=k, nrow=nsim)
    for(i in 1:nsim){
        #browser()
        Sigma[,,i] <- riwish(n, S)
        Mu[i,] <- mvrnorm(1, bar.x, Sigma[,,i]/n)
    }
    list(Sigma=Sigma, Mu=Mu)
}

ccc.jeffreys.multi <- function(X, alpha, nsim){
    mcmc.mvn <- simMvnJeffreys(X, nsim)
    Sigma <- mcmc.mvn$Sigma
    Mu <- mcmc.mvn$Mu
    k <- ncol(X)
    ccc.Bayes <- rep(0, nsim)
    for(iter in 1:nsim){
        Sig <- Sigma[,,iter]
        M <- Mu[iter,]
        ss <- 0
        sm <- 0
        for(i in 1:(k-1)){
            for(j in (i+1):k){
                ss <- ss + Sig[i,j]
                sm <- sm + (M[i] - M[j])^2
            }
        }
        ccc.Bayes[iter] <- 2*ss / ((k-1)*sum(diag(Sig)) + sm)
    }
    ccc.Bayes <- mcmc(ccc.Bayes)
    hpd <- HPDinterval(ccc.Bayes, prob=1-alpha)
    list(icl=hpd[1], icu=hpd[2])
    
}

ccc.mcmc.jeffreys <- function(X, alpha, nsim){
    if(ncol(X)==2)
        ccc.jeffreys.bi(X, alpha, nsim)
    else
        ccc.jeffreys.multi(X, alpha, nsim)
}

#-----------------------------------------------------------------------
agree.ccc <- function(ratings, conf.level=0.95,
                      method=c("jeffreys", "conjugate", "jackknifeZ", "jackknife"),
                      nsim=10000){
    if(!is.matrix(ratings) || ncol(ratings) < 2)
      stop("'ratings' has to be a matrix of at least two columns.")

    X <- ratings
    n <- nrow(X)

    ue <- getUEMVN(X)
    value <- getUECCC(ue[[1]], ue[[2]], n)

    method <- match.arg(method)

    alpha <- 1 - conf.level
    
    CI <- switch(method,
                 jeffreys = ccc.mcmc.jeffreys(X, alpha, nsim),
                 conjugate = ccc.mcmc.conjugate(X, alpha, nsim),
                 jackknifeZ = ccc.nonpara.jackknife(X, alpha)$CI.Z,
                 jackknife = ccc.nonpara.jackknife(X, alpha)$CI)
    list(value=value, lbound=CI[[1]], ubound=CI[[2]])
}
