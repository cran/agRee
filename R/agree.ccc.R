ccc.mm <- function(barY, S, n){
    if(!is.matrix(S) || nrow(S)!=ncol(S))
        stop("S has to be a square matrix.")
    if(nrow(S) != length(barY))
        stop("The leading dimenstion of 'S' has to be equal to the length of 'barY'.")

    if(! isSymmetric(S))
        stop("'S' has to be a symmetric matrix.")
    if(! all(eigen(S)$values > 0))
        stop("'S' has to be positive definite.")

    S <- S * (n-1) / n
    k <- nrow(S)
    #numerator
    nu <- 2 * sum(S[upper.tri(S)])
    #denominator
    d1 <- (k-1) * sum(diag(S))
    d2 <- sum(outer(barY, barY, "-")^2) / 2
    
    nu / (d1+ d2)
        
}


#leave one out jackknife based on biased estimator obtained using 'ccc.mm()'
ccc.nonpara.jackknife <- function(X, alpha){
    d <- dim(X)
    n <- d[1]
    k <- d[2]
    ue <- mvn.ub(X)
    hatTheta <- ccc.mm(ue[[1]], ue[[2]], n)
    hatTheta.Z <- 1/2 * log((1+hatTheta)/(1-hatTheta))
    
    hatThetaPartial <- hatThetaPseudo <- rep(0, n)
    hatThetaPartial.Z <- hatThetaPseudo.Z <- rep(0, n)
    for(i in 1:n){
        ue <- mvn.ub(X[-i,])
        hatThetaPartial[i] <- ccc.mm(ue[[1]], ue[[2]], n-1)
        hatThetaPseudo[i] <- n*hatTheta - (n-1)*hatThetaPartial[i]

        hatThetaPartial.Z[i] <- 1/2 *
          log((1+hatThetaPartial[i])/(1-hatThetaPartial[i]))
        hatThetaPseudo.Z[i] <- n*hatTheta.Z - (n-1)*hatThetaPartial.Z[i]
    }
    point <- mean(hatThetaPseudo)
    point.Z <- mean(hatThetaPseudo.Z)
    
    ST2 <- mean((hatThetaPseudo - point)^2) / (n-1)
    ST <- sqrt(ST2)
    z.alpha2 <- qt(1 - alpha/2, df=n-1)
    CI <- c(point - z.alpha2*ST, point + z.alpha2*ST)
    
    ST2.Z <- mean((hatThetaPseudo.Z - point.Z)^2) / (n-1)
    ST.Z <- sqrt(ST2.Z)
    CI.Z <- c(point.Z - z.alpha2*ST.Z, point.Z + z.alpha2*ST.Z)
    CI.Z <- c((exp(2*CI.Z[1])-1) / (exp(2*CI.Z[1])+1),
              (exp(2*CI.Z[2])-1) / (exp(2*CI.Z[2])+1))

    list(point=point, point.Z=(exp(2*point.Z[1])-1) / (exp(2*point.Z[1])+1),
         CI=CI, CI.Z=CI.Z)
    
}

#-----------------------------------------------------------------------------
#bootstrap method
ccc.boot.fun <- function(dd, i){
    X <- dd[i,]
    n <- nrow(X)
    ue <- mvn.ub(X)
    ccc.mm(ue[[1]], ue[[2]], n)
}

ccc.nonpara.bootstrap <- function(X, nboot, alpha){
    n <- nrow(X)
    point <- ccc.boot.fun(X, 1:n)
    index.array <- matrix(sample.int(n, size = n*nboot, replace = TRUE),
                          nboot, n)
    ccc.boot <- lapply(1:nboot, function(i) ccc.boot.fun(X, index.array[i,]))
    CI <- as.vector(quantile(unlist(ccc.boot), c(alpha/2, 1-alpha/2), type=6))
    list(point=point, icl=CI[1], icu=CI[2])
}

#-----------------------------------------------------------------------------
ccc.mvst.mcmc <- function(X, alpha, nmcmc=10000){

    mcmc.mvst <- mvst.mcmc(X, nmcmc=nmcmc*1.1)
    Mu <- mcmc.mvst$Mu
    D <- mcmc.mvst$Delta
    Sigma <- mcmc.mvst$Sigma
    nu <- mcmc.mvst$nu
    #DIC <- mcmc.mvst$DIC
    
    nmcmc <- nrow(Mu) 
    k <- ncol(X)
    ccc.Bayes <- rep(0, nmcmc)
    for(iter in 1:nmcmc){
        M <- Mu[iter,] + (nu[iter]/pi)^(1/2)*gamma((nu[iter]-1)/2) / gamma(nu[iter]/2) * D[iter,]
        DD <- D[iter,] %*% t(D[iter,])
        Sig <- (Sigma[iter,,] + diag(D[iter,])^2 + 2/pi*(DD-diag(D[iter,])^2)) * (nu[iter]/(nu[iter]-2)) -
               nu[iter] / pi * (gamma((nu[iter]-1)/2) / gamma(nu[iter]/2))^2 * DD

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
    hpd <- as.vector(HPDinterval(ccc.Bayes, 1-alpha))

    list(point=median(ccc.Bayes), icl=hpd[1], icu=hpd[2])

}
#-----------------------------------------------------------------------------
ccc.mvsn.mcmc <- function(X, alpha, nmcmc){

    mcmc.mvsn <- mvsn.mcmc(X, nmcmc=nmcmc*1.1)
    Mu <- mcmc.mvsn$Mu
    Delta <- mcmc.mvsn$Delta
    Sigma <- mcmc.mvsn$Sigma
    #DIC <- mcmc.mvsn$DIC

    nmcmc <- nrow(Mu) 
    k <- ncol(X)
    ccc.Bayes <- rep(0, nmcmc)
    for(iter in 1:nmcmc){
        Sig <- Sigma[iter,,] + (1-2/pi)*diag(Delta[iter,])^2
        M <- Mu[iter,] + sqrt(2/pi)*Delta[iter,]
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
    hpd <- as.vector(HPDinterval(ccc.Bayes, 1-alpha))

    list(point=median(ccc.Bayes), icl=hpd[1], icu=hpd[2])

}


#-----------------------------------------------------------------------------
ccc.lognormalNormal.mcmc <- function(X, alpha, nmcmc){
    mcmc.lognormalNormal <- lognormalNormal.mcmc(X=X, nmcmc=nmcmc)
    mu <- mcmc.lognormalNormal$mu
    sigma2.alpha <- mcmc.lognormalNormal$sigma2.alpha
    sigma2.e <- mcmc.lognormalNormal$sigma2.e
    betas <- mcmc.lognormalNormal$betas

    nrater <- ncol(X)
    diff.betas <- rep(0, length(mu))
    for(i in 1:(nrater-1)){
        for(j in 2:nrater){
            diff.betas <- (betas[,j] - betas[,i])^2 + diff.betas
        }
    }
    
    ccc.Bayes <- (exp(2*mu+sigma2.alpha) * (exp(sigma2.alpha)-1)) /
      (exp(2*mu+sigma2.alpha) * (exp(sigma2.alpha)-1) +
       diff.betas/(nrater *(nrater-1)) + sigma2.e)
    
    point <- median(ccc.Bayes, na.rm=TRUE)
    CI <- as.vector(quantile(ccc.Bayes, c(alpha/2, 1-alpha/2), na.rm=TRUE))
    list(point=point, icl=CI[1], icu=CI[2])

}
#-----------------------------------------------------------------------------
ccc.mvt.mcmc <- function(X, alpha, nmcmc,
                         prior.lower.v, prior.upper.v,
                         prior.Mu0, prior.Sigma0,
                         prior.p, prior.V,
                         initial.v, initial.Sigma){

    mcmc.mvt <- mvt.mcmc(X, prior.lower.v, prior.upper.v,
                         prior.Mu0, prior.Sigma0,
                         prior.p, prior.V,
                         initial.v, initial.Sigma,
			 nmcmc)
    Sigma <- mcmc.mvt$Sigma.save
    Mu <- mcmc.mvt$Mu.save
    v <- mcmc.mvt$v.save
    k <- ncol(X)
    ccc.Bayes <- rep(0, nmcmc)
    for(iter in 1:nmcmc){
        Sig <- Sigma[,,iter]*v[iter]/(v[iter]-2)
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
    list(point=median(ccc.Bayes), icl=hpd[1], icu=hpd[2])
}
#-----------------------------------------------------------------------------

ccc.mvn.mcmc <- function(X, alpha, nmcmc, method){
    mcmc.mvn <- mvn.bayes(X, nmcmc, prior=method)
    Sigma <- mcmc.mvn$Sigma.save
    Mu <- mcmc.mvn$Mu.save
    k <- ncol(X)
    ccc.Bayes <- rep(0, nmcmc)
    for(iter in 1:nmcmc){
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
    list(value=median(ccc.Bayes), icl=hpd[1], icu=hpd[2])
}
#-----------------------------------------------------------------------

agree.ccc <- function(ratings, conf.level=0.95,
                      method=c("jackknifeZ", "jackknife",
                               "bootstrap", "bootstrapBC",
                               "mvn.jeffreys", "mvn.conjugate",
                               "mvt", "lognormalNormal", "mvsn", "mvst"),
                      nboot=999, nmcmc=10000,
                      mvt.para=list(prior=list(lower.v=4, upper.v=25,
                                               Mu0=rep(0, ncol(ratings)),
                                               Sigma0=diag(10000, ncol(ratings)),
                                               p=ncol(ratings),
                                               V=diag(1, ncol(ratings))),
                                    initial=list(v=NULL, Sigma=NULL)),
                      NAaction=c("fail", "omit")){

    if(!is.matrix(ratings) || ncol(ratings) < 2|| nrow(ratings) < 3)
      stop("'ratings' has to be a matrix of at least two columns and three rows.")

    
    na <- match.arg(NAaction)
    ratings <- switch(na,
                      fail = na.fail(ratings),
                      omit = na.omit(ratings))
    if(!is.matrix(ratings) || ncol(ratings) < 2|| nrow(ratings) < 3)
      stop("'ratings' has to be a matrix of at least two columns and three rows after removing missing values.")

    X <- ratings
    n <- nrow(X)

    method <- match.arg(method)

    alpha <- 1 - conf.level
    
    estimate <- switch(method,
                       mvsn = ccc.mvsn.mcmc(X, alpha, nmcmc),
                       mvst = ccc.mvst.mcmc(X, alpha, nmcmc),
                       lognormalNormal = ccc.lognormalNormal.mcmc(X, alpha, nmcmc),
                       mvt = ccc.mvt.mcmc(X, alpha, nmcmc,
                           mvt.para$prior$lower.v, mvt.para$prior$upper.v,
                           mvt.para$prior$Mu0, mvt.para$prior$Sigma0,
                           mvt.para$prior$p, mvt.para$prior$V,
                           mvt.para$initial$v, mvt.para$initial$Sigma),
                       mvn.jeffreys = ccc.mvn.mcmc(X, alpha, nmcmc, "Jeffreys"),
                       mvn.conjugate = ccc.mvn.mcmc(X, alpha, nmcmc, "Conjugate"),
                       jackknifeZ = unlist(ccc.nonpara.jackknife(X, alpha)[c(2,4)]),
                       jackknife = unlist(ccc.nonpara.jackknife(X, alpha)[c(1,3)]),
                       bootstrap = ccc.nonpara.bootstrap(X, nboot, alpha),
                       bootstrapBC = sapply(ccc.nonpara.bootstrap(X, nboot, alpha), function(i) i*(1+10/n)))
    
    list(value=estimate[[1]], lbound=estimate[[2]], ubound=estimate[[3]])
}
