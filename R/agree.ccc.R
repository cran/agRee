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
    ue <- mvn.ub(X)
    hatTheta <- getBECCC(ue[[1]], ue[[2]], n)
    hatTheta.Z <- 1/2 * log((1+hatTheta)/(1-hatTheta))
    
    hatThetaPartial <- hatThetaPseudo <- rep(0, n)
    hatThetaPartial.Z <- hatThetaPseudo.Z <- rep(0, n)
    for(i in 1:n){
        ue <- mvn.ub(X[-i,])
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

ccc.mvn.mcmc <- function(X, alpha, nsim, method){
    mcmc.mvn <- mvn.bayes(X, nsim, prior=method)
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
    list(value=mean(ccc.Bayes), icl=hpd[1], icu=hpd[2])
}
#-----------------------------------------------------------------------

agree.ccc <- function(ratings, conf.level=0.95,
                      method=c("jackknifeZ", "jackknife",
                               "mvn.jeffreys", "mvn.conjugate"),
                      nsim=10000){
    if(!is.matrix(ratings) || ncol(ratings) < 2)
      stop("'ratings' has to be a matrix of at least two columns.")

    X <- ratings
    n <- nrow(X)
    

    method <- match.arg(method)

    alpha <- 1 - conf.level
    
    estimate <- switch(method,
                 mvn.jeffreys = ccc.mvn.mcmc(X, alpha, nsim, "Jeffreys"),
                 mvn.conjugate = ccc.mvn.mcmc(X, alpha, nsim, "Conjugate"),
                 jackknifeZ = unlist(ccc.nonpara.jackknife(X, alpha)[c(2,4)]),
                 jackknife = unlist(ccc.nonpara.jackknife(X, alpha)[c(1,3)]))
    list(value=estimate[[1]], lbound=estimate[[2]], ubound=estimate[[3]])
}
