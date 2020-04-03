lesionBurden.M <- matrix(c(20, 26.8, 9.6,
                           21.2, 26.5, 10.5,
                           20.8, 22.5, 10.6,
                           20.6, 23.1, 9.2,
                           20.2, 24.3, 10.4,
                           19.1, 24.1, 10.4,
                           21, 26, 10.1,
                           20.4, 26.8, 8,
                           19.2, 24.9, 10.1,
                           19.2, 27.7, 8.9),
                         nrow=3)

lesionBurden.G <- matrix(c(19.5, 22.1, 8.5,
                           19.5, 21.9, 8.5,
                           19.6, 22, 8.3,
                           19.7, 22.1, 8.3,
                           19.3, 21.9, 8.3,
                           19.1, 21.8, 8,
                           19.1, 21.7, 8,
                           19.3, 21.7, 8,
                           19.2, 21.7, 8,
                           19.5, 21.8, 8.1),
                           nrow=3)
lesionBurden <- array(0, dim=c(3, 10, 2))
lesionBurden[,,1] <- lesionBurden.M
lesionBurden[,,2] <- lesionBurden.G
dimnames(lesionBurden)[[3]] <- c("lesionBurden.M", "lesionBurden.G")
