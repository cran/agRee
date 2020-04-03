amygdala <- matrix(c(38.1, 35.7,
                  24.7,	21.9,
                  20.4,	19.5,
                  29.9,	29.8,
                  18.3,	19.4,
                  23.3,	25,
                  32.4,	28.2,
                  28,	27.8,
                  34.3,	33.2,
                  21.1,	23.8,
                  35.8,	35),
                ncol=2, byrow=TRUE)

cerebellum <- matrix(c(10.3, 10.8,
                      11.0, 11.8,
                      13.4, 13.5,
                      12.3, 13.1,
                      13.4, 14.7, 
                      11.9, 10.9),
                    ncol=2, byrow=TRUE)

brainStem <- matrix(c(0.42, 0.46,
                 0.9,  0.75,
                 1.05, 1.2,
                 0.53, 0.63,
                 0.63, 0.54,
                 0.48, 0.74,
                 1.22, 0.44,
                 0.86, 1.18,
                 0.94, 0.73,
                 0.88, 1.18), 
               ncol=2, byrow=TRUE)


petVT <- vector("list", 3)
petVT[[1]] <- amygdala
petVT[[2]] <- cerebellum
petVT[[3]] <- brainStem

names(petVT) <- c("amygdala", "cerebellum", "brainStem")
