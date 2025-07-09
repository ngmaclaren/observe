library(ppcor)

pID <- 27

datadir <- "~/Documents/reduction/data/TAPAS/"
datafiles <- list.files(path = datadir, pattern = ".csv")

tsdir <- "/projects/academic/naokimas/neil/brains-ns50/"
tsfiles <- list.files(path = tsdir, pattern = ".txt")

datafile <- datafiles[pID]
tsfile <- tsfiles[pID]

df <- read.table(paste0(tsdir, tsfile), sep = " ")
nr <- nrow(df)

X.train <- as.matrix(df[seq_len(nr-100), ])
X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

avgxi <- colMeans(X.test)
mu <- mean(avgxi) 
sigma <- sd(avgxi) 
threshold <- mu + 5*sigma 
test <- any(abs(avgxi) >= threshold) 
if(isFALSE(test)) flag <- "use" else flag <- "dontuse"

print(flag)

                                        # TAPAS sparse adjacency
A <- t(as.matrix(read.csv(paste0(datadir, datafile), header = FALSE)))
                                        # zero-order Pearson
Z <- cor(X.train, method = "pearson")
                                        # set diag to NA b/c not relevant for this analysis
diag(Z) <- NA
                                        # partial Pearson
P <- pcor(X.train, method = "pearson")$estimate # but which element?
                                        # same
diag(P) <- NA

print("TAPAS")
print(summary(as.numeric(A)))
print(sd(as.numeric(A), na.rm = TRUE))
print("Zero-order Pearson")
print(summary(as.numeric(Z)))
print(sd(as.numeric(Z), na.rm = TRUE))
print("Partial Pearson")
print(summary(as.numeric(P)))
print(sd(as.numeric(P), na.rm = TRUE))




