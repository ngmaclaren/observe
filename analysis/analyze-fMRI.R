                                        # 2025-05-18: Separate network creation?
                                        # I don't think it's necessary because we need X.test and computing the
                                        # network is comparatively trivial.
                                        # 2025-05-18: Remove old (commented) code
                                        # 2025-05-18: Do we want to do an ANOVA with this data, like for home/other?
                                        # 2025-05-18: Check the downstream code. What objects do I need?

library(parallel)
ncores <- detectCores()-2
library(igraph)
library(sdn)
library(optNS)

datadir <- "/projects/academic/naokimas/neil/brains-ns300/"
datafiles <- list.files(path = datadir, pattern = ".txt")

                                        # for testing
## datafile <- datafiles[10]
## datafile <- "./test/101006-50.txt"

make_nodesets <- function(datafile) {
    calc_obj <- function(z, y) sum((z - y)^2)/length(y)
    obj_fn <- function (vs, y, Y, optimize_weights = FALSE, ws = NULL, ...) {
        Z <- as.matrix(Y[, vs])
        if (optimize_weights) {
            if (is.null(ws)) 
                ws <- quadoptm(vs, y, Y)
            z <- apply(Z, 1, weighted.mean, ws)
        }
        else {
            z <- rowMeans(Z)
        }
        calc_obj(z, y)
    }

    df <- read.table(paste0(datadir, datafile), sep = " ")
    nr <- nrow(df)
    pos <- which(datafiles %in% datafile)

    X.train <- as.matrix(df[seq_len(nr-100), ])
    X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

                                        # 2025-05-18: Should I check X.test or X.train?
                                        # I guess it doesn't matter, but we said we used X.test so this is ok.
    avgxi <- colMeans(X.test)
    mu <- mean(avgxi)
    sigma <- sd(avgxi)
    threshold <- mu + 5*sigma # 4
    test <- any(abs(avgxi) >= threshold)
                                        # 2025-05-18: Take the opportunity to remove this silliness. The flag
                                        # should be stored in the data list, not as a separate object.
    if(isFALSE(test)) flag <- "use" else flag <- "dontuse"
                                        # uncomment here to return only flag
                                        # (e.g., creating check5.rds
    ##return(flag)
##}
    Cmat <- cor(X.train, method = "pearson")
    A <- Cmat
    A[which(A < 0, arr.ind = TRUE)] <- 0

    g <- graph_from_adjacency_matrix(A, "undirected", weighted = TRUE, diag = FALSE)
                                        # uncomment here to return the networks
##     return(g)
## }

    N <- vcount(g)
    n <- floor(log(N))

    times <- 0:15
    control <- list(times = times, ncores = ncores)
                                        # 2025-05-18: here. this needs to accept a shell arg
    params <- c(.doublewell, list(A = A)) 
    params$use.matrix <- TRUE
    Ds <- seq(0, 1, length.out = 100)
                                        # 2025-05-18: also here
    Y <- solve_in_range(Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode") 
    y <- rowMeans(Y)

                                        # 2025-05-18: include degree-preserving here ↓
    opts <- make_dataset(
        ntrials = 100, ns.type = "opt", ncores = ncores, n = n, g = g, y = y, Y = Y
    )
    rands <- make_dataset(
        ntrials = 100, ns.type = "rand", ncores = ncores, n = n, g = g, y = y, Y = Y
    )

                                        # The error stored in the node sets (opts, rands, fixeds) is based on
                                        # the "home" dynamics simulated above. The empirical test is on the
                                        # reserved test data stored in X.test. We can use obj_fn() to compute
                                        # the error on that data.
    testerror <- list(
        opt = sapply(opts, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test)),
        rand = sapply(rands, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test))
    )
    
    return(list(network = g, opts = opts, rands = rands, testerror = testerror))
}


datalist <- lapply(datafiles, make_nodesets)
                                        # 2025-05-18: Take the opportunity to clean up file names
saveRDS(datalist, "datalist-large.rds")

## check <- unlist(mclapply(datafiles, make_nodesets, mc.cores = ncores))
## saveRDS(check, "check5.rds")

## graphlist <- lapply(datafiles, make_nodesets)
## saveRDS(graphlist, "graphlist-large.rds")
