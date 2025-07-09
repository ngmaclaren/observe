
library(optparse)
optionlist <- list(
    make_option(
        "--dynamics", type = "character", default = "doublewell",
        help = "The dynamics on the network. Default is %default. Options: 'doublewell', 'SIS', 'genereg', 'mutualistic', and 'wilsoncowan'."
    ),
    make_option(
        "--test", action = "store_true", default = FALSE,
        help = "If flag is present, don't run the simulation for the whole data set. Instead, select an arbitrary individual's data and run on that data only."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

library(parallel)
ncores <- detectCores()-2
library(igraph)
library(sdn)
## library(optNS)
source("~/Documents/reduction/sims/local-optNS.R")

.wilsoncowan <- list(
    xinit.low = 0, xinit.high = 8,
    mu = 3, delta = 1,
    u = 0, use.matrix = TRUE, # because all fMRI networks are weighted
    D = 0.05, Ds = seq(0, 2, length.out = 100)
)
wilsoncowan <- function(t, x, params) {
    with(params, {
        coupling <- D*rowSums(A*outer(rep(1, length(x)), 1/(1 + exp(mu - delta*x))))
        dx <- -x + coupling + u
        return(list(c(dx)))
    })
}

if(interactive()) {
    args$dynamics <- "wilsoncowan"
    args$test <- TRUE
}

test <- args$test
dynamics <- args$dynamics

outfile <- paste0("datalist-fmri-", dynamics, ".rds")
print(outfile)

model <- get(dynamics)
params <- get(paste0(".", dynamics))

datadir <- "/projects/academic/naokimas/neil/brains-ns300/"
datafiles <- list.files(path = datadir, pattern = ".txt")

make_nodesets <- function(datafile) {
    df <- read.table(paste0(datadir, datafile), sep = " ")
    nr <- nrow(df)
    pos <- which(datafiles %in% datafile)

    X.train <- as.matrix(df[seq_len(nr-100), ])
    X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

    avgxi <- colMeans(X.test)
    mu <- mean(avgxi) 
    sigma <- sd(avgxi) 
    threshold <- mu + 5*sigma 
    test <- any(abs(avgxi) >= threshold) 
    if(isFALSE(test)) flag <- "use" else flag <- "dontuse"

    Cmat <- cor(X.train, method = "pearson")
    A <- Cmat
    A[which(A < 0, arr.ind = TRUE)] <- 0

    g <- graph_from_adjacency_matrix(A, "undirected", weighted = TRUE, diag = FALSE)
    N <- vcount(g)
    n <- floor(log(N))

    times <- 0:15
    control <- list(times = times, ncores = ncores) # 1
    params <- c(params, list(A = A))
    params$use.matrix <- TRUE
    Ds <- seq(0, 1, length.out = 100)

    xinit <- switch(dynamics, genereg = params$xinit.high, params$xinit.low)
    
    Y <- solve_in_range(
        Ds, "D", model, rep(xinit, N), params, control, "ode",
        method = "adams", maxsteps = 100000
    )
    y <- rowMeans(Y)
    print(datafile)
    print(dim(Y))

    opts <- make_dataset(
        ntrials = 100, ns.type = "opt", ncores = ncores, n = n, g = g, y = y, Y = Y
    )

    best <- opts[[which.min(get_error(opts))]]
    
    fixeds <- make_dataset(
        ntrials = 100, ns.type = "fixed", ncores = ncores, n = n, g = g, comps = best$vs, y = y, Y = Y
    )

    rands <- make_dataset(
        ntrials = 100, ns.type = "rand", ncores = ncores, n = n, g = g, y = y, Y = Y
    )
    
    testerror <- list(
        opt = sapply(opts, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test)),
        fixed = sapply(fixeds, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test)),
        rand = sapply(rands, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test))
    )
    
    return(
        list(
            ##network = g,
            opts = opts,
            fixeds = fixeds,
            rands = rands,
            testerror = testerror,
            flag = flag
        )
    )
}

if(test) {
    datafile <- datafiles[15]
    make_nodesets(datafile)
} else {
    datalist <- lapply(datafiles, make_nodesets)
    saveRDS(datalist, outfile)
}


    ## calc_obj <- function(z, y) sum((z - y)^2)/length(y)
    ## obj_fn <- function (vs, y, Y, optimize_weights = FALSE, ws = NULL, ...) {
    ##     Z <- as.matrix(Y[, vs])
    ##     if (optimize_weights) {
    ##         if (is.null(ws)) 
    ##             ws <- quadoptm(vs, y, Y)
    ##         z <- apply(Z, 1, weighted.mean, ws)
    ##     }
    ##     else {
    ##         z <- rowMeans(Z)
    ##     }
    ##     calc_obj(z, y)
    ## }

