                                        # optimization
calc_obj <- function(z, y) sum((z - y)^2)/length(y)#*mean(y)) # removing normalization

obj_fn <- function(vs, y, Y, optimize_weights = FALSE, ws = NULL, ...) {
    Z <- as.matrix(Y[, vs])

    if(optimize_weights) {
        if(is.null(ws)) ws <- quadoptm(vs, y, Y)
        z <- apply(Z, 1, weighted.mean, ws)
    } else {
        z <- rowMeans(Z)
    }

    calc_obj(z, y)
}

update_vs <- function(vs, g, y, Y, optimize_weights) {
    toreplace <- sample.local(1:length(vs), 1)
    replacewith <- sample.local(as.numeric(V(g)[-which(V(g) %in% vs)]), 1)
    vs[toreplace] <- replacewith
    return(vs)
}

                                        # make nodesets
sample.local <- function(x, ...) {
    x[sample.int(length(x), ...)]
}

sample.frombins <- function(bins) {
    sample.frombin <- function(bin, out) {
        avail <- bins[[i]][!(bins[[i]] %in% out)]
        if(length(avail) == 0) {
            return(NA)
        } else {
            return(sample.local(avail, size = 1))
        }
    }
    out <- numeric(length(bins))
    for(i in seq_along(bins)) out[i] <- sample.frombin(bins[[i]], out)

    return(out)
}

select_optimized <- function(n, g, y, Y,
                             optimize_weights = FALSE, sorted = TRUE, 
                             maxit = NULL, trace = FALSE) {
    if(is.null(maxit)) {
        maxit <- switch(
            optimize_weights+1, # convert from {0, 1} to {1, 2}
            50*ncol(Y), # more iterations if not optimizing weights
            25*ncol(Y)) # fewer iterations if optimizing weights
    }

    k <- degree(g)
    vs <- sample.local(as.numeric(V(g)), n)
    result <- optim(
        par = vs, fn = obj_fn, gr = update_vs, g = g, y = y, Y = Y,
        optimize_weights = optimize_weights,
        method = "SANN", control = list(trace = trace, maxit = maxit, temp = 10)
    )

    vs <- result$par
    
    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

select_randomized <- function(n, g, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)

    vs <- sample.local(as.numeric(V(g)), n)

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

select_bydegseq <- function(n, g, comps, y, Y, optimize_weights = FALSE, sorted = TRUE) {
    k <- degree(g)
    ks <- k[comps]
    poss <- lapply(ks, function(ki) as.numeric(V(g)[which(k == ki)]))
    ## vs <- sapply(poss, function(vec) sample.local(vec, size = 1)) # allows repetition
    vs <- sample.frombins(poss)

    make_dl(vs, g, y, Y, k, optimize_weights, sorted = sorted)
}

make_dl <- function(vs, g, y, Y, k, optimize_weights = FALSE, ws = NULL, sorted = TRUE) {
    if(optimize_weights & is.null(ws)) ws <- quadoptm(vs, y, Y)
    error <- obj_fn(vs, y, Y, optimize_weights = optimize_weights, ws = ws) 
    ks <- degree(g)[vs]

    if(sorted) {
        sorting <- order(ks)
        vs <- vs[sorting]
        ks <- ks[sorting]
        if(optimize_weights) ws <- ws[sorting]
    }
    
    return(list(vs = vs, error = error, ks = ks, ws = ws))
}

make_dataset <- function(ntrials,
                         ns.type = c("opt", "rand", "fixed"),
                         ncores = 1,
                         ...) {
    ns.type <- match.arg(ns.type)

    selector <- switch(
        ns.type,
        opt = select_optimized,
        rand = select_randomized,
        fixed = select_bydegseq,
    )

    with(list(args = list(...)), {
        mclapply(seq_len(ntrials), function(x) do.call(selector, args), mc.cores = ncores)
    })
}


                                        # retrievers
get_vs <- function(dl) sapply(dl, `[[`, "vs")
get_error <- function(dl) sapply(dl, `[[`, "error")
get_ks <- function(dl) sapply(dl, `[[`, "ks")
