#' Bayesian Nearest Neighours
#'
#' @description
#' Implementation of Nuti G's Algorithm for Bayesian Nearest Neighbours (see references) in R
#'
#' @param x A vector for categorical values, for class types of samples, ordered from the closest by distance to the target for prediction
#' @param alpha,beta The prior probablities for Beta distributions
#' @param p_gamma Proablity of streak extension, ie, not a breakpoint at the next position of the ordered samples
#' @param class1 The category in `x` for which the posterior probablities will be calculated
#' @param weights Weights for each sample; 
#'  - If `NA`, estimate weights using `get_sample_weight` from the input `x`;
#'  - If `NULL`, estimate weights using `get_sample_weight` from the input `x`;
#'  - Otherwise, values as a numeric vecotor are expected.
#' @param strata For stratified weights, provide strata groups for samples. Weights will be calculate per group if provided. Default: NULL, not stratified
#' @return A list containing changepoint position, and the posterior proablity
#'  - `pi`, `prob(x_t | k)`, A vector, proablities of the current position being class1 for different run lengths
#'  - `pk`, `prob(k | x)`, A vector, proablities of being a specific run length
#'  - `posterior`, posterior proablity of the current posistion is from class1
#' @references Nuti G. An efficient algorithm for bayesian nearest neighbours. Methodology and Computing in Applied Probability. 2019 Dec;21(4):1251-8.
#' @export
bnn <- function(x, alpha=10, beta=10, p_gamma = 0.05, class1=x[1], weights = NULL, strata = NULL) { # weight samples {{{
    # weights = 1 - unweighted
    #         = NULL - weight by sample frequency for unblanced data
    # x: classes for nearest neighours in the order from nearest, to next nearest, ....
    #           By default, include 5 * (alpha + beta) neighours;
    # alpha, beta: Beta prior proablity parameters;
    # p_gamma: proablity of streak extension, ie, not a breakpoint at the next position
    # output: a list containing changepoint position, and the posterior proablity
    #       - pi, prob(x_t | k), proablity of the current position being class1 for different run lengths
    #       - pk, prob(k | x), proablity of being a run lengths
    #       - posterior, posterior proablity of the current posistion is from class1
    x0 <- rev(x) # reorder: farest -> nearest
    classes <- c(class1, setdiff(x0, class1))
    x0 <- match(x0, classes)
    nc <- length(classes)
    if(is.null(weights)) {#{{{
        weights <-  get_sample_weight(x0, strata = strata)
    } else if(length(weights) == 1 && is.na(weights)) {
        weights <-  get_sample_weight(x0, strata = NA)
    } else if(length(weights) < length(x0)) {
        weights <- rep(weights, length = x0)
    }#}}}

    x <- 2 - x0 %in% 1 # class1 -> 1; class2 -> 2
    x <- c(x, 1) # add 1 to the last so as to calculate proablity of observing class1 at after finishing updating
    n <- length(x)
    x0 <- c(x0, 1)

    eta <- eta0 <- rbind(c(alpha, beta)) # priors
    eta.sum <- eta.sum0 <- sum(eta0)
    pk.prev <- pkx.prev <-  1 # pk = prob(k | x); pkx = prob(k, x)

    for(t in 1:n){
        # observe the next x_t
        #   pi(i): p(x_t | k_(t-1) = i, eta_i)
        #   pk(i): p(k_t = k_(t-1) + 1, x | k_(t-1) = i)

        # Compute predictive prob:
        #   pi_i = p(x_t | k_(t-1), eta_i)
        pi <- eta[,x[t]] / eta.sum

        if(t == n) {
            posterior = sum(pk.prev * pi, na.rm=TRUE)
            return(list(pi = pi, pk = pk.prev, posterior = posterior, bayes_factor = posterior / (1 - posterior), class1 = class1))
        }

        # Compute change-point prob:
        #   p(kt =0, x_0:t) = sum_{} 
        # Compute growth prob:
        #   p(k_t = k_(t-1) + 1, x_(0:t)) = p(k_(t-1), x_(0:t-1)) * pi_t * p_gamma (ERROR?)
        pkx <- c(sum(pkx.prev)  * p_gamma * pi[1], pkx.prev * pi * (1 - p_gamma))

        # Compute evidence
        #   p(x_0:t) = sum_kt { p(kt, x_0:t) }
        px <- sum(pkx)
        # Compute prob of k:
        #   p(k_t = i | x_0:t) = p(k_t = i, x_0:t) / p(x_0:t)
        pk <- pkx / px

        # Update distributions:
        #   eta_i <<- x_t
        incr <- weights[x0[t]]
        eta[, x[t]] <- eta[, x[t]] + incr
        eta <- rbind(eta0, eta)
        eta.sum <- c(eta.sum0, eta.sum + incr)
        pk.prev <- pk
        pkx.prev <- pkx
    }

} #}}}

#' Bayesian Nearest Neighours for a data matrix
#' @param data a matrix to calculate sample distance from, with columns corresponding to samples. Example: principle components for samples
#' @param sample_info a data.frame for sample information, columns corresponding to those of `dat`
#' @param query the name(s)/ID(s) of the samples for which BNN posteriors will be calculated
#' @param target the name(s)/ID(s) of the samples with which BNN posteriors will be calculated
#' @param k the maximum number of nearest neighours to be considered. Smaller `k` increases speed at the possible cost of missing info when remote samples are also predictive
#' @param alpha,beta The prior probablities for Beta distributions
#' @param p_gamma Proablity of streak extension, ie, not a breakpoint at the next position of the ordered samples
#' @param bnn_name either `NULL` or a string to name the columns of returned data.frame
#' @param strata For stratified weights, provide strata groups for samples. Weights will be calculate per group if provided. Default: NULL, not stratified
#'  - when `strata` is column name of the `sample_info` data.frame, the corresponding column will be used as `strata`
#' @param bnn_name either `NULL` or a string to name the columns of returned data.frame
#' @param sample_col the column name for sample names/IDs
#' @param status_col the column name for sample status, which is used in BNN prediction
#' @param class1 The category for which the posterior probablities will be calculated; See the `class1` argument of function `bnn`
#' @return Bayesian Nearest Neighours results
#'  - if `bnn_name` is `NULL`, returns a vector of Bayesian posterior probablities for each sample
#'  - if `bnn_name` is provided as a string, resturns a data.frame, with columns corresponding to the posterior probablities, and the most likely change points (`chngpnt`), with columns names as `paste0(c("bnn", "chngpnt), bnn_name)`
#' @import FNN
#' @export
calc_bnn <- function(dat, sample_info, query = NULL, target = NULL, k=200, alpha = 10, beta = 10, p_gamma = 0.05, bnn_name = NULL, strata = NULL, sample_col = "Name", status_col = "gender", class1 = "MALE"){ #{{{
   library(FNN)
    samples <- sample_info[[sample_col]]
    if(is.null(query)) query <- samples
    if(is.null(target)) target <- samples
    i_query <- match(query, samples)
    i_target <- match(target, samples)
    n_query <- nrow(query)
    n_target <- length(target)
    k <- min(k, n_target)
    x <- sample_info[[status_col]][i_target]

    dat.query <- t(dat[, i_query, drop = FALSE])
    dat.target <- t(dat[, i_target, drop = FALSE])
    nn.index <- get.knnx(dat.target, dat.query, k = k)$nn.index

    if(length(strata) == 1 && strata %in% colnames(sample_info)) strata <- sample_info[[strata]][i_target]
    weights <-  get_sample_weight(x, strata = strata)
    nn.gender <- matrix(x[nn.index], ncol=ncol(nn.index))
    bnn.ret <- lapply(1:nrow(nn.gender), function(x) bnn(nn.gender[x,] , alpha = alpha, beta = beta, p_gamma = p_gamma, class1 = class1, weights = weights[nn.index[x,]]))
    bnn_pp <- sapply(bnn.ret, function(x) x$posterior)
    # bnn.ret <- apply(nn.gender, 1, function(x) bnn(x , alpha = alpha, beta = beta, p_gamma = p_gamma, class1 = "MALE")$posterior)
    if(is.null(bnn_name)) {
        return(bnn_pp)
    } else {
        bnn_chngpnt <- sapply(bnn.ret, function(x) which.max(x$pk[-1]) - 1) # change point
        res <- data.frame(bnn_pp, bnn_chngpnt)
        names(res) <- paste0(c("bnn", "chngpnt"), bnn_name)
        return(res)
    }
} #}}}

#' Get Sample Weights
#' @param x A vector of categorical values, for class types of samples ordered by from the most to the least closest to the target for prediction
#' @param strata For stratified weights, provide strata groups for samples. Weights will be calculate per group if provided. Default: NULL, not stratified
#'  - if `NA`, use equal weigths for all non-NA samples
#'  - if `NULL`, calculate weights in all samples without stratification
#'  - otherwise, calculate weights seprately in each non-NA stratum; The total weights in a stratum is the total number of samples, the total weights of every unique category of are equal to each other.
#' @return A numeric vector for weights for each item in the input `x`
#' @export
get_sample_weight <- function(x, strata = NULL) {#{{{
    unique_x <- sort(setdiff(x, NA))
    res <- rep(1, length(x))
    res[is.na(x)] <-  NA
    if(is.null(strata)) strata <- res
    if(length(strata) == 1 && is.na(strata)) return(res)

    strata_groups <- setdiff(strata, NA)
    for(g in strata_groups) {#{{{
        i <- strata %in% g
        g_count <- sum(i)
        nc <- length(setdiff(x[i], NA))
        for(v in unique_x) {#{{{
            j <- i & x %in% v
            if(any(j)) res[j] <- g_count / nc / sum(j)
        }#}}}
    }#}}}
    res[!is.na(x) & is.na(strata)] <- 1

    return(res)
}#}}}

## internal functions
test_bnn__ <- function() {#{{{
    set.seed(0)
    res <- list(x = list(),
                prob1 = list(),
                prob2 = list(),
                n1 = list(),
                n2 = list(),
                bnn = list()
    )
    simulate <- function(i, ...) {#{{{
        res <- res
        input <- list(...)
        assign <- function(term) {#{{{
            if(is.null(input[[term]]) && i > 1) {
                res[[term]][[i]] <<- res[[term]][[i - 1]]
            } else {
                res[[term]][[i]] <<- input[[term]]
            }
        }#}}}
        for(term in c("n1", "n2", "prob1", "prob2", "x")) assign(term)
        if(any(names(input) %in% c("n1", "n2"))) res$x[[i]] <- c(sample(2, res$n1[[i]], TRUE, res$prob1[[i]]),
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
        if(any(names(input) %in% c("prob1"))) res$x[[i]][1:res$n1[[i]]] <- sample(2, res$n1[[i]], TRUE, res$prob1[[i]])
        if(any(names(input) %in% c("prob2"))) res$x[[i]][res$n1[[i]] + 1:res$n2[[i]]] <- sample(2, res$n2[[i]], TRUE, res$prob2[[i]])
        return(res)
    }#}}}

    # null case with equal chance of two groups
    i <- 1
    res <- simulate(i, prob1 = c(0.5, 0.5), prob2 = c(0.1, 0.9), n1 = 10, n2 = 20)
    res$bnn[[i]] <- bnn(res$x[[i]], 1, 1)

    # more chance of being 2
    i <- 2
    res <- simulate(i, prob1 = c(0.1, 0.9), prob2 = c(0.5, 0.5))
    res$bnn[[i]] <- bnn(res$x[[i]], 1, 1)

    # more chance of being 2, with biased background
    i <- 3
    res <- simulate(i, prob2 = c(0.8, 0.2))
    res$bnn[[i]] <- bnn(res$x[[i]], 1, 1)


    # more chance of being 2
    i <- 3
    res$prob1[[i]] <- res$prob1[[i-1]]
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn(res$x[[i]])

    # more chance of being 2
    i <- 2
    res$prob1[[i]] <- c(0.1, 0.2)
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn(res$x[[i]])

    # more chance of being 2
    i <- 2
    res$prob1[[i]] <- c(0.1, 0.2)
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn(res$x[[i]])

    # more chance of being 2
    i <- 2
    res$prob1[[i]] <- c(0.1, 0.2)
    res$prob2[[i]] <- c(0.5, 0.5)
    res$n1[[i]] <- 10; res$n2[[i]] <- 20
    res$x[[i]] <- c(sample(2, res$n1[[i]], replace = TRUE, res$prob1[[i]]), 
                    sample(2, res$n2[[i]], replace = TRUE, res$prob2[[i]]))
    res$bnn[[i]] <- bnn(res$x[[i]])

    return(res)
}#}}}
