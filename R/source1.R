
# gammas_init <- function(dat) { if
# (sum(dat$Mobs == 0) == 0) { out <- c(-Inf,
# 0) } else { ind_Del <- (dat$Mobs == 0) * 1
# logitmodel <- summary(glm(ind_Del ~ dat$X,
# family = 'binomial')) out <-
# logitmodel$coefficients[, 'Estimate'] }
# return(out) #c(gamma0,gamma1) }

gammas_init <- function(dat, num_Z = 0, Z_names = NULL) {
    if (sum(dat$Mobs == 0) == 0) {
        out <- c(-Inf, rep(0, num_Z + 1))
    } else {
        dat$ind_Del <- (dat$Mobs == 0) * 1
        fm <- as.formula(paste0("ind_Del~", paste0(c("X",
            Z_names), collapse = "+"), collapse = ""))
        logitmodel <- summary(glm(fm, family = "binomial",
            data = dat))
        out <- logitmodel$coefficients[, "Estimate"]
    }
    return(out)
}

trans <- function(dat, theta, K, xval, ...) {
    UseMethod("trans")
}

ComputeInit <- function(dat, K, ...) {
    UseMethod("ComputeInit")
}



G1_init <- function(dat, init, K, ...) {
    UseMethod("G1_init")
}

# fix parameter not estimating
negQ2_G1 <- function(dat, init, K, num_Z, Z_names,
    B = NULL, tauG1 = NULL, calculate_tau = F,
    calculate_ll = F) {
    initials <- G1_init(dat, init, K, num_Z, Z_names)
    out <- negQ_G1(dat, theta = initials, K, num_Z,
        Z_names, B, tauG1, calculate_tau, calculate_ll)
    return(out)
}

compare_mu <- function(dat, init2, group, ...) {
    UseMethod("compare_mu")
}

bounds <- function(dat, K, group, ...) {
    UseMethod("bounds")
}

ComputeInit2 <- function(dat, K, B, ...) {
    UseMethod("ComputeInit2")
}

# negative expectation of log-likelihood
# function with respect to conditional
# distribution of 1(Ci = k) given data and
# current estimates for group 2: -Q2
negQ_G2 <- function(dat, theta, K, ...) {
    UseMethod("negQ_G2")
}

# # Q for all negQ <- function(dat, theta, K,
# B, tau = NULL, calculate_tau = F,
# calculate_ll = F) { if (sum(dat$Mobs == 0)
# > 0) { if (!calculate_tau) { if
# (!calculate_ll) { out <- negQ_G1(dat,
# theta, K, B, tauG1 = tau$tauG1,
# calculate_tau, calculate_ll) + negQ_G2(dat,
# theta, K, B, tauG2 = tau$tauG2,
# calculate_tau, calculate_ll) } else { out
# <- negQ_G1(dat, theta, K, B, tauG1 = NULL,
# calculate_tau, calculate_ll) + negQ_G2(dat,
# theta, K, B, tauG2 = NULL, calculate_tau,
# calculate_ll) } # print(theta);print(out)
# if (is.nan(out)) # { print(theta)
# print(out) } } else { out <- list(tauG1 =
# negQ_G1(dat, theta, K, B = NULL, tauG1 =
# NULL, calculate_tau, calculate_ll), tauG2 =
# negQ_G2(dat, theta, K, B, tauG2 = NULL,
# calculate_tau, calculate_ll)) } } else { if
# (!calculate_tau) { if (!calculate_ll) { out
# <- negQ2_G1(dat, theta, K, B, tauG1 =
# tau$tauG1, calculate_tau, calculate_ll) }
# else { out <- negQ2_G1(dat, theta, K, B,
# tauG1 = NULL, calculate_tau, calculate_ll)
# } # print(theta);print(out) } else { out <-
# list(tauG1 = negQ2_G1(dat, theta, K, B =
# NULL, tauG1 = NULL, calculate_tau,
# calculate_ll)) } } return(out) } # for
# calculating observed I (EM approximation) g
# <- function(dat, x, y, K, B) { tau_f <-
# negQ(dat, y, K, B, calculate_tau = T) out
# <- negQ(dat, x, K, B, tau_f) return(out) }

# Q for all
negQ <- function(dat, theta, K, num_Z, Z_names,
    B, tau = NULL, calculate_tau = F, calculate_ll = F) {
    if (sum(dat$Mobs == 0) > 0) {
        if (!calculate_tau) {
            if (!calculate_ll) {
                out <- negQ_G1(dat, theta, K, num_Z,
                  Z_names, B, tauG1 = tau$tauG1,
                  calculate_tau, calculate_ll) +
                  negQ_G2(dat, theta, K, num_Z,
                    Z_names, B, tauG2 = tau$tauG2,
                    calculate_tau, calculate_ll)
            } else {
                out <- negQ_G1(dat, theta, K, num_Z,
                  Z_names, B, tauG1 = NULL, calculate_tau,
                  calculate_ll) + negQ_G2(dat,
                  theta, K, num_Z, Z_names, B,
                  tauG2 = NULL, calculate_tau,
                  calculate_ll)
            }
            # print(theta);print(out)
        } else {
            out <- list(tauG1 = negQ_G1(dat, theta,
                K, num_Z, Z_names, B = NULL, tauG1 = NULL,
                calculate_tau, calculate_ll), tauG2 = negQ_G2(dat,
                theta, K, num_Z, Z_names, B, tauG2 = NULL,
                calculate_tau, calculate_ll))
        }
    } else {
        if (!calculate_tau) {
            if (!calculate_ll) {
                out <- negQ2_G1(dat, theta, K,
                  num_Z, Z_names, B, tauG1 = tau$tauG1,
                  calculate_tau, calculate_ll)
            } else {
                out <- negQ2_G1(dat, theta, K,
                  num_Z, Z_names, B, tauG1 = NULL,
                  calculate_tau, calculate_ll)
            }
            # print(theta);print(out)
        } else {
            out <- list(tauG1 = negQ2_G1(dat, theta,
                K, num_Z, Z_names, B = NULL, tauG1 = NULL,
                calculate_tau, calculate_ll))
        }
    }
    return(out)
}

# for calculating observed I (EM
# approximation)
g <- function(dat, x, y, K, num_Z, Z_names, B) {
    tau_f <- negQ(dat, y, K, num_Z, Z_names, B,
        calculate_tau = T)
    out <- negQ(dat, x, K, num_Z, Z_names, B, tau_f)
    return(out)
}

effects <- function(dat, theta, x1, x2, K, ...) {
    UseMethod("effects")
}

parameter_names <- function(dat, K, ...) {
    UseMethod("parameter_names")
}

G1_index <- function(dat, K, ...) {
    UseMethod("G1_index")
}

# analysis <- function(dat, K, B, limits,
# seed) { cinit <- init2 <- ComputeInit2(dat,
# K, B, limits) # cinit <- init2 <- true2 if
# (sum(dat$Mobs == 0) == 0) { index <-
# G1_index(dat, K) cinit <- init2 <-
# init2[index] tau2 <- negQ(dat, init2, K, B,
# calculate_tau = T) p <- length(init2)
# countEM <- NA } else { tau2 <- negQ(dat,
# init2, K, B, calculate_tau = T) p <-
# length(init2) init <- rep(0, p) countEM <-
# 0 bd <- bounds(dat, K, group = 2) # M-step
# while (euclidean_dist(init, init2) >
# limits) { init <- init2 tau <- tau2 # m2 <-
# #
# nlminb(init,function(x)negQ(dat,x,K,B,tau),
# # lower = c(rep(-1000,6 + 2*K),rep(0, #
# K-1),-1000,-1000,1e-6,1e-6), upper = #
# c(rep(1000,6 + 2*K),rep(1, #
# K-1),rep(1000,4)), control = #
# list(eval.max=5000,iter.max=5000)) m <-
# constrOptim(init, function(x) negQ(dat, x,
# K, B, tau), grad = NULL, ui = bd$ui, ci =
# bd$ci, outer.iterations = 500, control =
# list(maxit = 50000)) #
# negQ(dat,true,K,B,negQ(dat,true,K,B,calculate_tau=T))
# if (m$convergence != 0) {
# print(paste0('seed=', seed, ', countEM=',
# countEM, ', convergence=', m$convergence))
# print(m) init2 <- NA break } else { init2
# <- m$par } switch <- compare_mu(dat, init2,
# group = 2, dat$X, K, B) init2 <-
# switch$init2 tau2 <- switch$tau2 countEM <-
# countEM + 1 } } est <- init2 n <-
# dim(dat)[1] negll <- negQ(dat, est, K, B,
# calculate_ll = T) AIC <- 2 * negll + 2 * p
# BIC <- 2 * negll + p * log(n)
# return(list(init = cinit, est = est,
# countEM = countEM, AIC = AIC, BIC = BIC,
# tau2 = tau2)) }

analysis <- function(dat, K, num_Z, Z_names, B,
    limits, seed) {
    cinit <- init2 <- ComputeInit2(dat, K, num_Z,
        Z_names, B, limits, explicit = T)
    # cinit <- init2 <- true2
    if (sum(dat$Mobs == 0) == 0) {
        index <- G1_index(dat, K)
        cinit <- init2 <- init2[index]
        tau2 <- negQ(dat, init2, K, num_Z, Z_names,
            B, calculate_tau = T)
        p <- length(init2)
        countEM <- NA
    } else {
        tau2 <- negQ(dat, init2, K, num_Z, Z_names,
            B, calculate_tau = T)
        p <- length(init2)
        init <- rep(0, p)
        countEM <- 0
        bd <- bounds(dat, K, group = 2, num_Z)
        # M-step
        while (euclidean_dist(init, init2) > limits) {
            init <- init2
            tau <- tau2

            # m2 <-
            # nlminb(init,function(x)negQ(dat,x,K,num_Z,Z_names,B,tau),
            # lower = c(rep(-1000,6 +
            # 2*K),rep(0,
            # K-1),-1000,-1000,1e-6,1e-6),
            # upper = c(rep(1000,6 +
            # 2*K),rep(1, K-1),rep(1000,4)),
            # control =
            # list(eval.max=5000,iter.max=5000))
            m <- constrOptim(init, function(x) negQ(dat,
                x, K, num_Z, Z_names, B, tau),
                grad = NULL, ui = bd$ui, ci = bd$ci,
                outer.iterations = 500, control = list(maxit = 50000))
            # negQ(dat,true,K,num_Z,Z_names,B,negQ(dat,true,K,num_Z,Z_names,B,calculate_tau=T))
            if (m$convergence != 0) {
                print(paste0("seed=", seed, ", countEM=",
                  countEM, ", convergence=", m$convergence))
                print(m)
                init2 <- NA
                break
            } else {
                init2 <- m$par
            }
            switch <- compare_mu(dat, init2, group = 2,
                K, num_Z, Z_names, B)
            init2 <- switch$init2
            tau2 <- switch$tau2
            countEM <- countEM + 1
            # print(countEM)
        }
    }
    est <- init2
    n <- dim(dat)[1]
    negll <- negQ(dat, est, K, num_Z, Z_names,
        B, calculate_ll = T)
    AIC <- 2 * negll + 2 * p
    BIC <- 2 * negll + p * log(n)
    out <- list(init = cinit, est = est, countEM = countEM,
        AIC = AIC, BIC = BIC, tau2 = tau2)
    return(out)
}

# analysis2 <- function(dat, K, B, est, tau2,
# x1, x2) { # method 1: hessian from
# log-likelihood function, # parameter
# estimation cannot do this hess <-
# pracma::hessian(function(x) negQ(dat, x, K,
# B, calculate_ll = T), x = est) vcovar <-
# solve(hess) var <- diag(vcovar) if (sum(var
# < 0) > 0) { hess <-
# numDeriv::hessian(function(x) negQ(dat, x,
# K, B, calculate_ll = T), x = est) vcovar <-
# solve(hess) var <- diag(vcovar) if (sum(var
# < 0) > 0) { # method 2: hessian matrix from
# EM # approximation h.1 <-
# pracma::hessian(function(x) negQ(dat, x, K,
# B, tau2), x = est) g2 <- function(y) {
# numDeriv::grad(function(x) g(dat, x, y, K,
# B), x = est) } h.2 <- jacobian(g2, est)
# hess <- h.1 + h.2 vcovar <- solve(hess) var
# <- diag(vcovar) if (sum(var < 0) > 0) {
# print('All methods var < 0.') } } } se <-
# sqrt(var) # mediation effects Group1 <-
# sum(dat$Mobs == 0) == 0 MedZIM <-
# effects(dat, est, x1, x2, K, calculate_se =
# T, vcovar, Group1) return(list(se = se,
# MedZIM = MedZIM)) }

analysis2 <- function(dat, K, B, est, tau2, x1,
    x2, num_Z, Z_names, zval) {
    # method 1: hessian from log-likelihood
    # function, parameter estimation cannot
    # do this
    hess <- pracma::hessian(function(x) negQ(dat,
        x, K, num_Z, Z_names, B, calculate_ll = T),
        x = est)
    vcovar <- solve(hess)
    var <- diag(vcovar)
    if (sum(var < 0) > 0) {
        hess <- numDeriv::hessian(function(x) negQ(dat,
            x, K, num_Z, Z_names, B, calculate_ll = T),
            x = est)
        vcovar <- solve(hess)
        var <- diag(vcovar)
        if (sum(var < 0) > 0) {
            # method 2: hessian matrix from
            # EM approximation
            h.1 <- pracma::hessian(function(x) negQ(dat,
                x, K, num_Z, Z_names, B, tau2),
                x = est)
            g2 <- function(y) {
                grad(function(x) g(dat, x, y, K,
                  num_Z, Z_names, B), x = est)
            }
            h.2 <- jacobian(g2, est)
            hess <- h.1 + h.2
            vcovar <- solve(hess)
            var <- diag(vcovar)
            if (sum(var < 0) > 0) {
                print("All methods var < 0.")
            }
        }
    }
    se <- sqrt(var)
    # mediation effects
    Group1 <- sum(dat$Mobs == 0) == 0
    MedZIM <- effects(dat, est, x1, x2, K, num_Z,
        zval, calculate_se = T, vcovar, Group1)
    out <- list(se = se, MedZIM = MedZIM)
    return(out)
}

realanalysis <- function(dat, distM_sequence, K_sequence,
    selection, x1, x2, zval, num_Z, Z_names, limits,
    B, seed, ncore) {
    set.seed(seed)
    iter <- data.frame(distM = rep(distM_sequence,
        each = length(K_sequence)), K = rep(K_sequence,
        length(distM_sequence)))
    niter <- nrow(iter)

    # t0 <- Sys.time()
    cl <- parallel::makeCluster(ncore, outfile = "")
    fun <- c("analysis", "bounds", "compare_mu",
        "ComputeInit", "ComputeInit2", "gammas_init",
        "euclidean_dist", "expit", "k_to_ik", "loghicpp_all",
        "logbinom", "loghik_zinbm", "loghik_zipm",
        "negQ", "negQ_G1", "negQ_G2", "negQ2_G1",
        "trans")
    parallel::clusterExport(cl = cl, varlist = c("nodes",
        "weights", fun), envir = parent.env(environment()))
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl = cl, seed)
    models <- foreach(i = seq_len(niter), .multicombine = T,
        .packages = c("stats", "flexmix", "MASS"),
        .combine = "list", .errorhandling = "pass") %dopar%
        {
            dat2 <- dat
            class(dat2) <- c(iter$distM[i], class(dat2))
            model_i <- tryCatch({
                analysis(dat2, K = iter$K[i], num_Z,
                  Z_names, B, limits, seed)
            }, error = function(e) {
                # print(e)
                list(init = NA, est = NA, countEM = NA,
                  AIC = Inf, BIC = Inf, tau2 = NA,
                  e = e)
            })
            # model_i <- analysis(dat2, K =
            # iter$K[i], num_Z, Z_names,B,
            # limits, seed)
            print(paste0("distM = ", iter$distM[i],
                ", K = ", iter$K[i], " completed."))
            return(model_i)
        }
    parallel::stopCluster(cl)
    if (niter == 1) {
        models <- list(models)
    }
    names(models) <- apply(iter, 1, function(t) paste0(t,
        collapse = "_K"))
    print(models)

    AICs <- sapply(models, function(x) x[["AIC"]])
    BICs <- sapply(models, function(x) x[["BIC"]])
    selection_criteria <- sapply(models, function(x) x[[selection]])
    selected_model <- models[[which.min(selection_criteria)]]
    selected_model_name <- names(models)[[which.min(selection_criteria)]]
    selected_dist <- strsplit(selected_model_name,
        "_K")[[1]][1]
    selected_K <- as.numeric(strsplit(selected_model_name,
        "_K")[[1]][2])
    KBIC <- K_sequence[which.min(BICs)]
    KAIC <- K_sequence[which.min(AICs)]

    class(dat) <- c(selected_dist, class(dat))
    # t1 <- Sys.time() out <- analysis2(dat,
    # K = selected_K, B, est =
    # selected_model$est, tau2 =
    # selected_model$tau2, x1, x2, num_Z,
    # Z_names, zval)
    analysis2_out <- tryCatch({
        analysis2(dat, K = selected_K, B, est = selected_model$est,
            tau2 = selected_model$tau2, x1, x2,
            num_Z, Z_names, zval)
    }, error = function(e) {
        print(e)
        list(se = NA, MedZIM = data.frame(eff = rep(NA,
            3), eff_se = rep(NA, 3)), e = e)
    })
    # t2 <- Sys.time() print(paste0('EM = ',
    # format(difftime(t1, t0)), ', hessian =
    # ', format(difftime(t2, t1))))

    results_parameters <- results(est = selected_model$est,
        se = analysis2_out$se, init = selected_model$init)
    paranames <- parameter_names(dat, K = selected_K,
        num_Z, Z_names)
    if (sum(dat$Mobs == 0) == 0) {
        index <- G1_index(dat, K = selected_K,
            num_Z)
        paranames <- paranames[index]
    }

    rownames(results_parameters) <- paranames

    # mediation and direct effects
    out_MedZIM <- analysis2_out$MedZIM
    results_effects <- results(out_MedZIM$eff,
        out_MedZIM$eff_se)
    rownames(results_effects) <- c("NIE1", "NIE2",
        "NIE")

    out <- list(results_effects = results_effects,
        results_parameters = results_parameters,
        selected_model_name = selected_model_name,
        BIC = BICs[which.min(selection_criteria)],
        AIC = AICs[which.min(selection_criteria)],
        models = models, analysis2_out = analysis2_out)
    return(out)
}
