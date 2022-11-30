# negative expectation of log-likelihood function with
# respect to conditional distribution of 1(Ci = k)
# given data and current estimates for group 1: -Q1
negQ_G1 <- function(dat, theta, K, ...) {
    UseMethod("negQ_G1")
}

negQ_G1.zilonm <- function(dat, theta, K, num_Z, Z_names,
    B = NULL, tauG1 = NULL, calculate_tau = F, calculate_ll = F) {
    # beta0-5, alpha0k, alpha1k, xi0, w0k, gamma0,
    # gamma1, eta, delta
    dat_g1 <- dat[which(dat$Mobs > 0), ]
    M_group1 <- dat_g1$Mobs
    Y_group1 <- dat_g1$Y
    X_group1 <- dat_g1$X
    if (num_Z == 0) {
        Z_group1 <- NULL
    } else {
        Z_group1 <- dat_g1[, Z_names]
    }
    theta_trans <- trans(dat, theta, K, X_group1, num_Z,
        Z_group1)
    # group 1
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (theta_trans[["eta"]] == Inf) {
            pf0 <- rep(0, length(M_group1))
        } else {
            pf0 <- exp(-theta_trans[["eta"]]^2 * M_group1)
        }
        pf0[M_group1 > B] <- 0
        l1_ik <- -log(theta_trans[["delta"]]) + log(1 -
            pf0) - log(M_group1) - log(2 * pi) - log(theta_trans[["sig"]]) -
            (log(M_group1) - theta_trans[["mu_ik"]])^2/(2 *
                theta_trans[["sig"]]^2) - (Y_group1 - theta_trans[["beta0"]] -
            theta_trans[["beta1"]] * M_group1 - theta_trans[["beta2"]] -
            (theta_trans[["beta3"]] + theta_trans[["beta4"]]) *
                X_group1 - theta_trans[["beta_T_Z"]])^2/(2 *
            theta_trans[["delta"]]^2)

        bigpsi_ik <- (1 - theta_trans[["Del_i"]]) * theta_trans[["psi_ik"]]
        pow <- log(bigpsi_ik) + l1_ik
        if (!calculate_ll) {
            # calculate -Q1
            Q1_ik <- tauG1 * pow
            out <- -sum(Q1_ik)
        } else {
            # calculate -l1 alternative method for
            # calculating hessian: from log likelihood
            # function parameter estimation cannot use
            # it (has many saddle points so need EM),
            # but hessian can l1_i <- log(
            # bigpsi1*exp(l1_ik1) +
            # bigpsi2*exp(l1_ik2) )
            pow_max <- apply(pow, 1, max)
            l1_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l1_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci =
        # k) given data and current estiamtes
        if (K == 1) {
            out <- matrix(1, nrow = length(M_group1), ncol = 1)
        } else {
            num <- theta_trans[["psi_ik"]] * dlnorm(M_group1,
                theta_trans[["mu_ik"]], theta_trans[["sig"]])
            out <- num/rowSums(num)
        }
    }
    return(out)
}

negQ_G1.zinbm <- function(dat, theta, K, num_Z, Z_names,
    B = NULL, tauG1 = NULL, calculate_tau = F, calculate_ll = F) {
    # beta0-5, alpha0k, alpha1k, xi0, w0k, gamma0,
    # gamma1, eta, delta
    dat_g1 <- dat[which(dat$Mobs > 0), ]
    M_group1 <- dat_g1$Mobs
    Y_group1 <- dat_g1$Y
    X_group1 <- dat_g1$X
    if (num_Z == 0) {
        Z_group1 <- NULL
    } else {
        Z_group1 <- dat_g1[, Z_names]
    }
    theta_trans <- trans(dat, theta, K, X_group1, num_Z,
        Z_group1)
    log_dnb_nz <- logbinom(theta_trans[["r"]], M_group1) +
        M_group1 * log(1 - theta_trans[["p_ik"]]) - log(theta_trans[["p_ik"]]^(-theta_trans[["r"]]) -
        1)
    # group 1
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (theta_trans[["eta"]] == Inf) {
            pf0 <- rep(0, length(M_group1))
        } else {
            pf0 <- exp(-theta_trans[["eta"]]^2 * M_group1)
        }
        pf0[M_group1 > B] <- 0
        l1_ik <- -log(theta_trans[["delta"]]) - 0.5 * log(2 *
            pi) + log(1 - pf0) - (Y_group1 - theta_trans[["beta0"]] -
            theta_trans[["beta1"]] * M_group1 - theta_trans[["beta2"]] -
            (theta_trans[["beta3"]] + theta_trans[["beta4"]]) *
                X_group1 - theta_trans[["beta_T_Z"]])^2/(2 *
            theta_trans[["delta"]]^2) + log_dnb_nz

        bigpsi_ik <- (1 - theta_trans[["Del_i"]]) * theta_trans[["psi_ik"]]
        pow <- log(bigpsi_ik) + l1_ik
        if (!calculate_ll) {
            # calculate -Q1
            Q1_ik <- tauG1 * pow
            out <- -sum(Q1_ik)
        } else {
            # calculate -l1 alternative method for
            # calculating hessian: from log likelihood
            # function parameter estimation cannot use
            # it (has many saddle points so need EM),
            # but hessian can l1_i <- log(
            # bigpsi1*exp(l1_ik1) +
            # bigpsi2*exp(l1_ik2) )
            pow_max <- apply(pow, 1, max)
            l1_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l1_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci =
        # k) given data and current estimates
        if (K == 1) {
            out <- rep(1, length(M_group1))
        } else {
            num <- theta_trans[["psi_ik"]] * exp(log_dnb_nz)
            out <- num/rowSums(num)
        }
    }
    return(out)
}

negQ_G1.zipm <- function(dat, theta, K, num_Z, Z_names,
    B = NULL, tauG1 = NULL, calculate_tau = F, calculate_ll = F) {
    # beta0-5, alpha0k, alpha1k, xi0, w0k, gamma0,
    # gamma1, eta, delta
    dat_g1 <- dat[which(dat$Mobs > 0), ]
    M_group1 <- dat_g1$Mobs
    Y_group1 <- dat_g1$Y
    X_group1 <- dat_g1$X
    if (num_Z == 0) {
        Z_group1 <- NULL
    } else {
        Z_group1 <- dat_g1[, Z_names]
    }
    theta_trans <- trans(dat, theta, K, X_group1, num_Z,
        Z_group1)
    log_dpois_nz <- M_group1 * theta_trans[["loglambda_ik"]] +
        theta_trans[["log_em1"]]

    # group 1
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (theta_trans[["eta"]] == Inf) {
            pf0 <- rep(0, length(M_group1))
        } else {
            pf0 <- exp(-theta_trans[["eta"]]^2 * M_group1)
        }
        pf0[M_group1 > B] <- 0

        l1_ik <- -log(theta_trans[["delta"]]) + log(1 -
            pf0) - 0.5 * log(2 * pi) + log_dpois_nz - sapply(M_group1,
            function(t) {
                sum(log(1:t))
            }) - (Y_group1 - theta_trans[["beta0"]] - theta_trans[["beta1"]] *
            M_group1 - theta_trans[["beta2"]] - (theta_trans[["beta3"]] +
            theta_trans[["beta4"]]) * X_group1 - theta_trans[["beta_T_Z"]])^2/(2 *
            theta_trans[["delta"]]^2)

        bigpsi_ik <- (1 - theta_trans[["Del_i"]]) * theta_trans[["psi_ik"]]
        pow <- log(bigpsi_ik) + l1_ik
        if (!calculate_ll) {
            # calculate -Q1
            Q1_ik <- tauG1 * pow
            out <- -sum(Q1_ik)
        } else {
            # calculate -l1 alternative method for
            # calculating hessian: from log likelihood
            # function parameter estimation cannot use
            # it (has many saddle points so need EM),
            # but hessian can l1_i <- log(
            # bigpsi1*exp(l1_ik1) +
            # bigpsi2*exp(l1_ik2) )
            pow_max <- apply(pow, 1, max)
            l1_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l1_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci =
        # k) given data and current estimates
        if (K == 1) {
            out <- matrix(1, nrow = length(M_group1), ncol = 1)
        } else {
            # num <-
            # theta_trans[['psi_ik']]*exp(log_dpois_nz)
            # out <- num/rowSums(num)
            pow <- log(theta_trans[["psi_ik"]]) + log_dpois_nz
            out <- NULL
            for (k in seq_len(K)) {
                set <- (1:K)[-k]
                pow2 <- pow[, set, drop = F] - pow[, k]
                out <- cbind(out, 1/(1 + rowSums(exp(pow2))))
            }
        }
    }
    return(out)
}
