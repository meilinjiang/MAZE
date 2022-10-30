
gammas_init <- function(dat) {
    if (sum(dat$Mobs == 0) == 0) {
        out <- c(-Inf, 0)
    } else {
        ind_Del <- (dat$Mobs == 0) * 1
        logitmodel <- summary(glm(ind_Del ~ dat$X, family = "binomial"))
        out <- logitmodel$coefficients[, "Estimate"]
    }
    return(out)  #c(gamma0,gamma1)
}

trans <- function(theta, x, K, distM) {
    beta0 <- theta[1]
    beta1 <- theta[2]
    beta2 <- theta[3]
    beta3 <- theta[4]
    beta4 <- theta[5]
    if (distM == "zilonm") {
        alpha_k <- theta[5 + 1:(2 * K)]
        xi0 <- theta[6 + 2 * K]
        if (K == 1) {
            psi_k <- 1
        } else {
            psi_k <- c(theta[6 + 2 * K + 1:(K - 1)], 1 - sum(theta[6 +
                2 * K + 1:(K - 1)]))
        }
        gamma0 <- theta[6 + 3 * K]
        gamma1 <- theta[7 + 3 * K]
        eta <- theta[8 + 3 * K]
        delta <- theta[9 + 3 * K]
        n0 <- length(x)
        alpha_ik <- k_to_ik(alpha_k, n0)
        mu_ik <- alpha_ik[, 1:K, drop = F] + alpha_ik[, K + (1:K)] *
            x
        sig <- exp(xi0)
        psi_ik <- k_to_ik(psi_k, n0)
        Del_i <- expit(gamma0 + gamma1 * x)
        theta_trans <- list(beta0 = beta0, beta1 = beta1, beta2 = beta2,
            beta3 = beta3, beta4 = beta4, psi_k = psi_k, eta = eta,
            delta = delta, mu_ik = mu_ik, sig = sig, psi_ik = psi_ik,
            Del_i = Del_i)
    } else if (distM == "zinbm") {
        r <- theta[6]
        alpha_k <- theta[6 + 1:(2 * K)]
        if (K == 1) {
            psi_k <- 1
        } else {
            psi_k <- c(theta[6 + 2 * K + 1:(K - 1)], 1 - sum(theta[6 +
                2 * K + 1:(K - 1)]))
        }
        gamma0 <- theta[6 + 3 * K]
        gamma1 <- theta[7 + 3 * K]
        eta <- theta[8 + 3 * K]
        delta <- theta[9 + 3 * K]
        n0 <- length(x)
        alpha_ik <- k_to_ik(alpha_k, n0)
        logmu_ik <- alpha_ik[, 1:K, drop = F] + alpha_ik[, K + (1:K)] *
            x
        mu_ik <- exp(logmu_ik)
        psi_ik <- k_to_ik(psi_k, n0)
        Delstar_i <- expit(gamma0 + gamma1 * x)
        p_ik <- r/(r + mu_ik)
        Del_i <- Delstar_i + (1 - Delstar_i) * rowSums(psi_ik * p_ik^r)
        theta_trans <- list(beta0 = beta0, beta1 = beta1, beta2 = beta2,
            beta3 = beta3, beta4 = beta4, r = r, alpha_k = alpha_k,
            psi_k = psi_k, gamma0 = gamma0, gamma1 = gamma1, eta = eta,
            delta = delta, alpha_ik = alpha_ik, logmu_ik = logmu_ik,
            mu_ik = mu_ik, psi_ik = psi_ik, Delstar_i = Delstar_i, p_ik = p_ik,
            Del_i = Del_i)
    } else if (distM == "zipm") {
        alpha_k <- theta[5 + 1:(2 * K)]
        if (K == 1) {
            psi_k <- 1
        } else {
            psi_k <- c(theta[5 + 2 * K + 1:(K - 1)], 1 - sum(theta[5 +
                2 * K + 1:(K - 1)]))
        }
        gamma0 <- theta[5 + 3 * K]
        gamma1 <- theta[6 + 3 * K]
        eta <- theta[7 + 3 * K]
        delta <- theta[8 + 3 * K]
        n0 <- length(x)
        alpha_ik <- k_to_ik(alpha_k, n0)
        loglambda_ik <- alpha_ik[, 1:K, drop = F] + alpha_ik[, K + (1:K)] *
            x
        lambda_ik <- exp(loglambda_ik)
        psi_ik <- k_to_ik(psi_k, n0)
        Delstar_i <- expit(gamma0 + gamma1 * x)
        Del_i <- Delstar_i + (1 - Delstar_i) * rowSums(psi_ik * exp(-lambda_ik))
        # -log[exp(lambda) - 1]: when lambda too large or small
        log_em1 <- -loglambda_ik
        log_em1[lambda_ik > 1e-10] <- (-lambda_ik - log(1 - exp(-lambda_ik)))[lambda_ik >
            1e-10]
        theta_trans <- list(beta0 = beta0, beta1 = beta1, beta2 = beta2,
            beta3 = beta3, beta4 = beta4, psi_k = psi_k, eta = eta,
            delta = delta, alpha_ik = alpha_ik, loglambda_ik = loglambda_ik,
            lambda_ik = lambda_ik, psi_ik = psi_ik, Delstar_i = Delstar_i,
            Del_i = Del_i, log_em1 = log_em1)
    }
    return(theta_trans)
}

ComputeInit <- function(dat, K, distM) {
    # group 1, flexmix
    M_group1 <- dat$Mobs[which(dat$Mobs > 0)]
    Y_group1 <- dat$Y[which(dat$Mobs > 0)]
    X_group1 <- dat$X[which(dat$Mobs > 0)]

    fit <- summary(lm(Y_group1 ~ M_group1 + X_group1))
    beta02 <- fit$coefficients["(Intercept)", "Estimate"]
    beta1 <- fit$coefficients[2, "Estimate"]
    beta34 <- fit$coefficients[3, "Estimate"]
    delta <- fit$sigma

    if (distM == "zilonm") {
        cl <- flexmix(log(M_group1) ~ X_group1, k = K)
        cl_size <- table(clusters(cl))
        while (length(cl_size) != K) {
            cl <- flexmix(log(M_group1) ~ X_group1, k = K)
            cl_size <- table(clusters(cl))
        }
        # method 3: distribution 1 as the one with smaller cluster
        # mean since x=age, intercepts are not accurate enough
        mu <- tapply(log(M_group1), clusters(cl), mean)
        # mu <-
        # parameters(cl)[1,]+parameters(cl)[2,]*mean(X_group1)
        ord <- order(mu)
        mix <- parameters(cl)[, ord, drop = F]
        # alpha_0k, alpha_1k, xi_0, psik
        mix_init <- c(mix[1, ], mix[2, ], log(mean(mix[3, ])))
        psi_k <- cl_size[ord]/length(M_group1)
        mix_init <- c(mix_init, psi_k[-K])

        init <- c(beta02, beta1, beta34, mix_init, delta)
    } else if (distM == "zinbm") {
        # cl <- flexmix(M_group1 ~ X_group1, k=K, model =
        # FLXMRnegbin())
        cl <- kmeans(M_group1, centers = K)
        para <- NULL
        # for (j in seq_len(K)) { fit2 <-
        # glm.nb(M_group1[cl$cluster==j] ~
        # X_group1[cl$cluster==j], control =
        # glm.control(maxit=1e15)) para <- cbind(para,
        # c(coef(fit2), fit2$theta )) }
        for (j in seq_len(K)) {
            m_j <- M_group1[cl$cluster == j]
            # r.init <- mean(m_j)^2/(var(m_j)-mean(m_j)) #
            # numerator larger, denominator smaller than true
            r.init <- mean(m_j)/2
            fit2 <- lm(log(M_group1[cl$cluster == j]) ~ X_group1[cl$cluster ==
                j])
            para <- cbind(para, c(coef(fit2), r.init))
        }

        # method 3: distribution 1 as the one with smaller cluster
        # mean since x=age, intercepts are not accurate enough
        mu <- cl$centers
        # mu <-
        # parameters(cl)[1,]+parameters(cl)[2,]*mean(X_group1)
        ord <- order(mu)
        mix <- para[, ord, drop = F]
        # r, alpha_0k, alpha_1k, xi_0, psik
        mix_init <- c(mean(mix[3, ]), mix[1, ], mix[2, ])
        psi_k <- cl$size[ord]/length(M_group1)
        mix_init <- c(mix_init, psi_k[-K])

        init <- c(beta02, beta1, beta34, mix_init, delta)
    } else if (distM == "zipm") {
        cl <- flexmix(M_group1 ~ X_group1, k = K, model = FLXMRglm(family = "poisson"))
        cl_size <- table(clusters(cl))
        while (length(cl_size) != K) {
            cl <- flexmix(M_group1 ~ X_group1, k = K, model = FLXMRglm(family = "poisson"))
            cl_size <- table(clusters(cl))
        }
        # method 3: distribution 1 as the one with smaller cluster
        # mean since x=age, intercepts are not accurate enough
        mu <- tapply(M_group1, clusters(cl), mean)
        # mu <-
        # parameters(cl)[1,]+parameters(cl)[2,]*mean(X_group1)
        ord <- order(mu)
        mix <- parameters(cl)[, ord, drop = F]
        # alpha_0k, alpha_1k, xi_0, psik
        mix_init <- c(mix[1, ], mix[2, ])
        psi_k <- cl_size[ord]/length(M_group1)
        mix_init <- c(mix_init, psi_k[-K])

        init <- c(beta02, beta1, beta34, mix_init, delta)
    }
    return(init)
}

logbinom <- function(r, M_i) {
    # sapply(M_i, function(m){lgamma(r+m) - lgamma(r) -
    # lgamma(m+1)})
    lgamma(r + M_i) - lgamma(r) - lgamma(M_i + 1)
}

# negative expectation of log-likelihood function with repect to
# conditional distribution of 1(Ci = k) given data and current
# estiamtes for gorup 1: -Q1
negQ_G1 <- function(theta, dat, distM, K, B = NULL, tauG1 = NULL, calculate_tau = F,
    calculate_ll = F) {
    # beta0-5, alpha0k, alpha1k, xi0, w0k, gamma0, gamma1, eta,
    # delta
    M_group1 <- dat$Mobs[which(dat$Mobs > 0)]
    Y_group1 <- dat$Y[which(dat$Mobs > 0)]
    X_group1 <- dat$X[which(dat$Mobs > 0)]
    theta_trans <- trans(theta, X_group1, K, distM)

    # group 1
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (theta_trans[["eta"]] == Inf) {
            pf0 <- rep(0, length(M_group1))
        } else {
            pf0 <- exp(-theta_trans[["eta"]]^2 * M_group1)
        }
        pf0[M_group1 > B] <- 0
        if (distM == "zilonm") {
            l1_ik <- -log(theta_trans[["delta"]]) + log(1 - pf0) - log(M_group1) -
                log(2 * pi) - log(theta_trans[["sig"]]) - (log(M_group1) -
                theta_trans[["mu_ik"]])^2/(2 * theta_trans[["sig"]]^2) -
                (Y_group1 - theta_trans[["beta0"]] - theta_trans[["beta1"]] *
                  M_group1 - theta_trans[["beta2"]] - (theta_trans[["beta3"]] +
                  theta_trans[["beta4"]]) * X_group1)^2/(2 * theta_trans[["delta"]]^2)
        } else if (distM == "zinbm") {
            log_dnb_nz <- logbinom(theta_trans[["r"]], M_group1) + M_group1 *
                log(1 - theta_trans[["p_ik"]]) - log(theta_trans[["p_ik"]]^(-theta_trans[["r"]]) -
                1)
            l1_ik <- -log(theta_trans[["delta"]]) - 0.5 * log(2 * pi) +
                log(1 - pf0) - (Y_group1 - theta_trans[["beta0"]] -
                theta_trans[["beta1"]] * M_group1 - theta_trans[["beta2"]] -
                (theta_trans[["beta3"]] + theta_trans[["beta4"]]) *
                  X_group1)^2/(2 * theta_trans[["delta"]]^2) + log_dnb_nz
        } else if (distM == "zipm") {
            log_dpois_nz <- M_group1 * theta_trans[["loglambda_ik"]] +
                theta_trans[["log_em1"]]
            l1_ik <- -log(theta_trans[["delta"]]) + log(1 - pf0) - 0.5 *
                log(2 * pi) + log_dpois_nz - sapply(M_group1, function(t) {
                sum(log(1:t))
            }) - (Y_group1 - theta_trans[["beta0"]] - theta_trans[["beta1"]] *
                M_group1 - theta_trans[["beta2"]] - (theta_trans[["beta3"]] +
                theta_trans[["beta4"]]) * X_group1)^2/(2 * theta_trans[["delta"]]^2)
        }

        bigpsi_ik <- (1 - theta_trans[["Del_i"]]) * theta_trans[["psi_ik"]]
        pow <- log(bigpsi_ik) + l1_ik
        if (!calculate_ll) {
            # calculate -Q1
            Q1_ik <- tauG1 * pow
            out <- -sum(Q1_ik)
        } else {
            # calculate -l1 alternative method for calculating
            # hessian: from log likelihood function parameter
            # estimation cannot use it (has many saddle points so
            # need EM), but hessian can l1_i <- log(
            # bigpsi1*exp(l1_ik1) + bigpsi2*exp(l1_ik2) )
            pow_max <- apply(pow, 1, max)
            l1_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l1_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci = k) given
        # data and current estimates
        if (K == 1) {
            out <- matrix(1, nrow = length(M_group1), ncol = 1)
        } else {
            if (distM == "zilonm") {
                num <- theta_trans[["psi_ik"]] * dlnorm(M_group1, theta_trans[["mu_ik"]],
                  theta_trans[["sig"]])
                out <- num/rowSums(num)
            } else if (distM == "zinbm") {
                log_dnb_nz <- logbinom(theta_trans[["r"]], M_group1) +
                  M_group1 * log(1 - theta_trans[["p_ik"]]) - log(theta_trans[["p_ik"]]^(-theta_trans[["r"]]) -
                  1)
                num <- theta_trans[["psi_ik"]] * exp(log_dnb_nz)
                out <- num/rowSums(num)
            } else if (distM == "zipm") {
                log_dpois_nz <- M_group1 * theta_trans[["loglambda_ik"]] +
                  theta_trans[["log_em1"]]
                # num <- theta_trans[['psi_ik']]*exp(log_dpois_nz)
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
    }
    return(out)
}

# fix parameter not estimating
negQ2_G1 <- function(init, dat, distM, K, B = NULL, tauG1 = NULL, calculate_tau = F,
    calculate_ll = F) {
    gammas <- gammas_init(dat)
    if (distM == "zilonm") {
        initials <- c(init[1:2], 0, init[3], 0, init[3 + 1:(3 * K)],
            gammas, 0.01, init[4 + 3 * K])
    } else if (distM == "zinbm") {
        initials <- c(init[1:2], 0, init[3], 0, init[3 + 1:(3 * K)],
            gammas, 0.01, init[4 + 3 * K])
    } else if (distM == "zipm") {
        initials <- c(init[1:2], 0, init[3], 0, init[3 + 1:(3 * K -
            1)], gammas, 0.01, init[3 + 3 * K])
    }
    out <- negQ_G1(theta = initials, dat, distM, K, B, tauG1, calculate_tau,
        calculate_ll)
    return(out)
}

compare_mu <- function(init2, dat, distM, group, xi, K, B = NULL) {
    x <- mean(xi)
    s <- ifelse(group == 1, 3, 5)
    if (distM == "zilonm") {
        mu_k <- init2[s + (1:K)] + init2[s + K + (1:K)] * x
        ord <- order(mu_k)
        if (sum(ord != 1:K) > 0) {
            init2[s + 1:(2 * K)] <- init2[s + 1:(2 * K)][c(ord, ord +
                K)]
            init2[s + 2 * K + 2:K] <- c(init2[s + 2 * K + 2:K], 1 -
                sum(init2[s + 2 * K + 2:K]))[ord][-K]
        }
    } else if (distM == "zinbm") {
        mu_k <- exp(init2[s + 1 + (1:K)] + init2[s + 1 + K + (1:K)] *
            x)
        ord <- order(mu_k)
        if (sum(ord != 1:K) > 0) {
            # print(init2);print(mu_k)
            init2[s + 1 + 1:(2 * K)] <- init2[s + 1 + 1:(2 * K)][c(ord,
                ord + K)]
            init2[s + 1 + 2 * K + 1:(K - 1)] <- c(init2[s + 1 + 2 *
                K + 1:(K - 1)], 1 - sum(init2[s + 1 + 2 * K + 1:(K -
                1)]))[ord][-K]
            # print(init2)
        }
    } else if (distM == "zipm") {
        mu_k <- exp(init2[s + (1:K)] + init2[s + K + (1:K)] * x)
        ord <- order(mu_k)
        if (sum(ord != 1:K) > 0) {
            # print(init2);print(mu_k)
            init2[s + 1:(2 * K)] <- init2[s + 1:(2 * K)][c(ord, ord +
                K)]
            init2[s + 2 * K + 1:(K - 1)] <- c(init2[s + 2 * K + 1:(K -
                1)], 1 - sum(init2[s + 2 * K + 1:(K - 1)]))[ord][-K]
            # print(init2)
        }
    }

    if (group == 1) {
        tau2 <- negQ2_G1(init2, dat, distM, K, calculate_tau = T)
    } else {
        tau2 <- negQ(init2, dat, distM, K, B, calculate_tau = T)
    }
    return(list(init2 = init2, tau2 = tau2))
}

bounds <- function(K, group, distM) {
    s <- ifelse(group == 1, 3, 5)
    e <- ifelse(group == 1, 1, 4)
    ul <- 1000
    if (distM == "zilonm") {
        if (K == 1) {
            ui <- diag(rep(1, s + 3 + e))
            ci <- c(rep(-1000, s + 3), c(1e-06, 1e-06, -1000, -1000)[e:1])
        } else {
            ui1 <- diag(rep(1, s + 3 * K + e))
            ui2 <- diag(c(rep(0, s + 2 * K + 1), rep(-1, K - 1), rep(0,
                e)))
            ui <- rbind(ui1, ui2[s + 2 * K + 1 + 1:(K - 1), ], c(rep(0,
                s + 2 * K + 1), rep(-1, K - 1), rep(0, e)))
            ci <- c(rep(-1000, s + 2 * K + 1), rep(1e-06, K - 1), c(1e-06,
                1e-06, -1000, -1000)[e:1], rep(-(1 - 1e-06), K - 1),
                -(1 - 1e-06))
        }
    } else if (distM == "zinbm") {
        if (K == 1) {
            ui <- diag(rep(1, s + 3 + e))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, 2), c(1e-06, 1e-06,
                -1000, -1000)[e:1])
        } else {
            ui1 <- diag(rep(1, s + 1 + (3 * K - 1) + e))
            ui <- rbind(ui1, -ui1[s + 1 + 2 * K + 1:(K - 1), ], c(rep(0,
                s + 1 + 2 * K), rep(-1, K - 1), rep(0, e)), c(rep(0,
                s), -1, rep(0, 3 * K - 1 + e)))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, 2 * K), rep(1e-06,
                K - 1), c(1e-06, 1e-06, -1000, -1000)[e:1], rep(-(1 -
                1e-06), K - 1), -(1 - 1e-06), -1000)
        }
    } else if (distM == "zipm") {
        if (K == 1) {
            ui <- diag(rep(1, s + 2 + e))
            ci <- c(rep(-1000, s + 2), c(1e-06, 1e-06, -1000, -1000)[e:1])
        } else {
            ui1 <- diag(rep(1, s + (3 * K - 1) + e))
            ui <- rbind(ui1, -ui1[s + 2 * K + 1:(K - 1), ], c(rep(0,
                s + 2 * K), rep(-1, K - 1), rep(0, e)))
            ci <- c(rep(-1000, s + 2 * K), rep(1e-06, K - 1), c(1e-06,
                1e-06, -1000, -1000)[e:1], rep(-(1 - 1e-06), K - 1),
                -(1 - 1e-06))
        }
    }
    return(list(ui = ui, ci = ci))
}

ComputeInit2 <- function(dat, distM, K, B, limits, explicit = T) {
    init2 <- ComputeInit(dat, K, distM)
    tau2 <- negQ2_G1(init2, dat, distM, K, calculate_tau = T)
    init <- rep(0, length(init2))
    countEM <- 0
    X_group1 <- dat$X[which(dat$Mobs > 0)]
    bd <- bounds(K, group = 1, distM)
    # M-step
    if (distM == "zilonm") {
        if (explicit) {
            # Method 1: explicit solution for group 1
            n1 <- length(X_group1)
            M_group1 <- dat$Mobs[which(dat$Mobs > 0)]
            Y_group1 <- dat$Y[which(dat$Mobs > 0)]
            while (euclidean_dist(init, init2) > limits) {
                init <- init2
                tau <- tau2

                psi_k <- colSums(tau)/n1
                alpha1_ik0 <- k_to_ik(init2[3 + K + 1:K], n1)
                alpha0_k <- colSums(tau * (log(M_group1) - alpha1_ik0 *
                  X_group1))/colSums(tau)
                alpha0_ik <- k_to_ik(alpha0_k, n1)
                alpha1_k <- colSums(tau * X_group1 * (log(M_group1) -
                  alpha0_ik))/colSums(tau * X_group1^2)
                alpha1_ik <- k_to_ik(alpha1_k, n1)
                sig2 <- sum(tau * (log(M_group1) - alpha0_ik - alpha1_ik *
                  X_group1)^2)/n1
                xi0 <- log(sqrt(sig2))
                beta02 <- sum(Y_group1 - init[2] * M_group1 - init[3] *
                  X_group1)/n1
                beta1 <- sum((Y_group1 - beta02 - init[3] * X_group1) *
                  M_group1)/sum(M_group1^2)
                beta34 <- sum((Y_group1 - beta02 - beta1 * M_group1) *
                  X_group1)/sum(X_group1^2)
                delta <- sqrt(sum((Y_group1 - beta02 - beta1 * M_group1 -
                  beta34 * X_group1)^2)/length(M_group1))
                init2 <- c(beta02, beta1, beta34, alpha0_k, alpha1_k,
                  xi0, psi_k[-K], delta)
                switch <- compare_mu(init2, dat, distM, group = 1, X_group1,
                  K)
                init2 <- switch$init2
                tau2 <- switch$tau2
                countEM <- countEM + 1
                # print(init2)
                if (K == 1) {
                  break
                }
            }
        } else {
            # Method 2: R optimization function
            while (euclidean_dist(init, init2) > limits) {
                init <- init2
                tau <- tau2

                # m <-
                # nlminb(init,function(x)negQ2_G1(x,dat,distM,K,B,tau),lower
                # = c(rep(-1000,9),1e-6))
                m <- constrOptim(init, function(x) negQ2_G1(x, dat,
                  distM, K, B, tau), grad = NULL, ui = bd$ui, ci = bd$ci,
                  outer.iterations = 500, control = list(maxit = 5000))
                if (m$convergence != 0) {
                  print(paste0("ComputeInit2=", m$convergence))
                  init2 <- NA
                  break
                } else {
                  init2 <- m$par
                }
                switch <- compare_mu(init2, dat, distM, group = 1, X_group1,
                  K)
                init2 <- switch$init2
                tau2 <- switch$tau2
                countEM <- countEM + 1
                if (K == 1) {
                  break
                }
            }
        }
        # beta0,beta1,beta3,beta5, alpha_0k, alpha_1k, xi_0, w_0k,
        # delta
        gammas <- gammas_init(dat)
        initials <- c(init2[1:2], 0, init2[3], 0, init2[3 + 1:(3 * K)],
            gammas, 0.01, init2[4 + 3 * K])
    } else if (distM == "zinbm") {
        while (euclidean_dist(init, init2) > limits) {
            init <- init2
            tau <- tau2

            # m <-
            # nlminb(init,function(x)negQ2_G1(x,dat,distM,K,B,tau),lower
            # = c(rep(-1000,9),1e-6))
            m <- constrOptim(init, function(x) negQ2_G1(x, dat, distM,
                K, B, tau), grad = NULL, ui = bd$ui, ci = bd$ci, outer.iterations = 500,
                control = list(maxit = 50000))
            if (m$convergence != 0) {
                print(paste0("ComputeInit2=", m$convergence))
                print(m)
                init2 <- NA
                break
            } else {
                init2 <- m$par
            }
            switch <- compare_mu(init2, dat, distM, group = 1, X_group1,
                K)
            init2 <- switch$init2
            tau2 <- switch$tau2
            countEM <- countEM + 1
            if (K == 1) {
                break
            }
        }
        # beta0,beta1,beta3, alpha_0k, alpha_1k, psi_k, delta
        gammas <- gammas_init(dat)
        # initials <- c(init2[1:2],0,init2[3],0, init2[3 +
        # 1:(3*K)],gammas,1e-2,init2[4 + 3*K] )
        initials <- c(init2[1:2], 0, init2[3], 0, init2[4]/2, init2[4 +
            1:(3 * K - 1)], gammas, 0.01, init2[4 + 3 * K])
    } else if (distM == "zipm") {
        while (euclidean_dist(init, init2) > limits) {
            init <- init2
            tau <- tau2

            # m <-
            # nlminb(init,function(x)negQ2_G1(x,dat,distM,K,B,tau),lower
            # = c(rep(-1000,9),1e-6))
            m <- constrOptim(init, function(x) negQ2_G1(x, dat, distM,
                K, B, tau), grad = NULL, ui = bd$ui, ci = bd$ci, outer.iterations = 500,
                control = list(maxit = 50000))
            if (m$convergence != 0) {
                print(paste0("ComputeInit2=", m$convergence))
                init2 <- NA
                break
            } else {
                init2 <- m$par
            }
            switch <- compare_mu(init2, dat, distM, group = 1, X_group1,
                K)
            init2 <- switch$init2
            tau2 <- switch$tau2
            countEM <- countEM + 1
            if (K == 1) {
                break
            }
        }
        # beta0,beta1,beta3, alpha_0k, alpha_1k, psi_k, delta
        gammas <- gammas_init(dat)
        initials <- c(init2[1:2], 0, init2[3], 0, init2[3 + 1:(3 * K -
            1)], gammas, 0.01, init2[3 + 3 * K])
    }

    return(initials)
}

# log(factorial(m)) = sum(log(m)) more stable using log
loghik_zinbm <- function(m, r, p, beta0, beta1, beta2, beta3, beta4,
    delta, eta, x, y) {
    logbinom(r, m) + m * log(1 - p) - ((y - beta0 - (beta1 * m) - beta2 -
        (beta3 + beta4) * x)^2)/(2 * delta^2) - eta^2 * m
}

# log(factorial(m)) = sum(log(m)) more stable using log
loghik_zipm <- function(m, alpha0, alpha1, beta0, beta1, beta2, beta3,
    beta4, delta, eta, x, y) {
    m * (alpha0 + alpha1 * x) - ((y - beta0 - (beta1 * m) - beta2 -
        (beta3 + beta4) * x)^2)/(2 * delta^2) - eta^2 * m - sum(log(1:m))
}

# negative expectation of log-likelihood function with repect to
# conditional distribution of 1(Ci = k) given data and current
# estiamtes for gorup 2: -Q2
negQ_G2 <- function(theta, dat, distM, K, B, tauG2 = NULL, calculate_tau = F,
    calculate_ll = F) {
    # group 2
    Y_group2 <- dat$Y[which(dat$Mobs == 0)]
    X_group2 <- dat$X[which(dat$Mobs == 0)]
    theta_trans <- trans(theta, X_group2, K, distM)

    if (distM == "zilonm") {
        # integration: output is the log(integrated value)
        LogInt_k <- loghicpp_all(X_group2, Y_group2, theta_trans[["mu_ik"]],
            theta_trans[["sig"]], theta_trans[["beta0"]], theta_trans[["beta1"]],
            theta_trans[["beta2"]], theta_trans[["beta3"]], theta_trans[["beta4"]],
            theta_trans[["delta"]], theta_trans[["eta"]], nodes, weights,
            B)

        l2_ik <- cbind(-log(theta_trans[["delta"]]) - (Y_group2 - theta_trans[["beta0"]] -
            theta_trans[["beta3"]] * X_group2)^2/(2 * theta_trans[["delta"]]^2) -
            0.5 * log(2 * pi), -log(theta_trans[["delta"]]) - log(theta_trans[["sig"]]) -
            log(2 * pi) + LogInt_k)


    } else if (distM == "zinbm") {
        # summation: output is the log[sum(hik)]
        loghik_m <- list(NULL)
        for (i in 1:B) {
            loghik_m[[i]] <- loghik_zinbm(m = i, theta_trans[["r"]],
                theta_trans[["p_ik"]], theta_trans[["beta0"]], theta_trans[["beta1"]],
                theta_trans[["beta2"]], theta_trans[["beta3"]], theta_trans[["beta4"]],
                theta_trans[["delta"]], theta_trans[["eta"]], x = X_group2,
                y = Y_group2)
        }
        loghik_mmax <- NULL
        for (k in seq_len(K)) {
            loghik_mmax <- cbind(loghik_mmax, apply(sapply(loghik_m,
                function(t) t[, k]), 1, max))
        }
        output <- log(Reduce("+", lapply(loghik_m, function(t) exp(t -
            loghik_mmax)))) + loghik_mmax

        l2_ik <- cbind(-(Y_group2 - theta_trans[["beta0"]] - theta_trans[["beta3"]] *
            X_group2)^2/(2 * theta_trans[["delta"]]^2), -log(theta_trans[["p_ik"]]^(-theta_trans[["r"]]) -
            1) + output) - log(theta_trans[["delta"]]) - 0.5 * log(2 *
            pi)
    } else if (distM == "zipm") {
        # summation: output is the log[sum(hik)]
        loghik_m <- list(NULL)
        for (i in 1:B) {
            loghik_m[[i]] <- loghik_zipm(m = i, alpha0 = theta_trans[["alpha_ik"]][,
                1:K, drop = F], alpha1 = theta_trans[["alpha_ik"]][,
                K + (1:K)], theta_trans[["beta0"]], theta_trans[["beta1"]],
                theta_trans[["beta2"]], theta_trans[["beta3"]], theta_trans[["beta4"]],
                theta_trans[["delta"]], theta_trans[["eta"]], x = X_group2,
                y = Y_group2)
        }
        loghik_mmax <- NULL
        for (k in seq_len(K)) {
            loghik_mmax <- cbind(loghik_mmax, apply(sapply(loghik_m,
                function(t) t[, k]), 1, max))
        }
        output <- log(Reduce("+", lapply(loghik_m, function(t) exp(t -
            loghik_mmax)))) + loghik_mmax

        # l2_ik <- cbind(-log(theta_trans[['delta']]) -
        # (Y_group2-theta_trans[['beta0']]-theta_trans[['beta3']]*X_group2)^2/(2*theta_trans[['delta']]^2)
        # - 0.5*log(2*pi), -log(theta_trans[['delta']]) -
        # 0.5*log(2*pi) + theta_trans[['log_em1']] + output )
        l2_ik <- cbind(-(Y_group2 - theta_trans[["beta0"]] - theta_trans[["beta3"]] *
            X_group2)^2/(2 * theta_trans[["delta"]]^2), theta_trans[["log_em1"]] +
            output) - log(theta_trans[["delta"]]) - 0.5 * log(2 * pi)
    }

    bigpsi_ik <- cbind(theta_trans[["Del_i"]], (1 - theta_trans[["Del_i"]]) *
        theta_trans[["psi_ik"]])
    pow <- log(bigpsi_ik) + l2_ik
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (!calculate_ll) {
            # calculate -Q2
            Q2_ik <- tauG2 * pow
            out <- -sum(Q2_ik)
        } else {
            # calculate -l2 alternative method for calculating
            # hessian: from log likelihood function parameter
            # estimation cannot use it (has many saddle points so
            # need EM), but hessian can l2_i <- log(
            # bigpsi0*exp(l2_ik0) + bigpsi1*exp(l2_ik1) +
            # bigpsi2*exp(l2_ik2) )
            pow_max <- apply(pow, 1, max)
            l2_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l2_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci = k) given
        # data and current estiamtes
        ll <- exp(pow)
        out <- ll/rowSums(ll)
    }
    return(out)
}

# Q for all
negQ <- function(theta, dat, distM, K, B, tau = NULL, calculate_tau = F,
    calculate_ll = F) {
    if (sum(dat$Mobs == 0) > 0) {
        if (!calculate_tau) {
            if (!calculate_ll) {
                out <- negQ_G1(theta, dat, distM, K, B, tauG1 = tau$tauG1,
                  calculate_tau, calculate_ll) + negQ_G2(theta, dat,
                  distM, K, B, tauG2 = tau$tauG2, calculate_tau, calculate_ll)
            } else {
                out <- negQ_G1(theta, dat, distM, K, B, tauG1 = NULL,
                  calculate_tau, calculate_ll) + negQ_G2(theta, dat,
                  distM, K, B, tauG2 = NULL, calculate_tau, calculate_ll)
            }
            # print(theta);print(out) if (is.nan(out)) {
            # print(theta) print(out) }
        } else {
            out <- list(tauG1 = negQ_G1(theta, dat, distM, K, B = NULL,
                tauG1 = NULL, calculate_tau, calculate_ll), tauG2 = negQ_G2(theta,
                dat, distM, K, B, tauG2 = NULL, calculate_tau, calculate_ll))
        }
    } else {
        if (!calculate_tau) {
            if (!calculate_ll) {
                out <- negQ2_G1(theta, dat, distM, K, B, tauG1 = tau$tauG1,
                  calculate_tau, calculate_ll)
            } else {
                out <- negQ2_G1(theta, dat, distM, K, B, tauG1 = NULL,
                  calculate_tau, calculate_ll)
            }
            # print(theta);print(out)
        } else {
            out <- list(tauG1 = negQ2_G1(theta, dat, distM, K, B = NULL,
                tauG1 = NULL, calculate_tau, calculate_ll))
        }
    }
    return(out)
}

# for calculating observed I (EM approximation)
g <- function(x, y, dat, distM, K, B) {
    tau_f <- negQ(y, dat, distM, K, B, calculate_tau = T)
    out <- negQ(x, dat, distM, K, B, tau_f)
    return(out)
}

f_NIE2_zinbm <- function(t, K, x1, x2) {
    beta2 <- t[1]
    beta4 <- t[2]
    r <- t[3]
    alpha0_k <- t[3 + 1:K]
    alpha1_k <- t[3 + K + 1:K]
    if (K == 1) {
        psi_k <- 1
    } else {
        psi_k <- c(t[3 + 2 * K + 1:(K - 1)], 1 - sum(t[3 + 2 * K + 1:(K -
            1)]))
    }
    gamma0 <- t[3 + 3 * K]
    gamma1 <- t[4 + 3 * K]
    mu_x1k <- exp(alpha0_k + alpha1_k * x1)
    mu_x2k <- exp(alpha0_k + alpha1_k * x2)
    Delstar_x1 <- expit(gamma0 + gamma1 * x1)
    Delstar_x2 <- expit(gamma0 + gamma1 * x2)
    out <- (beta2 + beta4 * x2) * ((1 - Delstar_x2) * (1 - sum(psi_k *
        (r/(r + mu_x2k))^r)) - (1 - Delstar_x1) * (1 - sum(psi_k * (r/(r +
        mu_x1k))^r)))
    return(out)
}

effects <- function(theta, x1, x2, K, distM, calculate_se = F, vcovar = NULL,
    Group1 = F) {
    if (Group1) {
        if (distM == "zilonm") {
            theta <- c(theta[1:2], 0, theta[3], 0, theta[3 + 1:(3 *
                K)], c(-Inf, 0), Inf, theta[4 + 3 * K])
        } else if (distM == "zinbm") {

        } else if (distM == "zipm") {
            theta <- c(theta[1:2], 0, theta[3], 0, theta[3 + 1:(3 *
                K - 1)], c(-Inf, 0), Inf, theta[3 + 3 * K])
        }

    }
    x12 <- c(x1, x2)
    theta_trans <- trans(theta, x12, K, distM)

    Del_x12 <- theta_trans[["Del_i"]]
    if (distM == "zilonm") {
        sig2 <- theta_trans[["sig"]]^2
        exp_x12 <- exp(theta_trans[["mu_ik"]] + sig2/2)
        if (K == 1) {
            exp_x12 <- as.matrix(exp_x12)
        }

        m_x12 <- rowSums(theta_trans[["psi_ik"]] * exp_x12)
        NIE1 <- theta_trans[["beta1"]] * diff((1 - Del_x12) * m_x12)
        # NDE <- (x2 - x1)*(theta_trans[['beta3']] +
        # theta_trans[['beta4']]*(1 - Del_x1))
    } else if (distM == "zinbm") {
        Delstar_x12 <- theta_trans[["Delstar_i"]]
        m_x12 <- rowSums(theta_trans[["psi_ik"]] * theta_trans[["mu_ik"]])
        NIE1 <- theta_trans[["beta1"]] * diff((1 - Delstar_x12) * m_x12)
        # f_NIE2(t, K, x1, x2) NDE <- (x2 -
        # x1)*(theta_trans[['beta3']] + theta_trans[['beta4']]*(1
        # - Del_x12[1]))
    } else if (distM == "zipm") {
        Delstar_x12 <- theta_trans[["Delstar_i"]]
        m_x12 <- rowSums(theta_trans[["psi_ik"]] * theta_trans[["lambda_ik"]])
        NIE1 <- theta_trans[["beta1"]] * diff((1 - Delstar_x12) * m_x12)
        # NDE <- (x2 - x1)*(theta_trans[['beta3']] +
        # theta_trans[['beta4']]*(1 - Del_x12[1]))
    }
    NIE2 <- (theta_trans[["beta2"]] + theta_trans[["beta4"]] * x2) *
        (-diff(Del_x12))
    NIE <- NIE1 + NIE2
    out <- data.frame(eff = c(NIE1, NIE2, NIE))

    if (calculate_se == T) {
        # asymptotic variance by delta method
        if (distM == "zilonm") {
            # beta1, alpha0_k, alpha1_k, xi0, psi_k-1, gamma0,
            # gamma1
            g_NIE1 <- c(diff((1 - Del_x12) * m_x12), theta_trans[["beta1"]] *
                theta_trans[["psi_k"]] * diff((1 - Del_x12) * exp_x12),
                theta_trans[["beta1"]] * theta_trans[["psi_k"]] * diff((1 -
                  Del_x12) * exp_x12 * x12), sig2 * NIE1, theta_trans[["beta1"]] *
                  diff((exp_x12[, 1:(K - 1)] - exp_x12[, K]) * (1 -
                    Del_x12)), theta_trans[["beta1"]] * (-diff(Del_x12 *
                  (1 - Del_x12) * m_x12)), theta_trans[["beta1"]] *
                  (-diff(Del_x12 * (1 - Del_x12) * m_x12 * x12)))
            if (K == 1) {
                g_NIE1 <- g_NIE1[-5]
            }

            # beta2, beta4, gamma0, gamma1
            g_NIE2 <- c(-diff(Del_x12), x2 * (-diff(Del_x12)), (theta_trans[["beta2"]] +
                theta_trans[["beta4"]] * x2) * (-diff(Del_x12 * (1 -
                Del_x12))), (theta_trans[["beta2"]] + theta_trans[["beta4"]] *
                x2) * (-diff(Del_x12 * (1 - Del_x12) * x12)))

            g_NIE <- c(g_NIE1[1], g_NIE2[1:2], g_NIE1[1 + 1:(3 * K)],
                g_NIE1[2 + 3 * K] + g_NIE2[3], g_NIE1[3 + 3 * K] + g_NIE2[4])

            if (Group1) {
                g_NIE1 <- g_NIE <- g_NIE1[1:(length(g_NIE1) - 2)]
                V1 <- V3 <- vcovar[c(2, 3 + 1:(3 * K)), c(2, 3 + 1:(3 *
                  K))]
                NIE2_se <- NA
            } else {
                V1 <- vcovar[c(2, 5 + 1:(3 * K + 2)), c(2, 5 + 1:(3 *
                  K + 2))]
                V2 <- vcovar[c(3, 5, 5 + 3 * K + 1:2), c(3, 5, 5 + 3 *
                  K + 1:2)]
                V3 <- vcovar[c(2, 3, 5, 5 + 1:(3 * K + 2)), c(2, 3,
                  5, 5 + 1:(3 * K + 2))]
                NIE2_se <- sqrt(c(g_NIE2 %*% V2 %*% g_NIE2))
            }
            NIE1_se <- sqrt(c(g_NIE1 %*% V1 %*% g_NIE1))
            NIE_se <- sqrt(c(g_NIE %*% V3 %*% g_NIE))
        } else if (distM == "zinbm") {
            # beta1, alpha0_k, alpha1_k, psi_k-1, gamma0, gamma1
            g_NIE1 <- c(diff((1 - Delstar_x12) * m_x12), theta_trans[["beta1"]] *
                theta_trans[["psi_k"]] * diff((1 - Delstar_x12) * theta_trans[["mu_ik"]]),
                theta_trans[["beta1"]] * theta_trans[["psi_k"]] * diff((1 -
                  Delstar_x12) * theta_trans[["mu_ik"]] * x12), theta_trans[["beta1"]] *
                  diff((theta_trans[["mu_ik"]][, 1:(K - 1)] - theta_trans[["mu_ik"]][,
                    K]) * (1 - Delstar_x12)), theta_trans[["beta1"]] *
                  (-1) * diff(Delstar_x12 * (1 - Delstar_x12) * m_x12),
                theta_trans[["beta1"]] * (-1) * diff(Delstar_x12 * (1 -
                  Delstar_x12) * m_x12 * x12))
            if (K == 1) {
                g_NIE1 <- g_NIE1[-4]
            }
            V1 <- vcovar[c(2, 6 + 1:(3 * K + 1)), c(2, 6 + 1:(3 * K +
                1))]
            NIE1_se <- sqrt(c(g_NIE1 %*% V1 %*% g_NIE1))

            t <- c(theta_trans[["beta2"]], theta_trans[["beta4"]], theta_trans[["r"]],
                theta_trans[["alpha_k"]], theta_trans[["psi_k"]][-K],
                theta_trans[["gamma0"]], theta_trans[["gamma1"]])
            g_NIE2 <- numDeriv::grad(function(x) f_NIE2_zinbm(x, K,
                x1, x2), t)
            V2 <- vcovar[c(3, 5, 5 + 1:(3 * K + 2)), c(3, 5, 5 + 1:(3 *
                K + 2))]
            NIE2_se <- sqrt(c(g_NIE2 %*% V2 %*% g_NIE2))

            # beta1 # beta2,4, r # alpha0_k, alpha1_k, psi_k-1,
            # gamma0, gamma1
            g_NIE <- c(g_NIE1[1], g_NIE2[1:3], g_NIE1[1 + 1:(3 * K +
                1)] + g_NIE2[3 + 1:(3 * K + 1)])
            V3 <- vcovar[c(2, 3, 5, 5 + 1:(3 * K + 2)), c(2, 3, 5, 5 +
                1:(3 * K + 2))]
            NIE_se <- sqrt(c(g_NIE %*% V3 %*% g_NIE))
        } else if (distM == "zipm") {
            # beta1, alpha0_k, alpha1_k, psi_k-1, gamma0, gamma1
            g_NIE1 <- c(diff((1 - Delstar_x12) * m_x12), theta_trans[["beta1"]] *
                theta_trans[["psi_k"]] * diff((1 - Delstar_x12) * theta_trans[["lambda_ik"]]),
                theta_trans[["beta1"]] * theta_trans[["psi_k"]] * diff((1 -
                  Delstar_x12) * theta_trans[["lambda_ik"]] * x12),
                theta_trans[["beta1"]] * diff((theta_trans[["lambda_ik"]][,
                  1:(K - 1)] - theta_trans[["lambda_ik"]][, K]) * (1 -
                  Delstar_x12)), theta_trans[["beta1"]] * (-1) * diff(Delstar_x12 *
                  (1 - Delstar_x12) * m_x12), theta_trans[["beta1"]] *
                  (-1) * diff(Delstar_x12 * (1 - Delstar_x12) * m_x12 *
                  x12))
            if (K == 1) {
                g_NIE1 <- g_NIE1[-4]
            }

            # beta2, beta4, alpha0_k, alpha1_k, psi_k-1, gamma0,
            # gamma1
            g_NIE2 <- c(-diff(Del_x12), x2 * (-diff(Del_x12)), (theta_trans[["beta2"]] +
                theta_trans[["beta4"]] * x2) * theta_trans[["psi_k"]] *
                diff((1 - Delstar_x12) * exp(-theta_trans[["lambda_ik"]]) *
                  theta_trans[["lambda_ik"]]), (theta_trans[["beta2"]] +
                theta_trans[["beta4"]] * x2) * theta_trans[["psi_k"]] *
                diff((1 - Delstar_x12) * exp(-theta_trans[["lambda_ik"]]) *
                  theta_trans[["lambda_ik"]] * x12), (theta_trans[["beta2"]] +
                theta_trans[["beta4"]] * x2) * (-1) * diff((1 - Delstar_x12) *
                (exp(-theta_trans[["lambda_ik"]])[, 1:(K - 1)] - exp(-theta_trans[["lambda_ik"]])[,
                  K])), (theta_trans[["beta2"]] + theta_trans[["beta4"]] *
                x2) * (-1) * diff(Delstar_x12 * (1 - Del_x12)), (theta_trans[["beta2"]] +
                theta_trans[["beta4"]] * x2) * (-1) * diff(Delstar_x12 *
                (1 - Del_x12) * x12))
            if (K == 1) {
                g_NIE2 <- g_NIE2[-5]
            }

            # beta1 # beta2,4 # alpha0_k, alpha1_k, psi_k-1,
            # gamma0, gamma1
            g_NIE <- c(g_NIE1[1], g_NIE2[1:2], g_NIE1[1 + 1:(3 * K +
                1)] + g_NIE2[2 + 1:(3 * K + 1)])

            if (Group1) {
                g_NIE1 <- g_NIE <- g_NIE1[1:(length(g_NIE1) - 2)]
                V1 <- V3 <- vcovar[c(2, 5 + 1:(3 * K - 1)), c(2, 5 +
                  1:(3 * K - 1))]
                NIE2_se <- NA
            } else {
                V1 <- vcovar[c(2, 5 + 1:(3 * K + 1)), c(2, 5 + 1:(3 *
                  K + 1))]
                V2 <- vcovar[c(3, 5, 5 + 1:(3 * K + 1)), c(3, 5, 5 +
                  1:(3 * K + 1))]
                V3 <- vcovar[c(2, 3, 5, 5 + 1:(3 * K + 1)), c(2, 3,
                  5, 5 + 1:(3 * K + 1))]
                NIE2_se <- sqrt(c(g_NIE2 %*% V2 %*% g_NIE2))
            }
            NIE1_se <- sqrt(c(g_NIE1 %*% V1 %*% g_NIE1))
            NIE_se <- sqrt(c(g_NIE %*% V3 %*% g_NIE))
        }

        out$eff_se <- c(NIE1_se, NIE2_se, NIE_se)
    }
    return(out)
}

analysis <- function(dat, distM, K, B, limits, seed) {
    cinit <- init2 <- ComputeInit2(dat, distM, K, B, limits)
    # cinit <- init2 <- true2
    if (sum(dat$Mobs == 0) == 0) {
        if (distM == "zilonm") {
            index <- -c(3, 5, 3 * K + (6:8))
        } else if (distM == "zinbm") {
            index <- -c(3, 5, 3 * K + (6:8))
        } else if (distM == "zipm") {
            index <- -c(3, 5, 3 * K + (5:7))
        }
        cinit <- init2 <- init2[index]
        tau2 <- negQ(init2, dat, distM, K, B, calculate_tau = T)
        p <- length(init2)
        countEM <- NA
    } else {
        tau2 <- negQ(init2, dat, distM, K, B, calculate_tau = T)
        p <- length(init2)
        init <- rep(0, p)
        countEM <- 0
        bd <- bounds(K, group = 2, distM)
        # M-step
        while (euclidean_dist(init, init2) > limits) {
            init <- init2
            tau <- tau2

            # m2 <-
            # nlminb(init,function(x)negQ(x,dat,distM,K,B,tau),
            # lower = c(rep(-1000,6 + 2*K),rep(0,
            # K-1),-1000,-1000,1e-6,1e-6), upper = c(rep(1000,6 +
            # 2*K),rep(1, K-1),rep(1000,4)), control =
            # list(eval.max=5000,iter.max=5000))
            m <- constrOptim(init, function(x) negQ(x, dat, distM, K,
                B, tau), grad = NULL, ui = bd$ui, ci = bd$ci, outer.iterations = 500,
                control = list(maxit = 50000))
            # negQ(true,dat,distM,K,B,negQ(true,dat,distM,K,B,calculate_tau=T))
            if (m$convergence != 0) {
                print(paste0("seed=", seed, ", countEM=", countEM, ", convergence=",
                  m$convergence))
                print(m)
                init2 <- NA
                break
            } else {
                init2 <- m$par
            }
            switch <- compare_mu(init2, dat, distM, group = 2, dat$X,
                K, B)
            init2 <- switch$init2
            tau2 <- switch$tau2
            countEM <- countEM + 1
        }
    }
    est <- init2
    n <- dim(dat)[1]
    negll <- negQ(est, dat, distM, K, B, calculate_ll = T)
    AIC <- 2 * negll + 2 * p
    BIC <- 2 * negll + p * log(n)
    return(list(init = cinit, est = est, countEM = countEM, AIC = AIC,
        BIC = BIC, tau2 = tau2))
}

analysis2 <- function(dat, distM, K, B, est, tau2, x1, x2) {
    # method 1: hessian from log-likelihood function, parameter
    # estimation cannot do this
    hess <- pracma::hessian(function(x) negQ(x, dat, distM, K, B, calculate_ll = T),
        x = est)
    vcovar <- solve(hess)
    var <- diag(vcovar)
    if (sum(var < 0) > 0) {
        hess <- numDeriv::hessian(function(x) negQ(x, dat, distM, K,
            B, calculate_ll = T), x = est)
        vcovar <- solve(hess)
        var <- diag(vcovar)
        if (sum(var < 0) > 0) {
            # method 2: hessian matrix from EM approximation
            h.1 <- pracma::hessian(function(x) negQ(x, dat, distM, K,
                B, tau2), x = est)
            g2 <- function(y) {
                numDeriv::grad(function(x) g(x, y, dat, distM, K, B),
                  x = est)
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
    MedZIM <- effects(est, x1, x2, K, distM, calculate_se = T, vcovar,
        Group1)

    return(list(se = se, MedZIM = MedZIM))
}

realanalysis <- function(dat, distM, x1, x2, select_K_range, limits,
    B, seed) {
    set.seed(seed)

    t0 <- Sys.time()
    # models <- list(NULL) for (j in seq_along(select_K_range)) {
    # models[[j]] <-
    # analysis(dat,distM,K=select_K_range[j],B,limits,seed)
    # print(paste0('K = ',select_K_range[j],' completed.')) }
    klength <- length(select_K_range)
    # cl <- parallel::makeCluster(klength,outfile='') fun <-
    # c('analysis', 'bounds','compare_mu', 'ComputeInit',
    # 'ComputeInit2','gammas_init', 'euclidean_dist','expit',
    # 'k_to_ik', 'loghicpp_all', 'logbinom', 'loghik_zinbm',
    # 'loghik_zipm', 'negQ', 'negQ_G1', 'negQ_G2', 'negQ2_G1',
    # 'trans') parallel::clusterExport(cl=cl,
    # varlist=c('nodes','weights',fun),
    # envir=parent.env(environment()))
    # doParallel::registerDoParallel(cl)
    # parallel::clusterSetRNGStream(cl=cl,seed) models <-
    # foreach(i = seq_along(select_K_range), .multicombine=T,
    # .packages = c('Rcpp','stats','flexmix','MASS'), .noexport =
    # 'loghicpp_all', .combine = 'list', .errorhandling='pass')
    # %dopar% { sourceCpp('main3.cpp') modelk <-
    # analysis(dat,distM,K=select_K_range[i],B,limits,seed)
    # print(paste0('seed = ',seed,', K = ',select_K_range[i],'
    # completed.')) return(modelk) } print(models)
    # parallel::stopCluster(cl)
    models <- analysis(dat, distM, K = select_K_range, B, limits, seed)

    if (klength == 1) {
        AIC <- models$AIC
        BIC <- models$BIC
        K <- KAIC <- select_K_range
        selected_model <- models
    } else {
        AIC <- sapply(models, function(x) x$AIC)
        BIC <- sapply(models, function(x) x$BIC)
        K <- select_K_range[which.min(BIC)]
        selected_model <- models[[which.min(BIC)]]
        KAIC <- select_K_range[which.min(AIC)]
    }

    # t1 <- Sys.time()
    out <- analysis2(dat, distM, K, B, est = selected_model$est, tau2 = selected_model$tau2,
        x1, x2)
    # t2 <- Sys.time() print(paste0('EM = ', format(difftime(t1,
    # t0)), ', hessian = ', format(difftime(t2, t1))))

    results_parameters <- results(est = selected_model$est, se = out$se,
        init = selected_model$init, d = 5)
    if (distM == "zilonm") {
        paranames <- c("beta0", "beta1", "beta2", "beta3", "beta4",
            paste0("alpha", 1:K, rep(0:1, each = K)), "xi0", paste0("psi",
                1:K)[-K], "gamma0", "gamma1", "eta", "delta")
        if (sum(dat$Mobs == 0) == 0) {
            index <- -c(3, 5, 3 * K + (6:8))
            paranames <- paranames[index]
        }
    } else if (distM == "zinbm") {
        paranames <- c("beta0", "beta1", "beta2", "beta3", "beta4",
            "r", paste0("alpha", 1:K, rep(0:1, each = K)), paste0("psi",
                1:K)[-K], "gamma0", "gamma1", "eta", "delta")
        if (sum(dat$Mobs == 0) == 0) {
            index <- -c(3, 5, 3 * K + (6:8))
            paranames <- paranames[index]
        }
    } else if (distM == "zipm") {
        paranames <- c("beta0", "beta1", "beta2", "beta3", "beta4",
            paste0("alpha", 1:K, rep(0:1, each = K)), paste0("psi",
                1:K)[-K], "gamma0", "gamma1", "eta", "delta")
        if (sum(dat$Mobs == 0) == 0) {
            index <- -c(3, 5, 3 * K + (5:7))
            paranames <- paranames[index]
        }
    }
    rownames(results_parameters) <- paranames

    # mediation and direct effects
    out_MedZIM <- out$MedZIM
    results_effects <- results(out_MedZIM$eff, out_MedZIM$eff_se, d = 5)
    rownames(results_effects) <- c("NIE1", "NIE2", "NIE")

    out <- list(results_effects = results_effects, results_parameters = results_parameters,
        BIC = min(BIC), AIC = min(AIC))
    return(out)
}

nodes <- list(c(0.999999999999957, 0.999999988875665, 0.999977477192462,
    0.997514856457224, 0.951367964072747, 0.674271492248436, 0), c(0.999999999952856,
    0.999999204737115, 0.999688264028353, 0.987040560507377, 0.859569058689897,
    0.377209738164034), c(0.999999999998232, 0.999999999142705, 0.999999892781612,
    0.99999531604122, 0.999909384695144, 0.999065196455786, 0.994055506631402,
    0.973966868195677, 0.914879263264575, 0.7806074389832, 0.539146705387968,
    0.194357003324935), c(0.999999999999708, 0.999999999990394, 0.999999999789733,
    0.999999996787199, 0.999999964239081, 0.999999698894153, 0.999998017140595,
    0.999989482014818, 0.999953871005628, 0.999828822072875, 0.999451434435275,
    0.998454208767698, 0.996108665437509, 0.991126992441699, 0.981454826677335,
    0.964112164223547, 0.935160857521985, 0.88989140278426, 0.823317005506402,
    0.731018034792561, 0.610273657500639, 0.461253543939586, 0.287879932742716,
    0.0979238852878323), c(0.999999999999886, 0.999999999999271, 0.999999999995825,
    0.999999999978455, 0.999999999899278, 0.999999999570788, 0.999999998323362,
    0.999999993964134, 0.999999979874503, 0.999999937554078, 0.999999818893713,
    0.999999507005719, 0.999998735471866, 0.99999693244919, 0.999992937876663,
    0.999984519902271, 0.999967593067943, 0.999935019925082, 0.99987486504878,
    0.999767971599561, 0.999584750351518, 0.999281111921792, 0.998793534298806,
    0.998033336315434, 0.996880318128192, 0.995176026155327, 0.992716997196827,
    0.989248431090134, 0.984458831167431, 0.977976235186665, 0.969366732896917,
    0.958136022710214, 0.943734786052757, 0.925568634068613, 0.903013281513574,
    0.875435397630409, 0.842219246350757, 0.802798741343241, 0.75669390863373,
    0.703550005147142, 0.643176758985205, 0.575584490635152, 0.501013389379309,
    0.419952111278447, 0.333142264577638, 0.241566319538884, 0.146417984290588,
    0.0490559673050779), c(0.99999999999993, 0.999999999999817, 0.999999999999537,
    0.999999999998861, 0.999999999997274, 0.999999999993646, 0.999999999985569,
    0.999999999968033, 0.999999999930888, 0.999999999854056, 0.999999999698754,
    0.999999999391777, 0.999999998797983, 0.999999997673237, 0.999999995585634,
    0.999999991786456, 0.999999985003076, 0.999999973113236, 0.999999952642664,
    0.999999918004795, 0.999999860371215, 0.999999766023332, 0.99999961398855,
    0.999999372707335, 0.99999899541069, 0.999998413810965, 0.999997529623805,
    0.999996203347166, 0.999994239627617, 0.999991368448345, 0.999987221282001,
    0.999981301270121, 0.999972946425232, 0.999961284807857, 0.999945180614459,
    0.999923170129289, 0.999893386547593, 0.999853472773111, 0.999800481431138,
    0.999730761519808, 0.9996398313456, 0.999522237651217, 0.999371401140938,
    0.999179448934886, 0.998937034833512, 0.998633148640677, 0.998254916171996,
    0.997787391958906, 0.997213347043469, 0.996513054640254, 0.995664076816953,
    0.994641055712511, 0.993415513169264, 0.991955663002678, 0.990226240467528,
    0.988188353800743, 0.985799363025283, 0.983012791481101, 0.979778275800616,
    0.976041560256577, 0.971744541565487, 0.966825370312356, 0.961218615151116,
    0.954855495805023, 0.947664190615153, 0.939570223933275, 0.930496937997153,
    0.920366053031953, 0.90909831816302, 0.896614254280076, 0.882834988244669,
    0.867683175775646, 0.851084007987849, 0.832966293919411, 0.813263608502974,
    0.791915492376142, 0.768868686768246, 0.744078383547347, 0.717509467487324,
    0.689137725061668, 0.65895099174335, 0.626950208051043, 0.593150353591953,
    0.557581228260778, 0.52028805069123, 0.481331846116905, 0.440789599033901,
    0.398754150467238, 0.355333825165075, 0.310651780552846, 0.264845076583448,
    0.218063473469712, 0.170467972382011, 0.122229122201558, 0.0735251229856713,
    0.0245397635746492), c(0.999999999999945, 0.999999999999911, 0.999999999999856,
    0.999999999999769, 0.999999999999632, 0.999999999999419, 0.999999999999088,
    0.99999999999858, 0.999999999997803, 0.999999999996623, 0.999999999994845,
    0.999999999992181, 0.999999999988217, 0.999999999982353, 0.999999999973737,
    0.99999999996115, 0.999999999942877, 0.999999999916506, 0.99999999987867,
    0.999999999824698, 0.999999999748146, 0.999999999640173, 0.999999999488718,
    0.999999999277422, 0.99999999898421, 0.999999998579458, 0.999999998023619,
    0.999999997264174, 0.999999996231727, 0.999999994835051, 0.999999992954804,
    0.999999990435633, 0.999999987076265, 0.999999982617173, 0.999999976725267,
    0.99999996897499, 0.999999958825101, 0.999999945590273, 0.999999928406514,
    0.999999906189268, 0.999999877582855, 0.999999840899736, 0.999999794047876,
    0.999999734444233, 0.999999658912159, 0.999999563560232, 0.999999443639729,
    0.999999293377668, 0.999999105781996, 0.999998872415162, 0.999998583131984,
    0.999998225777311, 0.999997785838635, 0.999997246048426, 0.999996585930571,
    0.999995781284928, 0.999994803603649, 0.999993619412539, 0.999992189530434,
    0.999990468239214, 0.999988402356833, 0.999985930205475, 0.999982980466756,
    0.999979470915751, 0.999975307025492, 0.999970380433589, 0.999964567262626,
    0.999957726286078, 0.999949696931689, 0.999940297114479, 0.999929320891884,
    0.999916535933951, 0.999901680802003, 0.999884462029766, 0.999864551001635,
    0.999841580623482, 0.999815141782273, 0.999784779591651, 0.999749989421649,
    0.999710212711751, 0.999664832567665, 0.999613169143347, 0.999554474811097,
    0.999487929123823, 0.999412633574952, 0.99932760616283, 0.999231775767899,
    0.99912397635239, 0.999002940993729, 0.998867295764364, 0.998715553472209,
    0.998546107277413, 0.998357224202664, 0.998147038555748, 0.997913545284577,
    0.997654593286378, 0.997367878694235, 0.997050938165609, 0.996701142198913,
    0.996315688505661, 0.9958915954671, 0.995425695705623, 0.994914629802639,
    0.994354840195921, 0.993742565290762, 0.993073833820584, 0.992344459493918,
    0.991550035965892, 0.990685932173615, 0.989747288075987, 0.988729010839602,
    0.987625771513511, 0.986432002236608, 0.985141894022398, 0.983749395166731,
    0.982248210324941, 0.980631800305449, 0.978893382627494, 0.977025932891058,
    0.975022187007306, 0.972874644337964, 0.9705755717919, 0.968117008926856,
    0.96549077410362, 0.962688471739078, 0.959701500703337, 0.956521063904581,
    0.953138179103394, 0.949543690995936, 0.945728284602632, 0.941682499995769,
    0.937396748395719, 0.93286132966124, 0.928066451194555, 0.923002248276554,
    0.917658805841549, 0.912026181694487, 0.906094431166404, 0.89985363319617,
    0.893293917818211, 0.886405495026979, 0.879178684979389, 0.871603949486378,
    0.863671924734105, 0.855373455164271, 0.846699628431464, 0.8376418113436,
    0.828191686679357, 0.818341290764098, 0.808083051673338, 0.797409827920312,
    0.78631494747182, 0.774792246924435, 0.762836110661392, 0.75044150979924,
    0.737604040722825, 0.724319962997407, 0.710586236438006, 0.696400557108459,
    0.681761392016427, 0.666668012265738, 0.651120524424304, 0.635119899864422,
    0.618668001832728, 0.601767610009633, 0.584422442322664, 0.566637173785019,
    0.548417452139792, 0.529769910101775, 0.510702174002558, 0.491222868660811,
    0.471341618317998, 0.451069043500452, 0.430416753691437, 0.409397335721529,
    0.388024337812118, 0.366312249234904, 0.344276475579705, 0.321933309653369,
    0.29929989806396, 0.276394203576179, 0.253234963356, 0.229841643254361,
    0.206234388311029, 0.182433969690289, 0.158461728289299, 0.134339515287672,
    0.110089629932628, 0.085734754877651, 0.06129788941366, 0.0368022809500251,
    0.0122713551180822))
weights = list(c(1.35817842745391e-12, 2.1431204556943e-07, 0.000266200513752717,
    0.0183431669899278, 0.230022394514789, 0.965976579412301, 1.5707963267949),
    c(1.16311658142558e-09, 1.19837013631707e-05, 0.00290251774790131,
        0.0763857435708323, 0.531078275428054, 1.38961475924726), c(4.93785387766319e-11,
        1.86872822687364e-08, 1.82633205937107e-06, 6.24825592407441e-05,
        0.000949946804283469, 0.00774260102606424, 0.0391750054936008,
        0.137422107733168, 0.360461418469344, 0.737437848361548, 1.19346302584916,
        1.52328371863471), c(8.6759314149796e-12, 2.52163479185301e-10,
        4.87600609742406e-09, 6.58351851271834e-08, 6.47775660359297e-07,
        4.82371820326155e-06, 2.81101643279401e-05, 0.0001320523412561,
        0.000513393824067903, 0.00169087399814264, 0.00481629814392846,
        0.012083543599158, 0.027133510013712, 0.0552896837422406, 0.103432154223333,
        0.179324412110728, 0.290240679312454, 0.440833236273858, 0.630405135164744,
        0.85017285645662, 1.08163498549007, 1.29747575042498, 1.46601442671697,
        1.55877335553333), c(3.48419376702611e-12, 2.09893354045115e-11,
        1.13060553474947e-10, 5.4828357797095e-10, 2.40917732564759e-09,
        9.64988889610896e-09, 3.54347771714219e-08, 1.19924427829028e-07,
        3.75954118623606e-07, 1.09688351259013e-06, 2.99166158781388e-06,
        7.65957585252031e-06, 1.84818135998792e-05, 4.21831838417576e-05,
        9.13908174907101e-05, 0.000188564429767003, 0.000371666936216778,
        0.000701859515684242, 0.00127332794470824, 0.00222508270647864,
        0.00375425097743183, 0.00613003763208303, 0.00970722373939169,
        0.0149378350960501, 0.0223794710636485, 0.032698732726609, 0.0466682080548466,
        0.0651555334325362, 0.0891031392409415, 0.119497411288696, 0.157326203484366,
        0.203523998858602, 0.258904639514054, 0.324082539611529, 0.399384741525717,
        0.484758091214755, 0.579678103087788, 0.683068516344264, 0.793242700820517,
        0.907879379154895, 1.02404493311181, 1.13827224337631, 1.24670120745186,
        1.34527888476625, 1.4300083548723, 1.49722622254104, 1.54388111617696,
        1.56778143130722), c(2.18359220992336e-12, 5.51823694681749e-12,
        1.35425129123363e-11, 3.23044643332524e-11, 7.49673975738182e-11,
        1.69394577894116e-10, 3.72995018430528e-10, 8.00997844797297e-10,
        1.67888976821619e-09, 3.43718567446501e-09, 6.8784610955899e-09,
        1.3464645522302e-08, 2.57995682295359e-08, 4.84209501980724e-08,
        8.90713951402424e-08, 1.60693945790762e-07, 2.84499236591598e-07,
        4.94582887027542e-07, 8.44737563848599e-07, 1.41830671554939e-06,
        2.34216672085281e-06, 3.80619832646449e-06, 6.0899100320949e-06,
        9.59819412837847e-06, 1.49085140318706e-05, 2.28321181090361e-05,
        3.44921247593432e-05, 5.14214974476588e-05, 7.56839965862015e-05,
        0.000110021128466667, 0.000158027884007012, 0.000224359652050085,
        0.000314972091860212, 0.000437394956159117, 0.000601039879911474,
        0.000817541013324695, 0.00110112611345194, 0.00146901435994298,
        0.00194183577598437, 0.00254406576752917, 0.00330446699403483,
        0.00425652959901786, 0.005438899797624, 0.006895785969066, 0.00867733074953918,
        0.0108399371682559, 0.0134465366052857, 0.0165667862542476,
        0.0202771838175001, 0.0246610873147533, 0.0298086281173101,
        0.0358165056041964, 0.0427876521577257, 0.0508307575725705,
        0.0600596423586363, 0.070592469906867, 0.0825507881107017, 0.0960583918651895,
        0.111239998988745, 0.128219733631201, 0.147119413257857, 0.168056637948269,
        0.191142684133427, 0.216480209117296, 0.24416077786984, 0.274262229689068,
        0.306845909417917, 0.341953795923017, 0.379605569386652, 0.419795668445015,
        0.462490398055368, 0.507625158831908, 0.555101878003633, 0.604786730578404,
        0.656508246131628, 0.710055901205469, 0.765179298908956, 0.821588035266965,
        0.878952345552782, 0.936904612745668, 0.995041804046133, 1.05292887995527,
        1.11010319396534, 1.16607986993243, 1.22035810957936, 1.27242834553786,
        1.32178011744377, 1.3679105116809, 1.41033297144626, 1.44858625496132,
        1.48224329788554, 1.51091972307417, 1.5342817381543, 1.55205316984541,
        1.56402140377323, 1.57004202927959), c(1.72376440360427e-12,
        2.76085876713983e-12, 4.38886518997793e-12, 6.92552801526814e-12,
        1.08491618343371e-11, 1.6874519343915e-11, 2.60619605028053e-11,
        3.99733895199303e-11, 6.08933870643807e-11, 9.21405642265189e-11,
        1.38502875258341e-10, 2.06842035390292e-10, 3.06927020787233e-10,
        4.52575761041538e-10, 6.63208634701621e-10, 9.65948518580993e-10,
        1.39844143125654e-09, 2.01262102188672e-09, 2.87970134644711e-09,
        4.09675791396852e-09, 5.79534957309661e-09, 8.15274648285765e-09,
        1.14064655548505e-08, 1.58729780909682e-08, 2.19716489170893e-08,
        3.02551964648995e-08, 4.14482335585054e-08, 5.64957638716706e-08,
        7.66238739748176e-08, 1.03415280407456e-07, 1.3890286996036e-07,
        1.85684913700251e-07, 2.47066245112497e-07, 3.27230373305176e-07,
        4.31448255871654e-07, 5.66330284020508e-07, 7.40128934909086e-07,
        9.63100521176592e-07, 1.24793551211775e-06, 1.61026800945077e-06,
        2.06927612573647e-06, 2.64838622540434e-06, 3.37609523480551e-06,
        4.28692649399982e-06, 5.42253589182403e-06, 6.83298627742185e-06,
        8.57820935370791e-06, 1.07296754068523e-05, 1.33722922845324e-05,
        1.66065559765083e-05, 2.05509759449342e-05, 2.53447989688501e-05,
        3.11510556774483e-05, 3.81599541202458e-05, 4.65926446304946e-05,
        5.67053798539329e-05, 6.87940931134964e-05, 8.31994172400365e-05,
        0.000100312164601124, 0.00012057928729057, 0.000144510334290946,
        0.000172684419885929, 0.000205757714679988, 0.000244471467286892,
        0.000289660561088772, 0.000342262606462977, 0.000403327564549541,
        0.000474027894018167, 0.000555669207425785, 0.00064970141867429,
        0.000757730357827359, 0.00088152982417297, 0.00102305404297454,
        0.00118445048589028, 0.0013680730096097, 0.00157649526191056,
        0.00181252429912886, 0.00207921435400867, 0.00237988068810044,
        0.00271811345835023, 0.00309779152330082, 0.00352309611044299,
        0.00399852426273353, 0.00452890197915629, 0.00511939696145378,
        0.00577553087680655, 0.00650319104428258, 0.00730864145131325,
        0.0081985330052612, 0.0091799129243109, 0.0102602331714108,
        0.0114473578347998, 0.0127495693587312, 0.0141755735283305,
        0.0157345031130598, 0.0174359200739775, 0.0192898162408456,
        0.0213066123661261, 0.0234971554639885, 0.0258727143436176,
        0.0284449732473342, 0.0312260235053317, 0.0342283531201757,
        0.037464834195629, 0.0409487081258673, 0.0446935684627655, 0.0487133413807014,
        0.0530222636602872, 0.0576348581146547, 0.0625659063844551,
        0.067830419030658, 0.0734436028576388, 0.0794208254030082, 0.0857775765352682,
        0.0925294271057785, 0.0996919846077886, 0.107280845802556, 0.115311546280937,
        0.123799506938413, 0.132759977352446, 0.142207976063402, 0.1521582277742,
        0.162625097499384, 0.173622521711631, 0.185163936552775, 0.197262203197437,
        0.209929530480208, 0.223177394922185, 0.237016458319419, 0.25145648308451,
        0.266506245563135, 0.282173447579596, 0.298464626499445, 0.315385064132678,
        0.332938694837748, 0.351128013224425, 0.369953981892114, 0.389415939679212,
        0.409511510938194, 0.430236516389789, 0.451584886147545, 0.473548575540614,
        0.496117484397316, 0.519279380484246, 0.543019827824843, 0.567322120646739,
        0.592167223728207, 0.617533719929909, 0.643397765708253, 0.669733055410291,
        0.696510795146547, 0.72369968702683, 0.751265924524357, 0.779173199704798,
        0.807382723018763, 0.835853256308267, 0.864541159619669, 0.893400452347196,
        0.922382889152455, 0.951438051016353, 0.980513451680855, 1.00955465962943,
        1.03850543563739, 1.0673078857975, 1.09590262979297, 1.1242289840506,
        1.15222515926255, 1.17982847161733, 1.20697556693131, 1.23360265672195,
        1.25964576511667, 1.28504098534673, 1.30972474443744, 1.33363407457576,
        1.35670688951561, 1.37888226427314, 1.40010071626944, 1.42030448599969,
        1.43943781524641, 1.45744722081255, 1.47428176172808, 1.4898932978833,
        1.50423673806368, 1.51727027540505, 1.52895560835458, 1.53925814531188,
        1.54814719123556, 1.55559611463166, 1.56158249349181, 1.56608823891746,
        1.56909969535167, 1.57060771653828))
