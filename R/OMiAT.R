OMiAT <-
function (Y, otu.tab, cov=NULL, tree, total.reads=NULL, model = c("gaussian", "binomial"), pow=c(1:4, Inf), g.unif.alpha=c(0.5), n.perm=3000) {

	if (!identical(colnames(otu.tab), tree$tip.label)) {
		stop("The OTU IDs and tip labels are not matched! Check if identical(colnames(otu.tab), tree$tip.label) is TRUE!")
	}
	if (!is.rooted(tree)) {
		stop("The tree is not rooted!")
	}
	
    if (length(Y) != nrow(otu.tab)) {
        otu.tab <- t(otu.tab)
    }
    if (is.null(total.reads)) {
		total.reads <- rowSums(otu.tab)
    }
    prop <- apply(otu.tab/total.reads, 2, scale)
	
    if (is.null(cov)) {
        r <- Y - mean(Y)
    } else {
        fit <- glm(Y ~ ., family = model, data = as.data.frame(cov))
        res <- Y - fitted.values(fit)
        r <- res - mean(res)
    }
	
    r <- jitter(r)
    r.s <- list()
    for (j in 1:n.perm) {
        r.s[[j]] <- r[shuffle(length(r))]
    }
    U <- as.vector(t(prop) %*% r)
    Ts = rep(NA, length(pow))
    Ts[which(pow == Inf)] <- max(U)
    Ts[which(pow != Inf)] <- unlist(lapply(as.list(pow[which(pow != Inf)]), function(x) return(sum(U^x))))
    U0 <- lapply(r.s, function(x) return(as.vector(t(prop) %*% x)))
    T0s <- list()
    pvs <- rep(NA, length(pow))
    for (j in 1:length(pow)) {
        T0s[[j]] <- T0s.e <- sapply(U0, function(x) if (pow[j] < Inf) return(sum(x^pow[j])) else return(max(x)))
        pvs[j] <- (sum(abs(T0s.e) > abs(Ts[j]))+1)/(n.perm+1)
    }
	
    T.aspu <- min(pvs)
    T0.aspu <- rep(NA, n.perm)
    for (l in 1:n.perm) {
        T0s.n <- list()
        for (m in 1:length(pow)) {
            T0s.n[[m]] <- T0s[[m]][-l]
        }
        a.U <- as.vector(t(prop) %*% r.s[[l]])
        a.pvs <- rep(NA, length(pow))
        a.Ts <- rep(NA, length(pow))
        a.Ts[which(pow == Inf)] <- max(a.U)
        a.Ts[which(pow != Inf)] <- unlist(lapply(as.list(pow[which(pow != Inf)]), function(x) return(sum(a.U^x))))
        a.pvs <- unlist(mapply(function(x, y) (sum(abs(x) > abs(y))+1)/(n.perm+1), T0s.n, a.Ts))
        T0.aspu[l] <- min(a.pvs)
    }
    p.aspu <- (sum(T0.aspu < T.aspu)+1)/(n.perm+1)
    Ts <- c(Ts, T.aspu)
    names(Ts) <- c(paste("SPU(", pow, ")", sep = ""), "aSPU")
    spu.pvs <- c(pvs, p.aspu)
    names(spu.pvs) <- c(paste("SPU(", pow, ")", sep = ""), "aSPU")
	
    if (is.null(total.reads)) {
        unifs <- GUniFrac(otu.tab, tree, alpha = c(g.unif.alpha, 1))$unifracs
		jac <- as.matrix(distance(otu.tab, method = "jaccard"))
        bray.curtis <- as.matrix(bcdist(otu.tab))
        u.unif <- unifs[, , "d_UW"]
        w.unif <- unifs[, , "d_1"]
        g.unif <- list()
    } else {
        unifs <- GUniFrac2(otu.tab, tree, alpha = c(g.unif.alpha, 1), total.reads = total.reads)$unifracs
		jac <- as.matrix(distance(otu.tab, method = "jaccard"))
        bray.curtis <- as.matrix(bcdist(otu.tab))
        u.unif <- unifs[, , "d_UW"]
        w.unif <- unifs[, , "d_1"]
        g.unif <- list()
    }
    for (k in 1:length(g.unif.alpha)) {
        g.unif[[k]] <- unifs[, , paste("d_", g.unif.alpha[k], sep = "")]
    }
	
    if (sum(is.na(u.unif)) > 0) {
        jac.kern <- D2K(jac)
		bray.curtis.kern <- D2K(bray.curtis)
        list.kernels <- list(jac.kern = jac.kern, bray.curtis.kern = bray.curtis.kern)
        Qs <- rep(NA, length(list.kernels))
        for (j in 1:length(list.kernels)) {
                Qs[j] <- t(r) %*% list.kernels[[j]] %*% r
        }
        Q0s <- list()
        for (j in 1:length(list.kernels)) {
            Q0s.inv <- rep(NA, n.perm)
            for (k in 1:n.perm) {
                  Q0s.inv[k] <- t(r.s[[k]]) %*% list.kernels[[j]] %*% r.s[[k]]
            }
            Q0s[[j]] <- Q0s.inv
        }
        mirkat.pvs <- rep(NA, length(list.kernels))
        for (j in 1:length(list.kernels)) {
            mirkat.pvs[j] <- (sum(abs(Q0s[[j]]) > abs(Qs[[j]]))+1)/(n.perm+1)
        }
		
        Q.omni <- min(mirkat.pvs)
        Q0.omni <- rep(NA, n.perm)
        for (l in 1:n.perm) {
            Q0s.n <- list()
            for (m in 1:length(list.kernels)) {
                Q0s.n[[m]] <- Q0s[[m]][-l]
            }
			a.Qs <- unlist(lapply(list.kernels, function(x) return(t(r.s[[l]]) %*% x %*% r.s[[l]])))
            a.pvs <- unlist(mapply(function(x, y) (sum(abs(x) > abs(y))+1)/(n.perm+1), Q0s.n, a.Qs))
            Q0.omni[l] <- min(a.pvs)
        }
        p.omni <- (sum(Q0.omni < Q.omni)+1)/(n.perm+1)
        Qs <- c(Qs, Q.omni)
        names(Qs) <- c("Jaccard", "Bray-Curtis", "Opt.MiRKAT")
        mirkat.pvs <- c(mirkat.pvs, p.omni)
        names(mirkat.pvs) <- c("Jaccard", "Bray-Curtis", "Opt.MiRKAT")
        M.omiat0 <- apply(cbind(T0.aspu, Q0.omni), 1, min)
        M.omiat <- min(T.aspu, Q.omni)
        p.omiat <- (sum(M.omiat0 < M.omiat)+1)/(n.perm+1)
        names(M.omiat) <- "OMiAT"
        names(p.omiat) <- "OMiAT"
    } else {
        jac.kern <- D2K(jac)
		bray.curtis.kern <- D2K(bray.curtis)
        u.unif.kern <- D2K(u.unif)
        w.unif.kern <- D2K(w.unif)
        g.unif.kern <- list()
        for (k in 1:length(g.unif.alpha)) {
            g.unif.kern[[k]] <- D2K(g.unif[[k]])
        }
        list.kernels <- c(list(jac.kern = jac.kern, bray.curtis.kern = bray.curtis.kern, u.unif.kern = u.unif.kern, w.unif.kern = w.unif.kern), g.unif.kern)
        Qs <- rep(NA, length(list.kernels))
        for (j in 1:length(list.kernels)) {
                Qs[j] <- t(r) %*% list.kernels[[j]] %*% r
        }
        Q0s <- list()
        for (j in 1:length(list.kernels)) {
            Q0s.inv <- rep(NA, n.perm)
            for (k in 1:n.perm) {
                  Q0s.inv[k] <- t(r.s[[k]]) %*% list.kernels[[j]] %*% r.s[[k]]
            }
            Q0s[[j]] <- Q0s.inv
        }
        mirkat.pvs <- rep(NA, length(list.kernels))
        for (j in 1:length(list.kernels)) {
            mirkat.pvs[j] <- (sum(abs(Q0s[[j]]) > abs(Qs[[j]]))+1)/(n.perm+1)
        }
		
        Q.omni <- min(mirkat.pvs)
        Q0.omni <- rep(NA, n.perm)
        for (l in 1:n.perm) {
            Q0s.n <- list()
            for (m in 1:length(list.kernels)) {
                Q0s.n[[m]] <- Q0s[[m]][-l]
            }
			a.Qs <- unlist(lapply(list.kernels, function(x) return(t(r.s[[l]]) %*% x %*% r.s[[l]])))
            a.pvs <- unlist(mapply(function(x, y) (sum(abs(x) > abs(y))+1)/(n.perm+1), Q0s.n, a.Qs))
            Q0.omni[l] <- min(a.pvs)
        }
        p.omni <- (sum(Q0.omni < Q.omni)+1)/(n.perm+1)
        Qs <- c(Qs, Q.omni)
        names(Qs) <- c("Jaccard", "Bray-Curtis", "U.UniFrac", "W.UniFrac", paste("G.UniFrac(", g.unif.alpha, ")", sep = ""), "Opt.MiRKAT")
        mirkat.pvs <- c(mirkat.pvs, p.omni)
        names(mirkat.pvs) <- c("Jaccard", "Bray-Curtis", "U.UniFrac", "W.UniFrac", paste("G.UniFrac(", g.unif.alpha, ")", sep = ""), "Opt.MiRKAT")
        M.omiat0 <- apply(cbind(T0.aspu, Q0.omni), 1, min)
        M.omiat <- min(T.aspu, Q.omni)
        p.omiat <- (sum(M.omiat0 < M.omiat)+1)/(n.perm+1)
        names(M.omiat) <- "OMiAT"
        names(p.omiat) <- "OMiAT"
    }
    return(list(SPU.pvs = spu.pvs, MiRKAT.pvs = mirkat.pvs, OMiAT.pvalue = p.omiat))
}
