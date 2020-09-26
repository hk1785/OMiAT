MiCAM <-
function(Y, otu.tab, cov=NULL, tax.tab, tree, model=c("gaussian", "binomial"), tax.levels=c("phylum to genus", "kingdom to species"),
multi.corr=c("BY", "BH"), pow=c(1:4, Inf), g.unif.alpha=c(0.5), n.perm=3000,
filename.rel.abundance="rel.abundance.txt", filename.unadj.pvs="unadj.pvalues.txt", 
filename.adj.pvs="adj.pvalues.txt", filename.sig.assoc.taxa="sig.assoc.taxa.txt", filename.graph="graph.pdf") {	

	if (!identical(colnames(otu.tab), tree$tip.label)) {
		stop("The OTU IDs and tip labels are not matched! Check if identical(colnames(otu.tab), tree$tip.label) is TRUE!")
	}
	if (!is.rooted(tree)) {
		stop("The tree is not rooted!")
	}
	
	if (length(Y) != nrow(otu.tab)) {
		otu.tab <- t(otu.tab)
	}
    if (!is.null(cov)) {
		cov <- as.data.frame(cov)
		rownames(cov) <- 1:nrow(otu.tab)
	}
	rownames(otu.tab) <- 1:nrow(otu.tab)
	colnames(otu.tab) <- sub('New.ReferenceOTU', 'N', colnames(otu.tab))
	rownames(tax.tab) <- colnames(otu.tab)
	tree$tip.label <- colnames(otu.tab)
	
	tax.tab[is.na(tax.tab[,1]),1] <- "k__"
	tax.tab[is.na(tax.tab[,2]),2] <- "p__"
	tax.tab[is.na(tax.tab[,3]),3] <- "c__"
	tax.tab[is.na(tax.tab[,4]),4] <- "o__"
	tax.tab[is.na(tax.tab[,5]),5] <- "f__"
	tax.tab[is.na(tax.tab[,6]),6] <- "g__"
	tax.tab[is.na(tax.tab[,7]),7] <- "s__"
	
	clade.names <- colnames(tax.tab)
	rownames(tax.tab) <- 1:nrow(tax.tab)
	pvs <- vector("list", length(clade.names)) 
	
	ind0 <- vector("list", length(names(table(tax.tab[,1]))))
	for (i in 1:length(names(table(tax.tab[,1])))) {
		ind <- which(tax.tab[,1] == names(table(tax.tab[,1]))[[i]])
		if (!is.matrix(tax.tab[,1])) ind0[i] <- list(as.numeric(names(tax.tab[ind,1]))) 
		if (is.matrix(tax.tab[,1])) ind0[i] <- list(as.numeric(rownames(tax.tab[,1][ind,])))
	}
	ind1 <- list()
	for (j in 1:length(ind0)) {
		for (i in 1:length(names(table(tax.tab[ind0[[j]],2])))) {
			ind <- which(tax.tab[ind0[[j]],2] == names(table(tax.tab[ind0[[j]],2]))[i])
			ind.out <- list(as.numeric(rownames(tax.tab[ind0[[j]],][ind,]))) 
			ind1 <- c(ind1, ind.out)
		}
	}
	ind2 <- list()
	for (j in 1:length(ind1)) {
		for (i in 1:length(names(table(tax.tab[ind1[[j]],3])))) {
			ind <- which(tax.tab[ind1[[j]],3] == names(table(tax.tab[ind1[[j]],3]))[i])
			ind.out <- list(as.numeric(rownames(tax.tab[ind1[[j]],][ind,]))) 
			ind2 <- c(ind2, ind.out)
		}
	}
	ind3 <- list()
	for (j in 1:length(ind2)) {
		for (i in 1:length(names(table(tax.tab[ind2[[j]],4])))) {
			ind <- which(tax.tab[ind2[[j]],4] == names(table(tax.tab[ind2[[j]],4]))[i])
			ind.out <- list(as.numeric(rownames(tax.tab[ind2[[j]],][ind,]))) 
			ind3 <- c(ind3, ind.out)
		}
	}
	ind4 <- list()
	for (j in 1:length(ind3)) {
		for (i in 1:length(names(table(tax.tab[ind3[[j]],5])))) {
			ind <- which(tax.tab[ind3[[j]],5] == names(table(tax.tab[ind3[[j]],5]))[i])
			ind.out <- list(as.numeric(rownames(tax.tab[ind3[[j]],][ind,]))) 
			ind4 <- c(ind4, ind.out)
		}
	}
	ind5 <- list()
	for (j in 1:length(ind4)) {
		for (i in 1:length(names(table(tax.tab[ind4[[j]],6])))) {
			ind <- which(tax.tab[ind4[[j]],6] == names(table(tax.tab[ind4[[j]],6]))[i])
			ind.out <- list(as.numeric(rownames(tax.tab[ind4[[j]],][ind,])))
			ind5 <- c(ind5, ind.out)
		}
	}
	ind6 <- list()
	for (j in 1:length(ind5)) {
		for (i in 1:length(names(table(tax.tab[ind5[[j]],7])))) {
			ind <- which(tax.tab[ind5[[j]],7] == names(table(tax.tab[ind5[[j]],7]))[i])
			ind.out <- list(as.numeric(rownames(tax.tab[ind5[[j]],][ind,]))) 
			names(ind.out[[1]]) <- colnames(otu.tab)[list(as.numeric(rownames(tax.tab[ind5[[j]],][ind,])))[[1]]]
			ind6 <- c(ind6, ind.out)
		}
	}
	ind.all <- list(ind0, ind1, ind2, ind3, ind4, ind5, ind6)  
	
	ind.s <- unlist(ind.all[[7]])
	direction <- rep(NA, length(ind.s))
	
    if (is.null(cov)) {
        r <- Y - mean(Y)
    } else {
        fit <- glm(Y ~ ., family = model, data = as.data.frame(cov))
        res <- Y - fitted.values(fit)
        r <- res - mean(res)
    }
	
	for (j in 1:length(ind.s)) {
		s.ind <- ind.s[j]
		sub.otu.tab <- otu.tab[,s.ind]
		if (is.null(cov)) {
			r <- Y - mean(Y)
			U <- as.vector(t(scale(sub.otu.tab)) %*% r)
		} else {	
			fit <- glm(Y ~ ., family = model, data=cov)
			res <- Y - fitted.values(fit)				
			r <- res - mean(res)
			U <- as.vector(t(scale(sub.otu.tab)) %*% r)
		}
		direction[j] <- U
	}
	ind.1 <- which(direction >= 0)
	ind.2 <- which(direction < 0)
	direction[ind.1] <- "(+)"
	direction[ind.2] <- "(-)"
	
	otu.numbers <- NULL
	for (i in 1:length(ind6)) {
		otu.numbers.out <- names(ind6[[i]])
		otu.numbers <- c(otu.numbers,otu.numbers.out)
	}
	otu.codes <- paste(otu.numbers, sep=" ")
	otu.codes <- paste(otu.codes, direction, sep=" ")
	
	otu.rr.tax <- matrix(NA, length(otu.numbers),3)
	for (i in 1:length(otu.numbers)) {
		ind <- which(colnames(otu.tab) == otu.numbers[i])
		relative.abundance <- round(sum(otu.tab[,ind])/sum(otu.tab),5)
		name <- paste(tax.tab[ind,1],tax.tab[ind,2],tax.tab[ind,3],tax.tab[ind,4],tax.tab[ind,5],tax.tab[ind,6],tax.tab[ind,7],sep="|")
		otu.rr.tax[i,] <- cbind(otu.numbers[i],relative.abundance,name)
	}
	colnames(otu.rr.tax) <- c("OTU ID", "Relative Abundance", "Taxonomy")
	write.table(as.data.frame(otu.rr.tax), file=paste(filename.rel.abundance), sep="\t", row.names=FALSE, quote=FALSE)
	
	if (is.null(cov)) {
		fit.glm <- glm(Y ~ 1, family = model)
	} else {
		fit.glm <- glm(Y ~ ., family = model, data = cov)
	}
	
	total.reads <- rowSums(otu.tab)
	prop <- apply(otu.tab/total.reads, 2, scale)
		
	res <- Y - fitted.values(fit.glm)
	r <- res - mean(res)
	r <- jitter(r)
	r.s <- list()
	for (j in 1:n.perm) {
		r.s[[j]] <- r[shuffle(length(r))]
	}
	
	if (tax.levels == "kingdom to species") {
		for (i in 1:length(clade.names)) {
			ind.s <- ind.all[[i]]
			for (j in 1:length(ind.s)) {
				s.ind <- ind.s[[j]]
				sub.otu.tab <- otu.tab[,s.ind]
				sub.prop <- prop[,s.ind]
				names.surv <- colnames(sub.prop)
				if (length(s.ind) == 1) {
					U <- as.vector(t(sub.prop) %*% r)
					U0 <- unlist(lapply(r.s,function(x) return(as.vector(t(sub.prop) %*% x))))
					omiat.out <- (sum(abs(U0) > abs(U))+1)/(n.perm+1)
				} else {
					sub.tree <- prune_taxa(names.surv, tree)
					U <- as.vector(t(sub.prop) %*% r)
					Ts = rep(NA, length(pow))
					Ts[which(pow==Inf)] <- max(U)
					Ts[which(pow!=Inf)] <- unlist(lapply(as.list(pow[which(pow!=Inf)]),function(x) return(sum(U^x))))
					U0 <- lapply(r.s,function(x) return(as.vector(t(sub.prop) %*% x)))
					T0s <- list()
					spu.pvs <- rep(NA, length(pow))
					for (j in 1:length(pow)) {
						T0s[[j]] <- T0s.e <- sapply(U0,function(x) if (pow[j] < Inf) return(sum(x^pow[j])) else return(max(x)))
						spu.pvs[j] <- (sum(abs(T0s.e) > abs(Ts[j]))+1)/(n.perm+1)
					}	
					T.aspu <- min(spu.pvs)	
					T0.aspu <- rep(NA, n.perm)
					for (l in 1:n.perm) {
						T0s.n <- list()
						for (m in 1:length(pow)) {
							T0s.n[[m]] <- T0s[[m]][-l]
						}
						a.U <-  as.vector(t(sub.prop) %*% r.s[[l]])
						a.pvs <- rep(NA, length(pow))
						a.Ts = rep(NA, length(pow))
						a.Ts[which(pow==Inf)] <- max(a.U)
						a.Ts[which(pow!=Inf)] <- unlist(lapply(as.list(pow[which(pow!=Inf)]),function(x) return(sum(a.U^x))))
						a.pvs <- unlist(mapply(function(x,y)(sum(abs(x) > abs(y))+1)/(n.perm+1),T0s.n,a.Ts))
						T0.aspu[l] <- min(a.pvs)
					}
					p.aspu <- (sum(T0.aspu < T.aspu)+1)/(n.perm+1)

					total.reads <- rowSums(otu.tab)
					unifs <- GUniFrac2(sub.otu.tab, sub.tree, alpha = c(g.unif.alpha, 1), total.reads = total.reads)$unifracs
					bray.curtis <- as.matrix(bcdist(sub.otu.tab))
					ind <- which(rowSums(sub.otu.tab) == 0)
					suppressWarnings(jac <- as.matrix(distance(sub.otu.tab, method = "jaccard")))
					u.unif <- unifs[,,"d_UW"]
					w.unif <- unifs[,,"d_1"]
					g.unif <- list()
					for (k in 1:length(g.unif.alpha)) {
						g.unif[[k]] <- unifs[,,paste("d_",g.unif.alpha[k],sep="")]
					}
					if (sum(is.na(u.unif)) > 0 | sum(is.na(jac)) > 0) {
						bray.curtis.kern <- D2K(bray.curtis)
						Q <- as.numeric(t(r)%*%bray.curtis.kern%*%r)
						Q0 <- rep(NA, n.perm)
						for (j in 1:n.perm) {
							Q0[j] <- t(r.s[[j]])%*%bray.curtis.kern%*%r.s[[j]]
						}
						Q.omni <- p.bc <- (sum(abs(Q0) > abs(Q))+1)/(n.perm+1)
						mirkat.pvs <- c(p.bc, p.bc)
						S.aomiat0 <- apply(cbind(T0.aspu, Q0),1,min)
						S.aomiat <- min(T.aspu, Q.omni)
						omiat.out <- (sum(S.aomiat0 < S.aomiat)+1)/(n.perm+1)	
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
							Qs[j] <- t(r)%*%list.kernels[[j]]%*%r
						}
						Q0s <- list()
						for (j in 1:length(list.kernels)) {
							Q0s.inv <- rep(NA, n.perm)
							for (k in 1:n.perm) {
								Q0s.inv[k] <- t(r.s[[k]])%*%list.kernels[[j]]%*%r.s[[k]]
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
							a.Qs <- unlist(lapply(list.kernels,function(x) return(t(r.s[[l]])%*%x%*%r.s[[l]])))
							a.pvs <- unlist(mapply(function(x,y) (sum(abs(x) > abs(y))+1)/(n.perm+1),Q0s.n,a.Qs))
							Q0.omni[l] <- min(a.pvs)
						}
						p.omni <- (sum(Q0.omni < Q.omni)+1)/(n.perm+1)
						S.aomiat0 <- apply(cbind(T0.aspu, Q0.omni),1,min)
						S.aomiat <- min(T.aspu, Q.omni)
						omiat.out <- (sum(S.aomiat0 < S.aomiat)+1)/(n.perm+1)					
					}
				}
				names(omiat.out) <- tax.tab[s.ind,i][1]
				pvs[[i]] <- c(pvs[[i]], omiat.out)
			}
		}	
		names(pvs) <- clade.names
	}
	
	if (tax.levels == "phylum to genus") {
		for (i in 1:length(clade.names)) {
			ind.s <- ind.all[[i]]
			for (j in 1:length(ind.s)) {
				s.ind <- ind.s[[j]]
				sub.otu.tab <- otu.tab[,s.ind]
				sub.prop <- prop[,s.ind]
				names.surv <- colnames(sub.prop)
				if (length(s.ind) == 1) {
					U <- as.vector(t(sub.prop) %*% r)
					U0 <- unlist(lapply(r.s,function(x) return(as.vector(t(sub.prop) %*% x))))
					omiat.out <- (sum(abs(U0) > abs(U))+1)/(n.perm+1)
				} else {
					sub.tree <- prune_taxa(names.surv, tree)
					U <- as.vector(t(sub.prop) %*% r)
					Ts = rep(NA, length(pow))
					Ts[which(pow==Inf)] <- max(U)
					Ts[which(pow!=Inf)] <- unlist(lapply(as.list(pow[which(pow!=Inf)]),function(x) return(sum(U^x))))
					U0 <- lapply(r.s,function(x) return(as.vector(t(sub.prop) %*% x)))
					T0s <- list()
					spu.pvs <- rep(NA, length(pow))
					for (j in 1:length(pow)) {
						T0s[[j]] <- T0s.e <- sapply(U0,function(x) if (pow[j] < Inf) return(sum(x^pow[j])) else return(max(x)))
						spu.pvs[j] <- (sum(abs(T0s.e) > abs(Ts[j]))+1)/(n.perm+1)
					}	
					T.aspu <- min(spu.pvs)	
					T0.aspu <- rep(NA, n.perm)
					for (l in 1:n.perm) {
						T0s.n <- list()
						for (m in 1:length(pow)) {
							T0s.n[[m]] <- T0s[[m]][-l]
						}
						a.U <-  as.vector(t(sub.prop) %*% r.s[[l]])
						a.pvs <- rep(NA, length(pow))
						a.Ts = rep(NA, length(pow))
						a.Ts[which(pow==Inf)] <- max(a.U)
						a.Ts[which(pow!=Inf)] <- unlist(lapply(as.list(pow[which(pow!=Inf)]),function(x) return(sum(a.U^x))))
						a.pvs <- unlist(mapply(function(x,y)(sum(abs(x) > abs(y))+1)/(n.perm+1),T0s.n,a.Ts))
						T0.aspu[l] <- min(a.pvs)
					}
					p.aspu <- (sum(T0.aspu < T.aspu)+1)/(n.perm+1)

					total.reads <- rowSums(otu.tab)
					unifs <- GUniFrac2(sub.otu.tab, sub.tree, alpha = c(g.unif.alpha, 1), total.reads = total.reads)$unifracs
					bray.curtis <- as.matrix(bcdist(sub.otu.tab))
					suppressWarnings(jac <- as.matrix(distance(sub.otu.tab, method = "jaccard")))
					u.unif <- unifs[,,"d_UW"]
					w.unif <- unifs[,,"d_1"]
					g.unif <- list()
					for (k in 1:length(g.unif.alpha)) {
						g.unif[[k]] <- unifs[,,paste("d_",g.unif.alpha[k],sep="")]
					}
					if (sum(is.na(u.unif)) > 0 | sum(is.na(jac)) > 0) {
						bray.curtis.kern <- D2K(bray.curtis)
						Q <- as.numeric(t(r)%*%bray.curtis.kern%*%r)
						Q0 <- rep(NA, n.perm)
						for (j in 1:n.perm) {
							Q0[j] <- t(r.s[[j]])%*%bray.curtis.kern%*%r.s[[j]]
						}
						Q.omni <- p.bc <- (sum(abs(Q0) > abs(Q))+1)/(n.perm+1)
						mirkat.pvs <- c(p.bc, p.bc)
						S.aomiat0 <- apply(cbind(T0.aspu, Q0),1,min)
						S.aomiat <- min(T.aspu, Q.omni)
						omiat.out <- (sum(S.aomiat0 < S.aomiat)+1)/(n.perm+1)	
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
							Qs[j] <- t(r)%*%list.kernels[[j]]%*%r
						}
						Q0s <- list()
						for (j in 1:length(list.kernels)) {
							Q0s.inv <- rep(NA, n.perm)
							for (k in 1:n.perm) {
								Q0s.inv[k] <- t(r.s[[k]])%*%list.kernels[[j]]%*%r.s[[k]]
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
							a.Qs <- unlist(lapply(list.kernels,function(x) return(t(r.s[[l]])%*%x%*%r.s[[l]])))
							a.pvs <- unlist(mapply(function(x,y) (sum(abs(x) > abs(y))+1)/(n.perm+1),Q0s.n,a.Qs))
							Q0.omni[l] <- min(a.pvs)
						}
						p.omni <- (sum(Q0.omni < Q.omni)+1)/(n.perm+1)
						S.aomiat0 <- apply(cbind(T0.aspu, Q0.omni),1,min)
						S.aomiat <- min(T.aspu, Q.omni)
						omiat.out <- (sum(S.aomiat0 < S.aomiat)+1)/(n.perm+1)					
					}
				}
				names(omiat.out) <- tax.tab[s.ind,i][1]
				pvs[[i]] <- c(pvs[[i]], omiat.out)
			}
		}	
		names(pvs) <- clade.names
	}
	
	if (tax.levels == "kingdom to species") {
		if (multi.corr == "BH") {
			adj.pvs <- pvs
			for (i in 1:length(pvs)) {
				p <- pvs[[i]]
				adj.pvs[[i]] <- p.adjust(p, method = "BH", n = length(p))
			}
		}
		if (multi.corr == "BY") {
			adj.pvs <- pvs
			for (i in 1:length(pvs)) {
				p <- pvs[[i]]
				adj.pvs[[i]] <- p.adjust(p, method = "BY", n = length(p))
			}
		}
	}
	
	if (tax.levels == "phylum to genus") {
		if (multi.corr == "BH") {
			adj.pvs <- pvs
			for (i in 2:(length(pvs)-1)) {
				p <- pvs[[i]]
				adj.pvs[[i]] <- p.adjust(p, method = "BH", n = length(p))
			}
		}
		if (multi.corr == "BY") {
			adj.pvs <- pvs
			for (i in 2:(length(pvs)-1)) {
				p <- pvs[[i]]
				adj.pvs[[i]] <- p.adjust(p, method = "BY", n = length(p))
			}
		}
	}

	count.table <- vector("list", length(clade.names)) 
	for (i in 1:length(ind0)) {
		count <- length(ind0[[i]])
		count.table[[1]] <- c(count.table[[1]], count)
	}
	for (i in 1:length(ind1)) {
		count <- length(ind1[[i]])
		count.table[[2]] <- c(count.table[[2]], count)
	}
	for (i in 1:length(ind2)) {
		count <- length(ind2[[i]])
		count.table[[3]] <- c(count.table[[3]], count)
	}
	for (i in 1:length(ind3)) {
		count <- length(ind3[[i]])
		count.table[[4]] <- c(count.table[[4]], count)
	}
	for (i in 1:length(ind4)) {
		count <- length(ind4[[i]])
		count.table[[5]] <- c(count.table[[5]], count)
	}
	for (i in 1:length(ind5)) {
		count <- length(ind5[[i]])
		count.table[[6]] <- c(count.table[[6]], count)
	}
	for (i in 1:length(ind6)) {
		count <- length(ind6[[i]])
		count.table[[7]] <- c(count.table[[7]], count)
	}

	adj.count.table <- matrix(NA, length(clade.names), nrow(tax.tab))
	for (i in 1:length(clade.names)) {
		adj.count.table[i,] <- c(count.table[[i]],rep(0,nrow(tax.tab)-length(count.table[[i]])))
	}

	if (tax.levels == "kingdom to species") {
		name.table <- vector("list", length(clade.names)) 
		for (i in 1:length(clade.names)) {
			name.table[[i]] <- names(adj.pvs[[i]])
		}
		names(name.table) <- clade.names
		
		adj.name.table <- matrix(NA, length(clade.names), nrow(tax.tab))
		for (i in 1:length(clade.names)) {
			adj.name.table[i,] <- c(name.table[[i]],rep(0,nrow(tax.tab)-length(name.table[[i]])))
		}
		
		col.table <-  matrix("black", length(clade.names), nrow(tax.tab))
		for (i in 1:length(clade.names)) {
			for (j in 1:length(adj.pvs[[i]])) {
				if (adj.pvs[[i]][j] < 0.05) {
					col.table[i,j] <- "red2"
				} else if (adj.pvs[[i]][j] >= 0.05) {
					col.table[i,j] <- "lightgrey"
				} else if (is.na(adj.pvs[[i]][j])) {
					col.table[i,j] <- "white"
				}
			}
		}
		unadj.pvs.table <-  matrix(NA, length(clade.names)*2, nrow(tax.tab))
		for (i in 1:length(clade.names)) {
			for (j in 1:length(pvs[[i]])) {
				unadj.pvs.table[2*i-1,j] <- names(pvs[[i]][j])
				unadj.pvs.table[2*i,j] <- pvs[[i]][j]
			}
		}
		rownames.unadj.pvs.table <- rep(NA, length(clade.names)*2)
		for (i in 1:length(clade.names)) {
			rownames.unadj.pvs.table[2*i-1] <- clade.names[i]
			rownames.unadj.pvs.table[2*i] <- "Unadj. pvs"
		}
		rownames(unadj.pvs.table) <- rownames.unadj.pvs.table
		write.table(as.data.frame(unadj.pvs.table), file=paste(filename.unadj.pvs), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

		adj.pvs.table <-  matrix(NA, length(clade.names)*2, nrow(tax.tab))
		for (i in 1:length(clade.names)) {
			for (j in 1:length(adj.pvs[[i]])) {
				adj.pvs.table[2*i-1,j] <- names(adj.pvs[[i]][j])
				adj.pvs.table[2*i,j] <- adj.pvs[[i]][j]
			}
		}
		rownames.adj.pvs.table <- rep(NA, length(clade.names)*2)
		for (i in 1:length(clade.names)) {
			rownames.adj.pvs.table[2*i-1] <- clade.names[i]
			rownames.adj.pvs.table[2*i] <- "Adj. pvs"
		}
		rownames(adj.pvs.table) <- rownames.adj.pvs.table
		write.table(as.data.frame(adj.pvs.table), file=paste(filename.adj.pvs), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)	

		assoc.taxa.names.table <- matrix(NA, 7, 1)
		ind <- which(as.numeric(adj.pvs.table[2,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[1,1] <- paste(adj.pvs.table[1,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[4,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[2,1] <- paste(adj.pvs.table[3,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[6,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[3,1] <- paste(adj.pvs.table[5,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[8,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[4,1] <- paste(adj.pvs.table[7,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[10,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[5,1] <- paste(adj.pvs.table[9,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[12,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[6,1] <- paste(adj.pvs.table[11,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[14,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[7,1] <- paste(adj.pvs.table[13,ind], collapse = "|")
		assoc.taxa.names.table <- as.data.frame(assoc.taxa.names.table)
		rownames(assoc.taxa.names.table) <- clade.names
		colnames(assoc.taxa.names.table) <- "Significantly associated taxa"
		write.table(assoc.taxa.names.table, file=paste(filename.sig.assoc.taxa), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	
		pdf(file=paste(filename.graph), width=3.4, height=round((ncol(otu.tab)*0.12),1))
		layout(matrix(c(1,2,3,4,5,6,7,8), ncol=8, byrow=FALSE), widths=c(5,5,5,5,5,5,5,9))
		par(mai=c(0.3,0.05,0.15,0))
		par(oma=c(0,0,0,0))
		for (i in 1:length(clade.names)) {
			barplot(t(t(adj.count.table[i,])), names.arg = clade.names[i], cex.names=0.8, las=1, horiz=FALSE, axes=FALSE, col=col.table[i,])
		}
		axis(4, seq(0.5,nrow(tax.tab)-0.5,1), otu.codes, las=2, cex.axis=0.6)
		par(mai=c(0,0,0,0))
		plot.new()
		legend(x="bottomleft", c(expression("adj. P" < "0.05"), expression("adj. P " >= " 0.05")), fill=c("red2", "lightgrey"), cex=0.65, horiz=F, bty='n')
		graphics.off() 
	}
	
	if (tax.levels == "phylum to genus") {
		name.table <- vector("list", length(clade.names)) 
		for (i in 2:(length(clade.names)-1)) {
			name.table[[i]] <- names(adj.pvs[[i]])
		}
		names(name.table) <- clade.names
		
		adj.name.table <- matrix(NA, length(clade.names), nrow(tax.tab))
		for (i in 2:(length(clade.names)-1)) {
			adj.name.table[i,] <- c(name.table[[i]],rep(0,nrow(tax.tab)-length(name.table[[i]])))
		}

		col.table <-  matrix("black", length(clade.names), nrow(tax.tab))
		for (i in 2:(length(clade.names)-1)) {
			for (j in 1:length(adj.pvs[[i]])) {
				if (adj.pvs[[i]][j] < 0.05) {
					col.table[i,j] <- "red2"
				} else if (adj.pvs[[i]][j] >= 0.05) {
					col.table[i,j] <- "lightgrey"
				} else if (is.na(adj.pvs[[i]][j])) {
					col.table[i,j] <- "white"
				}
			}
		}	
		
		unadj.pvs.table <-  matrix(NA, length(clade.names)*2, nrow(tax.tab))
		for (i in 2:(length(clade.names)-1)) {
			for (j in 1:length(pvs[[i]])) {
				unadj.pvs.table[2*i-1,j] <- names(pvs[[i]][j])
				unadj.pvs.table[2*i,j] <- pvs[[i]][j]
			}
		}
		rownames.unadj.pvs.table <- rep(NA, length(clade.names)*2)
		for (i in 2:(length(clade.names)-1)) {
			rownames.unadj.pvs.table[2*i-1] <- clade.names[i]
			rownames.unadj.pvs.table[2*i] <- "Unadj. pvs"
		}
		rownames(unadj.pvs.table) <- rownames.unadj.pvs.table
		unadj.pvs.table <- unadj.pvs.table[-c(1,2,13,14),]
		ind <- which(apply(unadj.pvs.table,2,function(x)sum(is.na(x)))==10)
		unadj.pvs.table <- unadj.pvs.table[,-ind]
		write.table(as.data.frame(unadj.pvs.table), file=paste(filename.unadj.pvs), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

		adj.pvs.table <-  matrix(NA, length(clade.names)*2, nrow(tax.tab))
		for (i in 2:(length(clade.names)-1)) {
			for (j in 1:length(adj.pvs[[i]])) {
				adj.pvs.table[2*i-1,j] <- names(adj.pvs[[i]][j])
				adj.pvs.table[2*i,j] <- adj.pvs[[i]][j]
			}
		}
		rownames.adj.pvs.table <- rep(NA, length(clade.names)*2)
		for (i in 2:(length(clade.names)-1)) {
			rownames.adj.pvs.table[2*i-1] <- clade.names[i]
			rownames.adj.pvs.table[2*i] <- "Adj. pvs"
		}
		rownames(adj.pvs.table) <- rownames.adj.pvs.table
		adj.pvs.table <- adj.pvs.table[-c(1,2,13,14),]
		ind <- which(apply(adj.pvs.table,2,function(x)sum(is.na(x)))==10)
		adj.pvs.table <- adj.pvs.table[,-ind]
		write.table(as.data.frame(adj.pvs.table), file=paste(filename.adj.pvs), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
		
		assoc.taxa.names.table <- matrix(NA, 5, 1)
		ind <- which(as.numeric(adj.pvs.table[2,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[1,1] <- paste(adj.pvs.table[1,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[4,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[2,1] <- paste(adj.pvs.table[3,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[6,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[3,1] <- paste(adj.pvs.table[5,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[8,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[4,1] <- paste(adj.pvs.table[7,ind], collapse = "|")
		ind <- which(as.numeric(adj.pvs.table[10,1:ncol(adj.pvs.table)])<0.05)
		assoc.taxa.names.table[5,1] <- paste(adj.pvs.table[9,ind], collapse = "|")
		rownames(assoc.taxa.names.table) <- clade.names[-c(1,7)]
		colnames(assoc.taxa.names.table) <- "Significantly associated taxa"
		write.table(assoc.taxa.names.table, file=paste(filename.sig.assoc.taxa), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
		
		pdf(file=paste(filename.graph), width=3.4, height=round((ncol(otu.tab)*0.12),1))
		layout(matrix(c(1,2,3,4,5,6), ncol=6, byrow=FALSE), widths=c(5,5,5,5,5,9))
		par(mai=c(0.3,0.05,0.15,0))
		par(oma=c(0,0,0,0))
		for (i in 2:(length(clade.names)-1)) {
			barplot(t(t(adj.count.table[i,])), names.arg = clade.names[i], cex.names=0.8, las=1, horiz=FALSE, axes=FALSE, col=col.table[i,])
		}
		axis(4, seq(0.5,nrow(tax.tab)-0.5,1), otu.codes, las=2, cex.axis=0.6)
		par(mai=c(0,0,0,0))
		plot.new()
		legend(x="bottomleft", c(expression("adj. P" < "0.05"), expression("adj. P " >= " 0.05")), fill=c("red2", "lightgrey"), cex=0.65, horiz=F, bty='n')
		graphics.off()
	}		
}
