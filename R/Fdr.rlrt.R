Fdr.rlrt <- function(rlrt.obj, threshold) {
	rstat = pmax(0, rlrt.obj$stat)
	rsim = rlrt.obj$sim
	tbar = mean(rsim)
	p.mom = 1 - 3 * tbar^2 / mean(rsim^2)
	a.mom = tbar / (1 - p.mom)
	pi0 = mean(rstat == 0) / p.mom
	cat("Estimated null proportion", pi0, "\n")
	Fdr1 = pi0 * (1 - p.mom) * (1 - pchisq(threshold / a.mom, 1)) / mean(rstat >= threshold)
	Fdr2 = pi0 * mean(rsim >= threshold) / mean(rstat >= threshold)
	otpt = c(Fdr1, Fdr2)
	names(otpt) = c("MoM", "ML")
	otpt
}

