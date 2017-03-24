library(MCMCpack)
library(DAAG)

sim_expression_poisson = function(n_cell, n, lambda) {
	expression_per_cell = matrix(rpois(n_cell*length(lambda),lambda),ncol=length(lambda),byrow=TRUE)
	reads_per_cell = matrix(nrow=n_cell,ncol=length(lambda))
	for (i in 1:n_cell) {
		reads_per_cell[i,] = rmultinom(1,n[i],expression_per_cell[i,])
	}
	return(list(expression=expression_per_cell,reads=reads_per_cell)) }


sim_expression_arbitrary = function(n_cell, n, pi_k) {
	expression_per_cell = matrix(nrow=n_cell,ncol=nrow(pi_k))
	for (i in 1:nrow(pi_k)) {
		expression_per_cell[,i] = sample(0:(ncol(pi_k)-1),n_cell,replace=TRUE,prob=pi_k[i,])
	}
	reads_per_cell = matrix(nrow=n_cell,ncol=length(lambda))
	for (i in 1:n_cell) {
		reads_per_cell[i,] = rmultinom(1,n[i],expression_per_cell[i,])
	}
	return(list(expression=expression_per_cell,reads=reads_per_cell)) }

approximate_EM = function(dat, num_iter = 10,k_plus = 100, lambda = rep(1,ncol(dat)),eps = .1) {
	num_cell = nrow(dat)
	num_gene = ncol(dat)
	num_read = rowSums(dat)
	#initialize
	#cur_pi = rdirichlet(num_gene,rep(alpha,k_plus+1))
	cur_pi = matrix(dpois(sapply(0:k_plus,rep,num_gene),lambda),ncol=k_plus+1)
	#iterate
	k = matrix(0:k_plus,byrow=TRUE,nrow=num_cell,ncol=k_plus+1)
	T = 0
	varT = 0
	for (iter in 1:num_iter) {
		pi_list = list()
		pi_list[[1]] = cur_pi
		for (EMstep in 2:3) {
			#approximate the total number of transcripts using the pis
			E = cur_pi%*%k[1,] #expectation of expression for every gene
			E2 = cur_pi%*%k[1,]^2 #expectation of expression squared for every gene
			varE = E2-E^2 #variance of expression for every gene
			varT = sum(varE)
			T = sum(E)
			#loop over genes
			for (j in 1:num_gene) {
				if (j%%(num_gene/20) == 1) {
					cat("Iteration", iter, "EM step", EMstep-1, "gene", j, "\r")
					flush.console()
				}
				#to deal with dumb vectorization...
				pi_mat = matrix(log(cur_pi[j,]), nrow=num_cell, ncol = k_plus+1, byrow=TRUE)

				#fix the expression of the current gene
				Tjk = T- E[j] + k
				varTjk = varT - varE[j]
				
				#This version has no penalty
				#logLike = dat[,j]*log(k)+(num_read-dat[,j])*log(T-k)+pi_mat
				
				#This penalty comes from Taylor expansion on E(log(T-K))
				#logLike = dat[,j]*log(k)+(num_read-dat[,j])*log(T-k)+pi_mat-(num_read-dat[,j])/(2*(T-k)^2)*varT
				
				#This penalty comes from Taylor expansion on E((T-k)^(n-r)) and then assuming the log goes through
				#logLike = dat[,j]*log(k)+(num_read-dat[,j])*log(T-k)+pi_mat+(num_read-dat[,j]-2)*log((T-k))
				
				#this formula comes from actually explicitly including E(1/T) and approximating with Taylor series
				#ETinv = 1/T+1/T^3*varT
				#logLike = dat[,j]*log(k*ETinv)+(num_read-dat[,j])*log(1-k*ETinv)+pi_mat
				
				#This version comes from expanding E(binomial)
				Ck= exp(log(dat[,j]) + log(Tjk^2-4*Tjk*k+2*k^2+(Tjk-2*k)^2*dat[,j]) + log(varTjk) - log(2) - 2*log(Tjk) - 2*log(Tjk-k) )
				penalty = log(1+Ck)
				logLike = dat[,j]*log(k/Tjk) + (num_read-dat[,j])*log(1-k/Tjk)+penalty+pi_mat
				
				#replace the nans with 0s which they should be
				bad_0 = which(is.na(logLike[,1]))
				logLike[bad_0] = (num_read*log(T))[bad_0]
				#compute posterior
				post = exp(logLike-apply(logLike,1,max))
				post = post/rowSums(post)
				cur_pi[j,] = colSums(post)/sum(post)
			}
			print(cur_pi)
			pause()
			pi_list[[EMstep]] = cur_pi
		}
		r = pi_list[[2]]-pi_list[[1]]
		v = pi_list[[3]]-2*pi_list[[2]]+pi_list[[1]]
		step_size = sum(r^2)/sum(v^2)
		cur_pi = pi_list[[1]]-step_size*(pi_list[[2]]-pi_list[[1]])
		print(pi_list[[1]])
		print(pi_list[[2]])
		print(cur_pi)
		pause()
	}
	return(cur_pi)
}
