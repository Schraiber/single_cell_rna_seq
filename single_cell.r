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

#Simplex projection method of Condat (2016) "Fast projection onto the simplex and the l_1 ball
proj_simplex = function(y) {
	v = list()
	v = y[1]
	vtilde = c()
	rho = y[1] - 1
	for (i in 2:length(y)) {
		rho = rho + (y[i]-rho)/(length(v) + 1)
		if (rho > y[i] - 1) {
			v = unlist(list(v,y[i]))
		} else {
			vtilde = unlist(list(vtilde,v))
			v = y[i]
			rho = y[i] - 1
		}
	}
	if (length(vtilde) > 0) {
		for (i in 1:length(vtilde)) {
			if (vtilde[i] > rho) {
				v = unlist(list(v,vtilde[i]))
				rho = rho + (vtilde[i]-rho)/length(v)
			}
		}
	}
	newV = v
	while (TRUE) {
		oldV = newV
		newV = c()
		num_removed = 0
		for (i in 1:length(oldV)) {
			if (oldV[i] <= rho) {
				num_removed = num_removed + 1
				rho = rho + (rho - oldV[i])/(length(oldV)-num_removed)
			} else {
				newV = unlist(list(newV, oldV[i]))
			}
		}
		if (!num_removed) { break } 
		
	}
	tau = rho
	K = length(newV)
	x = numeric(length(y))
	for (i in 1:length(y)) {
		x[i] = max(y[i] - tau, 0)
	}
	return(x)
}

approximate_EM = function(dat, num_iter = 10,k_plus = 100, accel_iter = .5*num_iter,lambda = rep(1,ncol(dat)),eps = .1) {
	pi_per_iteration = list()
	num_cell = nrow(dat)
	num_gene = ncol(dat)
	num_read = rowSums(dat)
	#initialize
	#cur_pi = rdirichlet(num_gene,rep(alpha,k_plus+1))
	#cur_pi = matrix(dpois(sapply(0:k_plus,rep,num_gene),lambda),ncol=k_plus+1)
	cur_pi = matrix(1/(k_plus+1),ncol=k_plus+1,nrow=num_gene)
	#iterate
	k = matrix(0:k_plus,byrow=TRUE,nrow=num_cell,ncol=k_plus+1)
	T = 0
	varT = 0
	for (iter in 1:num_iter) {
		pi_list = list()
		pi_list[[1]] = cur_pi
		#TODO: for the first few iterations, just do a normal EM?
		if (iter > accel_iter) {
			extra_EM = 1
		} else {
			extra_EM = 0
		}
		for (EMstep in 2:(2+extra_EM)) {
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
			pi_list[[EMstep]] = cur_pi
		}
		r = pi_list[[2]]-pi_list[[1]]
		if (extra_EM == 1) {
			v = pi_list[[3]]-2*pi_list[[2]]+pi_list[[1]]
			step_size = sum(r^2)/sum(v^2)
			if (step_size < 1) step_size = 1 #never take a step smaller than the EM would take!
		} else {
			v = 0
			step_size = .5 #to counteract the 2 in the SQUAREM update
		}
		cat("\n")
		print(c("current step is", step_size))
		#looks like + gives the right thing, b/c if stepsize = 1, you get cur_pi = pi_list[[2]]
		#cur_pi = pi_list[[1]]+step_size*(pi_list[[2]]-pi_list[[1]])
		cur_pi = pi_list[[1]]+2*step_size*r+step_size^2*v
		print(sum(cur_pi<0))
		if (extra_EM == 1) {
			cur_pi = t(apply(cur_pi,1,proj_simplex))
		}
		print(sum(cur_pi<0))
		#Shrink the step til none escape the bounds
		#while (sum(cur_pi<0) > 0) {
		#	step_size = (step_size+1)/2
			#print(c("current step is", step_size))
			#cur_pi = pi_list[[1]]+step_size*(pi_list[[2]]-pi_list[[1]])
		#	cur_pi = pi_list[[1]]+2*step_size*r+step_size^2*v
		#} 
		#print(c("current step is", step_size))
		#some may have overstepped...
		#TODO: should they overstep? Is this the right thing to do?
		#TODO: 0 is a bad idea. Makes it impossible to get info for that one...
		#TODO: maybe make it 1e-100?
		#cur_pi[cur_pi<0] = 1e-100
		pi_per_iteration[[iter]] = cur_pi
	}
	cat("\n")
	return(pi_per_iteration)
}
