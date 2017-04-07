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

approximate_EM = function(dat, num_iter = 10,k_plus = 100, accel_iter = .5*num_iter,lambda = rep(1,ncol(dat)),eps = .1,start_pi = NULL) {
	pi_per_iteration = list()
	num_cell = nrow(dat)
	num_gene = ncol(dat)
	num_read = rowSums(dat)
	#initialize
	if (is.null(start_pi)) { 
		#alpha = 1
		#start_pi = rdirichlet(num_gene,rep(alpha,k_plus+1))
		#start_pi = matrix(dpois(sapply(0:k_plus,rep,num_gene),lambda),ncol=k_plus+1)
		#start_pi = matrix(1/(k_plus+1),ncol=k_plus+1,nrow=num_gene)
		start_pi = matrix(dnbinom(sapply(0:k_plus,rep,num_gene),size=1,mu=2),ncol=k_plus+1)
	} else {
		if (any(dim(start_pi) != c(num_gene,k_plus+1))) { 
			stop("Dimension of start_pi not = (num_gene, k_plus+1)")
		}
	}
	cur_pi = start_pi	

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
				if (j%%(num_gene/20) == 0) {
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
				#ETinv = 1/Tjk+1/T^3*varTjk
				#logLike = dat[,j]*log(k*ETinv)+(num_read-dat[,j])*log(1-k*ETinv)
				#penalty = 0
				
				#This version comes from expanding E(binomial)
				penalty = log( 1 +  (k^2*num_read*(num_read+1)-2*Tjk*k*num_read*(dat[,j]+1) + Tjk^2*dat[,j]*(1+dat[,j]))/ 
					(2*Tjk^2*(Tjk-k)^2)*varTjk)
					
				logLike = dbinom(dat[,j], num_read,k/Tjk,log=TRUE) + penalty

				#compute posterio
				logPost = logLike + pi_mat
					
				#compute posterior
				post = exp(logPost-apply(logPost,1,max))
				post = post/rowSums(post)
	
				#update pi
				cur_pi[j,] = colSums(post)/num_cell
			}
			pi_list[[EMstep]] = cur_pi
		}
		r = pi_list[[2]]-pi_list[[1]]
		if (extra_EM == 1) {
			v = pi_list[[3]]-2*pi_list[[2]]+pi_list[[1]]
			step_size = sum(r^2)/sum(v^2)
			#step_size[step_size<1] = 1 #never take a step smaller than the EM would take!
			#step_size[step_size>5] = 5 #don't go toooooooo far
		} else {
			v = 0
			step_size = .5 #to counteract the 2 in the SQUAREM update
		}
		cat("\n")
		print(c("current step summary is", summary(step_size)))
		#cur_pi = pi_list[[1]]+step_size*(pi_list[[2]]-pi_list[[1]])
		cur_pi = pi_list[[1]]+2*step_size*r+step_size^2*v
		if (extra_EM == 1) {
			#Shrink the step til none escape the bounds
			#print(pi_list[[1]][1,])
			#print(cur_pi[1,])
			#while (sum(cur_pi<0) > 0) {
			#	step_size = (step_size+1)/2
			#	#cur_pi = pi_list[[1]]+step_size*(pi_list[[2]]-pi_list[[1]])
			#	cur_pi = pi_list[[1]]+2*step_size*r+step_size^2*v
			#}
			#project onto Simplex
			print(cur_pi[626,])
			cur_pi = t(apply(cur_pi,1,proj_simplex))
			print(cur_pi[626,])
			#set any that are 0 to 1e-100
			#NB: this loses the normalization...
			#cur_pi[cur_pi == 0] = 1e-100
		}
		print(c("current mean step is", mean(step_size)))
		pi_per_iteration[[iter]] = cur_pi
		if (iter > 1) {
			par_change = sqrt(sum(pi_per_iteration[[iter]]-pi_per_iteration[[iter-1]])^2)
		} else {
			par_change = sqrt(sum(pi_per_iteration[[iter]]-start_pi)^2)
		}
		print(c("par change is", par_change))
		if (par_change<1e-50) {
			break
		}
	}
	cat("\n")
	return(pi_per_iteration)
}
