##### Code Used for Simulation Settings from Table 3 #####

source('main.R')

##### Network Settings ######
                  
rho = 0.6
comm.size.vec = c(60,90,120,150)
K = length(comm.size.vec)
Theta = matrix(0.05, K, K) + diag(rep(0.30, K))
Theta[K,1:K] = 0.35
Theta[1:K,K] = 0.35

totsim = 3
cbic.res = rep(0,totsim)
bay.res = rep(0,totsim); comm.bay = vector("list",totsim)

for (p in 1:totsim){
    set.seed(p)
    #Generating correlated adjacency matrix
    A = Blockwise.Correlated.A(comm.size.vec, Theta, WC = TRUE, WCC = "decaying", rho, BC = FALSE)

    r1 = sb.BIC(A, 1:12, model = "SBM", DetectionAlg = "spectral", composite = T); cbic.res[p] = which(r1$obj.fun==min(r1$obj.fun))

    xout.bay = mixer(A,qmin=1,qmax=12,method="bayesian"); cand.bay = NULL
    for(j in 1:12){
       cand.bay = c(cand.bay,xout.bay$output[[j]]$criterion)
    }
    bay.res[p] = which(cand.bay==max(cand.bay)) 
      
}

cbic.res
bay.res
