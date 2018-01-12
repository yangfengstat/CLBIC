##### Code Used for Simulation Settings from Table 5 #####

source('main.R')

##### Network Settings ######

rho = 0.2; alpha = 0.8; rho_n = 0.03
comm.size.vec = c(60,90,120,150,60,90)  
K = length(comm.size.vec); N = sum(comm.size.vec)                

totsim = 10
scorecbic.res = rep(0,totsim); comm.scorecbic = NULL; gfcbic = rep(0,totsim); medcbic = rep(0,totsim)
scorebic.res = rep(0,totsim); comm.scorebic = NULL; gfbic = rep(0,totsim); medbic = rep(0,totsim)
error.rate = rep(0,totsim); ratio.oracle = rep(0,totsim); ratio.sc = rep(0,totsim)

for (p in 1:totsim){
    set.seed(p)
    C = NULL
    for(k in 1:K){ C = c(C,rep(k, comm.size.vec[k])) }
    theta = runif(N, min = 1/5, max = 9/5) 
    A.aux = Global.Correlated.ADCBM(comm.size.vec, alpha, rho_n, GC = "constant", rho, flag.theta = TRUE, value.theta = theta)
    A = A.aux$A   
    Omega = A.aux$Omega          

    RN = graph.adjacency(A, mode="undirected")
    conn.comp = clusters(RN)
    aux = which(conn.comp$csize == max(conn.comp$csize))
    maximal.conn.comp = which(conn.comp$membership == aux)  
    A = A[maximal.conn.comp,maximal.conn.comp]
    Omega = Omega[maximal.conn.comp,maximal.conn.comp]
    removed.labels = setdiff(1:N,maximal.conn.comp)
    if(length(removed.labels) != 0) C = C[-removed.labels]
    RN = graph.adjacency(A, mode="undirected")

    r1 = sb.BIC(A, 1:18, model="DCBM", DetectionAlg = "SCORE", composite = T); scorecbic.res[p] = which(r1$obj.fun==min(r1$obj.fun))
    comm.scorecbic = SCORE(A, which(r1$obj.fun==min(r1$obj.fun)))
    r2 = sb.BIC(A, 1:18, model="DCBM", DetectionAlg = "SCORE", composite = F); scorebic.res[p] = which(r2$obj.fun==min(r2$obj.fun))
    comm.scorebic = SCORE(A, which(r2$obj.fun==min(r2$obj.fun)))

    fit = SCORE(A, K)
    label = fit$Membership.Vector
    error.rate[p] = Hamming.Error(label, C, K)$HE/dim(A)[1]

    #SCORE Estimate#
    label.SC=SCORE(A, K)$Membership.Vector
    label.pop=Hamming.Error(label.SC, C, K)$new.label
    aux.label = unique(label.pop)

    for(i in 1:length(C)) label.SC[i] = which(aux.label == label.SC[i])
    
    est.theta = estThetaDCBM(A, label.SC, K)
    est.omega = estOmega(A, label.SC, K)
    
    Omega.SC = matrix(0, dim(Omega)[1], dim(Omega)[1])
    for(i in 1:dim(Omega)[1]){
        for(j in i:dim(Omega)[1]){
            Omega.SC[i,j] = est.omega[i]*est.omega[j]*est.theta[label.SC[i],label.SC[j]]
            Omega.SC[j,i] = Omega.SC[i,j]
        }
    } 

    #Oracle Estimate#
    est.theta = estThetaDCBM(A, C, K)
    est.omega = estOmega(A, C, K)
    
    Omega.Oracle = matrix(0, dim(Omega)[1], dim(Omega)[1])
    for(i in 1:dim(Omega)[1]){
        for(j in i:dim(Omega)[1]){
            Omega.Oracle[i,j] = est.omega[i]*est.omega[j]*est.theta[C[i],C[j]]
            Omega.Oracle[j,i] = Omega.Oracle[i,j]
        }
    } 
 
    diag(Omega) = diag(Omega.SC) = diag(Omega.Oracle) = rep(0, dim(Omega)[1])
    ratio.oracle[p] = norm(Omega.Oracle - Omega, type="F")/norm(Omega, type="F")
    ratio.sc[p] = norm(Omega.SC - Omega, type="F")/norm(Omega, type="F")

    aux = Goodness(A=A, v0=C, v1=comm.scorecbic, v2=comm.scorebic)
    gfcbic[p] = aux$gfcbic
    gfbic[p] = aux$gfbic
    medcbic[p] = aux$ratiocbicmedian
    medbic[p] = aux$ratiobicmedian
 
}       

output = cbind(error.rate,ratio.oracle,ratio.sc,scorecbic.res,scorebic.res,gfcbic,medcbic,gfbic,medbic)
output
