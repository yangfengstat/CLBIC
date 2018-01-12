##### Code Used for Simulation Settings from Table 4 #####

source('main.R')

##### Network Settings ######

rho = 0.2; alpha = 0.8; rho_n = 0.03
comm.size.vec = c(60,90,120,150)
K = length(comm.size.vec)                  

totsim = 5
scorecbic.res = rep(0,totsim)
scorebic.res = rep(0,totsim)

for (p in 1:totsim){
    set.seed(p)
    A.aux = Global.Correlated.ADCBM(comm.size.vec, alpha, rho_n, GC = "constant", rho)         
    A = A.aux$A
    RN = graph.adjacency(A, mode="undirected")
    conn.comp = clusters(RN)
    aux = which(conn.comp$csize == max(conn.comp$csize))
    maximal.conn.comp = which(conn.comp$membership == aux) 
    A = A[maximal.conn.comp,maximal.conn.comp]

    r1 = sb.BIC(A, 1:12, model="DCBM", DetectionAlg = "SCORE", composite = T); scorecbic.res[p] = which(r1$obj.fun==min(r1$obj.fun))
    r2 = sb.BIC(A, 1:12, model="DCBM", DetectionAlg = "SCORE", composite = F); scorebic.res[p] = which(r2$obj.fun==min(r2$obj.fun))
      
}

scorecbic.res
scorebic.res
