##### Loading Required Packages #####

library(igraph)
library(stringr)
library(Matrix)
library(MASS)
library(mixer)
library(combinat)

comTheta = function(A, label, K){
#### Computes the maximum number of possible edges between communities for a given label assignment ####
    temp = matrix(0, K, K)
    n.com = c()
    for(i in 1: K){
        n.com[i] = sum(label == i)
    }
    for(i in 1: K){
        for(j in i: K){
            if(j == i){
                if(n.com[i] >= 2) temp[i, j] = choose(n.com[i], 2)
            } 
            if(j > i){
                temp[i, j] = n.com[i] * n.com[j]
            }
            temp[j, i] = temp[i, j]
        }
    }
    return(temp)
}

estTheta = function(A, label, K){
#### Estimates the parameter theta under the regular Stochastic Blockmodel ####
    n.A = dim(A)[1]
    com = comTheta(A = A, label = label, K = K)
    temp = matrix(0, K, K)
    for (i in 1:K){
        for (j in i:K){
            ind1 = which(label == i)
            ind2 = which(label == j)
            tempA = A
            diag(tempA) = rep(0, n.A)
            if(com[i, j] + com[j, i] == 0) temp[i, j] == 0
            if(com[i, j] + com[j, i] != 0){
               if(i ==j) 
                  temp[i, j] = sum(tempA[ind1, ind2])/(com[i, j] + com[j, i])   
               else
                  temp[i, j] = sum(tempA[ind1, ind2])/(com[i, j])   
            }
            temp[j, i] = temp[i, j] 
        }
    }
    return(temp)
}

estThetaDCBM = function(A, label, K){
#### Estimates the parameter theta under the Degree-Corrected Blockmodel ####
    n.A = dim(A)[1]
    temp = matrix(0, K, K)
    for (i in 1:K){
        for (j in i:K){
            ind1 = which(label == i)
            ind2 = which(label == j)
            tempA = A
            diag(tempA) = rep(0, n.A)
            temp[i, j] = sum(tempA[ind1, ind2])
            temp[j, i] = temp[i, j] 
        }
    }
    return(temp)
}

estOmega = function(A, label, K){
#### Estimates the parameter lowercase omega under the Degree-Corrected Blockmodel ####
    n.A = dim(A)[1]
    deg.seq = rowSums(A)
    tot.deg.comm = rep(0, K)
    temp = rep(0, n.A)
    for (k in 1:K){
         ind = which(label == k)
         tot.deg.comm[k] = sum(deg.seq[ind])
    }
    for (i in 1:n.A) temp[i] = deg.seq[i]/tot.deg.comm[label[i]]
    return(temp)
}

varThetajack = function(A, label, model = c("SBM","DCBM"), est, K){
#### Computation of the Jackknife covariance matrix estimator ####
    n.A = dim(A)[1]
    temp = matrix(0, 0.5*K*(K+1), 0.5*K*(K+1))
    id = matrix(upper.tri(est, diag=TRUE), K, K, byrow=T)
    est = est[id] 
    
    for (i in 1:n.A){
         A.minus.i = A[-i,-i]
         l.minus.i = label[-i]
         if (model == "SBM"){
             estTilde = estTheta(A.minus.i, l.minus.i, K)
             estTilde = estTilde[id] 
         }
         if (model == "DCBM"){
             estTilde = estThetaDCBM(A.minus.i, l.minus.i, K)
             estTilde = estTilde[id] 
         }         
         temp  = temp + outer(estTilde - est, estTilde - est)
    }
    temp = ((n.A - 1)/(n.A))*temp
    return(temp)
}

lik = function(A, R, K, model = c("SBM","DCBM"), DetectionAlg = c("spectral","SCORE","population"), normalized=TRUE, magnitude=TRUE, true.label=NULL){
#### Composite Likelihood calculations under Stochastic Blockmodels and Degree-Corrected Blockmodels ####

    #### the following lines are to get the community labels
    if (model=="SBM" && DetectionAlg=="spectral"){
       fit = spectral.clustering.fast(R = R, K = K)
       label = fit$Membership.Vector
    }    
    if (model=="DCBM" && DetectionAlg=="SCORE"){
       fit = SCORE.fast(R = R, K = K)
       label = fit$Membership.Vector
    }     
    if (DetectionAlg=="population"){
       label = true.label
    }    
    
    n = dim(A)[1]
    diag(A) = rep(0, n)
    val = 0
    
    if (model == "SBM"){
        est.theta = estTheta(A = A, label = label, K = K)
        for(i in 1: (n-1)){
            for(j in (i+1): n){
                l1 = label[i]
                l2 = label[j]
                est.temp = est.theta[l1, l2]
                if(est.temp < 0.000001) temp1 = (1 - A[i, j]) * log(1 - est.temp)
                if(est.temp > 0.999999) temp1 = A[i, j] * log(est.temp)
                if(est.temp >= 0.000001 && est.temp <= 0.999999) temp1 = A[i, j] * log(est.temp) + (1 - A[i, j]) * log(1 - est.temp)
                val = val + temp1
            }
        }
    }
    if (model == "DCBM"){    
        est.theta = estThetaDCBM(A = A, label = label, K = K)
        est.omega = estOmega(A = A, label = label, K = K)
        deg.seq = rowSums(A)
            s1 = 2*sum(deg.seq*log(est.omega))
            s2 = 0
            for (a in 1:K){
                 for (b in 1:K){
                      if(est.theta[a,b] > 0.000001) s2 = s2 + est.theta[a,b]*log(est.theta[a,b]) - est.theta[a,b]
                 }
            }
            val = s1 + s2  
    }    
    return(list(val = val, label = label, est = est.theta))
}

sb.BIC = function(A, Kseq, model = c("SBM","DCBM"), DetectionAlg = c("spectral","SCORE","population"), normalized=TRUE, magnitude=TRUE, composite = FALSE, true.label=NULL){
#### Model Selection through the CL-BIC criterion for a sequence of candidate values of K ####
    model = match.arg(model)
    DetectionAlg = match.arg(DetectionAlg)
    
    len.Kseq = length(Kseq)
    A.aux = as(A, "dgCMatrix")
    n.A = dim(A)[1]
    if (model == "SBM"){
        # Compute Laplacian Matrix
        if (normalized==TRUE){
            D = as(diag(1/sqrt(rowSums(A.aux))), "dgCMatrix")
            L =  D%*% A.aux %*%D
        }
        if (normalized==FALSE){
            D = as(diag(rowSums(A.aux)), "dgCMatrix")
            L =  D - A.aux
        }
        # Calculate Eigenvectors
        EV = eigen(L, symmetric=TRUE)
        if (magnitude==TRUE){
            EV.sort = sort(abs(EV$values), decreasing=TRUE, index.return=TRUE)
        }
        if (magnitude==FALSE){
            EV.sort = sort(EV$values, decreasing=FALSE, index.return=TRUE)
        }   
    }
    if (model == "DCBM"){
        # Calculate Eigenvectors
        EV = eigen(A.aux, symmetric=TRUE)
        EV.sort = sort(abs(EV$values), decreasing=TRUE, index.return=TRUE)
    }    
    # Create Matrix R
    R = EV$vectors[,EV.sort$ix[1:max(Kseq)]]    
    obj.fun = c()
    mod.dim = c()
    for(i in 1:len.Kseq){
        if(composite == FALSE){
            if(Kseq[i] == 1) d = 1
            if(Kseq[i] > 1) d = choose(Kseq[i], 2) + Kseq[i]
            # Likelihood computations for each candidate value of K
            l = lik(A, R, Kseq[i], model, DetectionAlg, normalized, magnitude, true.label)
            if(Kseq[i] == 1) obj.fun[i] = -2*l$val + d*log(n.A * (n.A - 1)/2)
            if(Kseq[i] > 1) obj.fun[i] = -2*l$val + d*log(n.A * (n.A - 1)/2) 
        }
        if(composite == TRUE){
            # Composite Likelihood computations for each candidate value of K
            l = lik(A, R, Kseq[i], model, DetectionAlg, normalized, magnitude, true.label)
            label = l$label
            est = l$est
            Thetajack = varThetajack(A, label, model, est, Kseq[i])
            # Calculation of model complexity term for CL-BIC
            d = c.d(A, label, model, est, Thetajack)
            if(Kseq[i] == 1) obj.fun[i] = -2*l$val + d*log(n.A * (n.A - 1)/2)
            if(Kseq[i] > 1) obj.fun[i] = -2*l$val + d*log(n.A * (n.A - 1)/2) 
        }
        mod.dim[i] = d
    }
    return(list(obj.fun=obj.fun, mod.dim=mod.dim))
}

c.d = function(A, label, model = c("SBM","DCBM"), est, Thetajack){
#### Computes model complexity term for CL-BIC under Stochastic Blockmodels and Degree-Corrected Blockmodels ####
    K = max(label)
    temp.aux = 0
    n = dim(A)[1]
    diag(A) = rep(0, n)
    
    VJack = matrix(0,K,K)
    id  = matrix(upper.tri(VJack, diag=TRUE), K, K, byrow=T)
    VJack[id] = diag(Thetajack)
    VJack[upper.tri(VJack)] = t(VJack)[upper.tri(VJack)]
    
    if(model == "SBM"){
       for(i in 1: K){
           for(j in i: K){
               correction = 0 
               ind1 = which(label == i)
               ind2 = which(label == j)
               if(est[i,j] > 0.000001){
                  if(i == j)
                     temp11 = sum(A[ind1, ind2])/(2*est[i,j]^2)
                  else
                     temp11 = sum(A[ind1, ind2])/est[i,j]^2           
               }
               else{
                  temp11 = 0
               }
               if(est[i,j] <= 0.999999){
                  if(i == j){
                     correction = length(ind1)
                     temp12 = (sum(1 - A[ind1, ind2]) - correction)/(2*(1 - est[i,j])^2)   
                  }
                  else{
                     temp12 = sum(1 - A[ind1, ind2])/(1 - est[i,j])^2          
                  }
               }
               else{
                  temp12 = 0
               }                                       
               temp1 = temp11 + temp12            
               temp.aux  = temp.aux + VJack[i,j]*temp1
           }
       }
    }
    if(model == "DCBM"){
       for(i in 1: K){
           for(j in i: K){
               if(est[i,j] > 0.000001) temp.aux = temp.aux + VJack[i,j]/est[i,j]
           }
       }        
    
    }
    return(temp.aux)
}

Re.Index.Labels <- function(z,K){
class.index = 1:K; new.z = rep(0,length(z)) #node.index = z;
for (j in 1:K){
    p = min(which(z %in% class.index))
    aux = which(z==z[p]) 
    new.z[aux] = j
    out = z[p]
    class.index = setdiff(class.index,out)   
}
return(new.z)
}

spectral.clustering.fast <- function(R, K){
# Create Matrix R for K-means
R = R[,1:K]
# Classify using K-means
k.means = kmeans(R, centers = K, iter.max = 150, nstart = 50)
# Preparing Output
Membership.Vector = Re.Index.Labels(k.means$cluster, K)
comm = vector("list",K)
for (j in 1:K) comm[[j]]= which(Membership.Vector == j)
return(list(Membership.Vector=Membership.Vector, Clusters=comm))
}

spectral.clustering <- function(A, K, normalized=TRUE, magnitude=TRUE){
# creates normalized graph Laplacian
if (normalized==TRUE){
    A = as(A, "dgCMatrix")
    D = as(diag(1/sqrt(rowSums(A))), "dgCMatrix")
    L =  D%*%A%*%D
}
# creates un-normalized graph Laplacian
if (normalized==FALSE){
    A = as(A, "dgCMatrix")
    D = as(diag(rowSums(A)), "dgCMatrix")
    L =  D - A
}

# Step 1
EV = eigen(L,symmetric=TRUE)
if (magnitude==TRUE){
    EV.sort = sort(abs(EV$values), decreasing=TRUE, index.return=TRUE)
}
if (magnitude==FALSE){
    EV.sort = sort(EV$values, decreasing=FALSE, index.return=TRUE)
}    
X = EV$vectors[,EV.sort$ix[1:K]]

# Step 2
k.means = kmeans(X, centers = K, iter.max = 150, nstart = 50)

# Preparing Output
Membership.Vector = Re.Index.Labels(k.means$cluster, K)
comm = vector("list",K)
for (j in 1:K) comm[[j]]= which(Membership.Vector == j)
return(list(Membership.Vector=Membership.Vector, Clusters=comm))
}


SCORE.fast <- function(R, K){
# Create Matrix R for K-means
R = R[,1:K]
if (K > 1){
    N = dim(R)[1]
    R = R/R[,1]
    R = R[,-1]
}
else{
    N = length(R)
    R = rep(1, N)
}
R[R >= log(N)] = log(N)
R[R <= -log(N)] = -log(N)
# Classify using K-means
k.means = kmeans(R, centers = K, iter.max = 150, nstart = 50)
# Preparing Output
Membership.Vector = Re.Index.Labels(k.means$cluster, K)
comm = vector("list",K)
for (j in 1:K) comm[[j]]= which(Membership.Vector == j)
return(list(Membership.Vector=Membership.Vector, Clusters=comm))
}

SCORE <- function(A, K){
# creates sparse version of adjacency matrix A
A = as(A, "dgCMatrix")
N = dim(A)[1]
# Calculate Eigenvectors
EV = eigen(A, symmetric=TRUE)
EV.sort = sort(abs(EV$values), decreasing=TRUE, index.return=TRUE)
# Create Matrix R
R = EV$vectors[,EV.sort$ix[1:K]]
if (K > 1){
    R = R/R[,1]
    R = R[,-1]
}
else{
    R = rep(1, N)
}
R[R >= log(N)] = log(N)
R[R <= -log(N)] = -log(N)
# Classify using K-means
k.means = kmeans(R, centers = K, iter.max = 150, nstart = 50)
# Preparing Output
Membership.Vector = Re.Index.Labels(k.means$cluster, K)
comm = vector("list",K)
for (j in 1:K) comm[[j]]= which(Membership.Vector == j)
return(list(Membership.Vector=Membership.Vector, Clusters=comm))
}

perm.label.vec <- function(perm, C){
result = c()
for(i in 1:length(C)) result[i] = perm[C[i]]
return(result)
}

Hamming.Error <- function(label, C, K){
#### Computes number of misclustered nodes after permutation of the labels ####
cand.labels = permn(1:K); obj = c(); candidate = vector("list",length(cand.labels))   
for(i in 1:length(cand.labels)){
    candidate[[i]] = perm.label.vec(cand.labels[[i]], C)
    obj[i] = length(which(label != candidate[[i]]))
}
ind = which.min(obj)
winner = candidate[[ind]]
HE = obj[ind]
return(list(HE=HE, new.label=winner))
}

Goodness <- function(A, v0, v1, v2){
#### Computes goodness-of-fit measures ####
n = length(v0); K1 = max(v1$Membership.Vector); K2 = max(v2$Membership.Vector)
a1 = a2 = 0

for (i in 1:(n-1))
for (j in (i+1):n){
     a1 = a1 + as.numeric((v0[i]==v0[j])==(v1$Membership.Vector[i]==v1$Membership.Vector[j]))
     a2 = a2 + as.numeric((v0[i]==v0[j])==(v2$Membership.Vector[i]==v2$Membership.Vector[j]))
}

p1 = a1/(n*(n-1)/2)
p2 = a2/(n*(n-1)/2)

cocbic = matrix(0,K1,K1)
for (i in 1:K1)
for (j in 1:K1){
     s1 = v1$Clusters[[i]]
     s2 = v1$Clusters[[j]]
     cocbic[i,j] = sum(A[s1,s2])
}

cobic = matrix(0,K2,K2)
for (i in 1:K2)
for (j in 1:K2){
     s1 = v2$Clusters[[i]]
     s2 = v2$Clusters[[j]]
     cobic[i,j] = sum(A[s1,s2])
}

vcocbic = as.vector(cocbic[upper.tri(cocbic)])
vcobic = as.vector(cobic[upper.tri(cobic)])

ratiocbicmedian = median(diag(cocbic))/median(vcocbic)
ratiobicmedian = median(diag(cobic))/median(vcobic)

return(list(gfcbic=p1,gfbic=p2,ratiocbicmedian=ratiocbicmedian,ratiobicmedian=ratiobicmedian))
}

Blockwise.Correlated.A <- function(comm.size.vec, Theta, WC = TRUE, WCC = "decaying", rho, BC = FALSE, BCC = "decaying", rhob){
#### Generates correlated adjacency matrices according to the settings described in Tables 1--3 ####
# Fix parameters
K = length(comm.size.vec); N = sum(comm.size.vec)
A = matrix(0, ncol = N, nrow = N); C = NULL
for(i in 1:K){ C = c(C,rep(i, comm.size.vec[i])) }

# Generate Independent Bernoulli Data
for(i in 1:N){
    for(j in i:N){
        A[i,j] = rbinom(1, 1, Theta[C[i],C[j]])
        A[j,i] = A[i,j]
    }
}  

if(WC == TRUE && WCC == "constant"){   
    # form mean vector and covariance matrix for CONSTANT CORRELATION diagonal blocks (row-wise generation of the data)
    for(k in 1:K){
        n = comm.size.vec[k]; s = sum(comm.size.vec[1:k]) - comm.size.vec[k] + 1
        for(i in 1:(n-1)){        
            fr = c(1,rep(rho,n-i-1)); sigma = toeplitz(fr); cholmat = chol(sigma)
            mu = rep(qnorm(Theta[k,k]), n-i)
            z = rnorm(n-i, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
            j = i + s - 1
            A[j,(j+1):sum(comm.size.vec[1:k])] = z
        }
        A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])][lower.tri(A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])])] = t(A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])])[lower.tri(t(A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])]))]
    }     
}

if(WC == TRUE && WCC == "decaying"){   
    # form mean vector and covariance matrix for DECAYING CORRELATION in diagonal blocks
    for(k in 1:K){
        n = comm.size.vec[k]; s = sum(comm.size.vec[1:k]) - comm.size.vec[k] + 1
        for(i in 1:(n-1)){
            fr = sapply(0:(n-i-1),function(x) rho^x); sigma = toeplitz(fr); cholmat = chol(sigma)
            mu = rep(qnorm(Theta[k,k]), n-i)
            z = rnorm(n-i, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
            j = i + s - 1
            A[j,(j+1):sum(comm.size.vec[1:k])] = z
        }
        A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])][lower.tri(A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])])] = t(A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])])[lower.tri(t(A[s:sum(comm.size.vec[1:k]),s:sum(comm.size.vec[1:k])]))]
    }     
}

if(BC == TRUE && BCC == "constant"){   
    # form mean vector and covariance matrix for CONSTANT CORRELATION in off-diagonal blocks (row-wise generation of the data)
    for(t1 in 1:(K-1)){
        for(t2 in (t1+1):K){
            nr = comm.size.vec[t1]; nc = comm.size.vec[t2]
            sr = sum(comm.size.vec[1:t1]) - comm.size.vec[t1] + 1; sc = sum(comm.size.vec[1:t2]) - comm.size.vec[t2] + 1
            fr = c(1,rep(rhob,nc - 1)); sigma = toeplitz(fr); cholmat = chol(sigma)
            mu = rep(qnorm(Theta[t1,t2]), nc)
            for(i in 1:nr){            
                z = rnorm(nc, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
                j = i + sr -1
                A[j,sc:sum(comm.size.vec[1:t2])] = z
            }
            A[sc:sum(comm.size.vec[1:t2]),sr:sum(comm.size.vec[1:t1])] = t(A[sr:sum(comm.size.vec[1:t1]),sc:sum(comm.size.vec[1:t2])])
        }
    }       
}

if(BC == TRUE && BCC == "decaying"){   
    # form mean vector and covariance matrix for DECAYING CORRELATION in off-diagonal blocks
    for(t1 in 1:(K-1)){
        for(t2 in (t1+1):K){
            nr = comm.size.vec[t1]; nc = comm.size.vec[t2]
            sr = sum(comm.size.vec[1:t1]) - comm.size.vec[t1] + 1; sc = sum(comm.size.vec[1:t2]) - comm.size.vec[t2] + 1
            fr = sapply(0:(nc-1),function(x) rhob^x); sigma = toeplitz(fr); cholmat = chol(sigma)
            mu = rep(qnorm(Theta[t1,t2]), nc)
            for(i in 1:nr){
                z = rnorm(nc, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
                j = i + sr -1
                A[j,sc:sum(comm.size.vec[1:t2])] = z
            }
            A[sc:sum(comm.size.vec[1:t2]),sr:sum(comm.size.vec[1:t1])] = t(A[sr:sum(comm.size.vec[1:t1]),sc:sum(comm.size.vec[1:t2])])
        }    
    }     
}
diag(A) = rep(0, N)
return(A=A)
}

Global.Correlated.A <- function(comm.size.vec, Theta, GC = "constant", rho){
#### Generates correlated adjacency matrices according to the settings described in Tables 1--3 ####
# Fix parameters
K = length(comm.size.vec); N = sum(comm.size.vec)
A = matrix(0, ncol = N, nrow = N); C = NULL
for(i in 1:K){ C = c(C,rep(i, comm.size.vec[i])) }

# Generate Independent Bernoulli Data
for(i in 1:N){
    for(j in i:N){  
        A[i,j] = rbinom(1, 1, Theta[C[i],C[j]])
        A[j,i] = A[i,j]
    }
}  

if(GC == "constant"){   
    # form mean vector and covariance matrix for CONSTANT CORRELATION globally
    p = 1
    for(i in 1:(N-1)){
        if(i >= sum(comm.size.vec[1:p])){ p = p + 1 }
        fr = c(1,rep(rho, N-i-1)); sigma = toeplitz(fr); cholmat = chol(sigma); mu = NULL
        for(j in 1:(K-p+1)){
            if(p > 1){
               if(i == sum(comm.size.vec[1:(p-1)])){
                  if(j == 1){ mu = c(mu, rep(qnorm(Theta[p-1,p]), sum(comm.size.vec[1:p])-i)) }
                  else{ mu = c(mu, rep(qnorm(Theta[p-1,j+1]),comm.size.vec[j+p-1])) }               
               }
               else{
                  if(j == 1){ mu = c(mu, rep(qnorm(Theta[p,p]), sum(comm.size.vec[1:p])-i)) }
                  else{ mu = c(mu, rep(qnorm(Theta[p,j+p-1]),comm.size.vec[j+p-1])) }
               }
            }
            else{           
               if(j == 1){ mu = c(mu, rep(qnorm(Theta[p,p]), sum(comm.size.vec[1:p])-i)) }
               else{ mu = c(mu, rep(qnorm(Theta[p,j+p-1]),comm.size.vec[j+p-1])) }            
            }
        }
        z = rnorm(N-i, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
        A[i,(i+1):N] = z
        A[(i+1):N,i] = t(A[i,(i+1):N])
    } 
}

if(GC == "decaying"){   
    # form mean vector and covariance matrix for DECAYING CORRELATION globally
    p = 1
    for(i in 1:(N-1)){
        if(i >= sum(comm.size.vec[1:p])){ p = p + 1 }
        fr = sapply(0:(N-i-1),function(x) rho^x); sigma = toeplitz(fr); cholmat = chol(sigma); mu = NULL
        for(j in 1:(K-p+1)){
            if(p > 1){
               if(i == sum(comm.size.vec[1:(p-1)])){
                  if(j == 1){ mu = c(mu, rep(qnorm(Theta[p-1,p]), sum(comm.size.vec[1:p])-i)) }
                  else{ mu = c(mu, rep(qnorm(Theta[p-1,j+1]),comm.size.vec[j+p-1])) }               
               }
               else{
                  if(j == 1){ mu = c(mu, rep(qnorm(Theta[p,p]), sum(comm.size.vec[1:p])-i)) }
                  else{ mu = c(mu, rep(qnorm(Theta[p,j+p-1]),comm.size.vec[j+p-1])) }
               }
            }
            else{           
               if(j == 1){ mu = c(mu, rep(qnorm(Theta[p,p]), sum(comm.size.vec[1:p])-i)) }
               else{ mu = c(mu, rep(qnorm(Theta[p,j+p-1]),comm.size.vec[j+p-1])) }            
            }
        }
        z = rnorm(N-i, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
        A[i,(i+1):N] = z
        A[(i+1):N,i] = t(A[i,(i+1):N])
    } 
}
diag(A) = rep(0, N)
return(A=A)
}


Global.Correlated.ADCBM <- function(comm.size.vec, alpha, rho_n, GC = "constant", rho, flag.theta = FALSE, value.theta = NULL){
#### Generates correlated adjacency matrices according to the settings described in Tables 4--5 ####
# Fix parameters
K = length(comm.size.vec); N = sum(comm.size.vec); theta = rep(0, N)
A = matrix(0, ncol = N, nrow = N); C = NULL
for(k in 1:K){ C = c(C,rep(k, comm.size.vec[k])) }

# Generate Latent Variable Vector of Theta Parameters
nu = runif(N, 0, 2); prob.alpha  = runif(N, 0, 1)
cum = cumsum(c(alpha, 0.5*(1-alpha), 0.5*(1-alpha)))
Dist = cbind(nu, rep(2/11,N), rep(20/11,N))
for(i in 1:N){
    aux = min(which(prob.alpha[i] <= cum))
    theta[i] = Dist[i,aux]
}
if (flag.theta == TRUE) theta = value.theta

# Form symmetric matrix P and Expected Degree Parameter rho_n
P = matrix(1, K, K) + diag(rep(6, K))
if(rho_n > 0.0357) rho_n = 0.03
P = rho_n*P

CapTheta = matrix(0, N, N)
# Generate Independent Bernoulli Data
for(i in 1:N){
    for(j in i:N){
        CapTheta[i,j] = theta[i]*theta[j]*P[C[i],C[j]]
        CapTheta[j,i] = CapTheta[i,j]
        A[i,j] = rbinom(1, 1, theta[i]*theta[j]*P[C[i],C[j]])
        A[j,i] = A[i,j]
    }
} 

deg.in = NULL; deg.out = NULL
for(k1 in 1:K){
    n1 = comm.size.vec[k1]; s1 = sum(comm.size.vec[1:k1]) - comm.size.vec[k1] + 1
    for(k2 in k1:K){
        n2 = comm.size.vec[k2]; s2 = sum(comm.size.vec[1:k2]) - comm.size.vec[k2] + 1
        A.aux = CapTheta[s1:sum(comm.size.vec[1:k1]), s2:sum(comm.size.vec[1:k2])]
        if(k1 == k2) deg.in = c(deg.in, mean(mean(A.aux[upper.tri(A.aux)])))
        else deg.out = c(deg.out, mean(A.aux))
    }
}   

if(GC == "constant"){   
    # form mean vector and covariance matrix for CONSTANT CORRELATION globally
    for(i in 1:N){
        fr = c(1,rep(rho, N-i)); sigma = toeplitz(fr); cholmat = chol(sigma); mu = NULL
        for(j in i:N){ mu = c(mu, qnorm(theta[i]*theta[j]*P[C[i],C[j]])) }
        z = rnorm(N-i+1, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
        A[i,i:N] = z
        A[i:N,i] = t(A[i,i:N])
    } 
}

if(GC == "decaying"){   
    # form mean vector and covariance matrix for DECAYING CORRELATION globally
    for(i in 1:N){
        fr = sapply(0:(N-i),function(x) rho^x); sigma = toeplitz(fr); cholmat = chol(sigma); mu = NULL
        for(j in i:N){ mu = c(mu, qnorm(theta[i]*theta[j]*P[C[i],C[j]])) }        
        z = rnorm(N-i+1, mean=0, sd=1); z = mu + z%*%cholmat; z = as.numeric(z >= 0)
        A[i,i:N] = z
        A[i:N,i] = t(A[i,i:N])
    } 
}
diag(A) = rep(0, N)
return(list(A=A, Omega=CapTheta, deg.in=mean(deg.in), deg.out=mean(deg.out)))
}
