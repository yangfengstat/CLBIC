##### Code Used for Real Data Analysis Examples 4.2.1 and 4.2.2 #####

source('main.R')

########################## Trade Data Largest Network #################################

#Read Edge Data from Trade.csv file
#trade = read.csv("Trade.csv", header = TRUE) 
year = 15

trade = trade[which(trade$t==year),]; trade = trade[order(trade$exporter),]; 
trade = subset(trade,select=-c(number,unknown,t,exporter,importer,log_gdp_exp,log_gdp_imp,log_distance,polity_exp,polity_imp,cc))
colnames(trade) = c("from","to","weight")

# Read Attributes Data
node.attr = read.csv("Trade_Attributes.csv", header = TRUE) 

# generate Adjacency and Weight Matrices
N = max(c(trade$from,trade$to)); W = matrix(0,N,N)
for (i in 1:dim(trade)[1]){ W[trade$from[i],trade$to[i]] = trade$weight[i] }

for(i in 1:N){
    for(j in i:N){
       if(j>i){ 
          W[i,j] = W[i,j] + W[j,i]
          W[j,i] = W[i,j] 
       }
    }
}
A = matrix(as.numeric(W >=quantile(W,0.50)), N, N)

# extract maximum connected component
RN = graph.adjacency(A, mode="undirected")
conn.comp = clusters(RN)
maximal.conn.comp = which(conn.comp$membership==1)

A = A[maximal.conn.comp,maximal.conn.comp]
node.attr = node.attr[maximal.conn.comp,]
diag(A) = 0
RN = graph.adjacency(A, mode="undirected")
V(RN)$names = as.character(node.attr$country)
V(RN)$names

# choose number of communities K 
r1 = sb.BIC(A, 1:12, model = "SBM", DetectionAlg = "spectral", composite = T); cbic.res = which(r1$obj.fun==min(r1$obj.fun))
comm.cbic = spectral.clustering(A, which(r1$obj.fun==min(r1$obj.fun)))
    
r2 = sb.BIC(A, 1:12, model = "SBM", DetectionAlg = "spectral", composite = F); bic.res = which(r2$obj.fun==min(r2$obj.fun))
comm.bic = spectral.clustering(A, which(r2$obj.fun==min(r2$obj.fun)))

print(cbic.res)
print(r1)
print(comm.cbic)

print(bic.res)
print(r2)
print(comm.bic)

# graph results
par(mfrow=c(1,2))
lay1 = layout.auto(RN)
E(RN)$width = 0.1

V(RN)[unlist(comm.cbic$Clusters[1])]$color = "green"
V(RN)[unlist(comm.cbic$Clusters[2])]$color = "SkyBlue2"
V(RN)[unlist(comm.cbic$Clusters[3])]$color = "yellow"

par(ask=TRUE)
plot(RN, layout=lay1, vertex.size=8, vertex.label=V(RN)$names, vertex.frame.color=NA, vertex.label.cex=0.80)

V(RN)[unlist(comm.bic$Clusters[1])]$color = "darkgoldenrod1"
V(RN)[unlist(comm.bic$Clusters[2])]$color = "mediumseagreen"
V(RN)[unlist(comm.bic$Clusters[3])]$color = "green"
V(RN)[unlist(comm.bic$Clusters[4])]$color = "thistle1"
V(RN)[unlist(comm.bic$Clusters[5])]$color = "SkyBlue2"
V(RN)[unlist(comm.bic$Clusters[6])]$color = "salmon"
V(RN)[unlist(comm.bic$Clusters[7])]$color = "magenta"
V(RN)[unlist(comm.bic$Clusters[8])]$color = "yellow"
V(RN)[unlist(comm.bic$Clusters[9])]$color = "cyan"
V(RN)[unlist(comm.bic$Clusters[10])]$color = "red"

par(ask=TRUE)
plot(RN, layout=lay1, vertex.size=8, vertex.label=V(RN)$names, vertex.frame.color=NA, vertex.label.cex=0.80)

####### National Longitudinal Study of Adolescent Health - School #7 ########

#Read Edge Data from School7_Arcs.txt file
#g = scan("School7_Arcs.txt",n=-1); totdat = length(g)/3
from = rep(0,totdat); to = rep(0,totdat); weight = rep(0,totdat)
for(i in 1:totdat){ from[i] = g[3*(i-1) + 1] }
for(i in 1:totdat){ to[i] = g[3*(i-1) + 2] }
for(i in 1:totdat){ weight[i] = g[3*(i-1) + 3] }
N = max(c(from,to))

# Read Attributes Data
node.attr = read.csv("School7_Attributes.csv", header = TRUE) 

# generate Adjacency and Weight Matrices
A = matrix(0,N,N); W = matrix(0,N,N)
for (i in 1:totdat){
    A[from[i],to[i]] = 1
    W[from[i],to[i]] = weight[i]
    A[to[i],from[i]] = A[from[i],to[i]]
    W[to[i],from[i]] = W[from[i],to[i]]
}
diag(A) = 0; diag(W) = 0

# extract maximum connected component
RN = graph.adjacency(A, mode="undirected")
plot(RN, vertex.size=7, vertex.frame.color=NA, vertex.label.cex=0.45)
conn.comp = clusters(RN)
maximal.conn.comp = which(conn.comp$membership==1)
no.grade = which(node.attr$grade == 0)
maximal.conn.comp = setdiff(maximal.conn.comp,no.grade)
A.LC = A[maximal.conn.comp,maximal.conn.comp]
node.attr = node.attr[maximal.conn.comp,]
rownames(node.attr) = NULL
dim(A)
dim(A.LC)
A = A.LC
diag(A) = 0
RN = graph.adjacency(A, mode="undirected")
V(RN)$names = maximal.conn.comp

r1 = sb.BIC(A, 1:15, model = "DCBM", DetectionAlg = "SCORE", composite = T); scorecbic.res = which(r1$obj.fun==min(r1$obj.fun))
comm.scorecbic = SCORE(A, which(r1$obj.fun==min(r1$obj.fun)))
set.seed(3)
r2 = sb.BIC(A, 1:15, model = "DCBM", DetectionAlg = "SCORE", composite = F); scorebic.res = which(r2$obj.fun==min(r2$obj.fun))
comm.scorebic = SCORE(A, which(r2$obj.fun==min(r2$obj.fun)))

print(scorecbic.res)
print(r1)
print(comm.scorecbic)
black = unlist(comm.scorecbic$Clusters[4])

print(scorebic.res)
print(r2)
print(comm.scorebic)

##### Proportion Pairs #####

aux = Goodness(A=A, v0=node.attr$grade, v1=comm.scorecbic, v2=comm.scorebic)
gfcbic = aux$gfcbic
gfbic = aux$gfbic
medcbic = aux$ratiocbicmedian
medbic = aux$ratiobicmedian

inspect(gfcbic)
inspect(gfbic)
inspect(medcbic)
inspect(medbic)

# layouts
lay6 = layout.drl(RN, options=list(simmer.attraction=0))

# graph results by grade
par(mfrow=c(1,3))

V(RN)[unlist(comm.scorecbic$Clusters[1])]$color = "green"
V(RN)[unlist(comm.scorecbic$Clusters[2])]$color = "purple"
V(RN)[unlist(comm.scorecbic$Clusters[3])]$color = "red"
V(RN)[unlist(comm.scorecbic$Clusters[4])]$color = "gray20"
V(RN)[unlist(comm.scorecbic$Clusters[5])]$color = "yellow"
V(RN)[unlist(comm.scorecbic$Clusters[6])]$color = "SkyBlue2"
par(ask=TRUE)
plot(RN, layout=lay6, vertex.size=7, vertex.label=V(RN)$names, vertex.frame.color=NA, vertex.label.cex=0.45)

V(RN)[unlist(comm.scorebic$Clusters[1])]$color = "purple"
V(RN)[unlist(comm.scorebic$Clusters[2])]$color = "green"
V(RN)[unlist(comm.scorebic$Clusters[3])]$color = "cyan"
V(RN)[unlist(comm.scorebic$Clusters[4])]$color = "red"
V(RN)[unlist(comm.scorebic$Clusters[5])]$color = "yellow"
V(RN)[unlist(comm.scorebic$Clusters[6])]$color = "gray20"
V(RN)[unlist(comm.scorebic$Clusters[7])]$color = "magenta"
V(RN)[unlist(comm.scorebic$Clusters[8])]$color = "SkyBlue2"
V(RN)[unlist(comm.scorebic$Clusters[9])]$color = "salmon"
par(ask=TRUE)
plot(RN, layout=lay6, vertex.size=7, vertex.label=V(RN)$names, vertex.frame.color=NA, vertex.label.cex=0.45)

V(RN)[node.attr$grade == 7]$color = "SkyBlue2"
V(RN)[node.attr$grade == 8]$color = "yellow"
V(RN)[node.attr$grade == 9]$color = "green"
V(RN)[node.attr$grade == 10]$color = "purple"
V(RN)[node.attr$grade == 11]$color = "red"
V(RN)[node.attr$grade == 12]$color = "gray20"
par(ask=TRUE)
plot(RN, layout=lay6, vertex.size=7, vertex.label=V(RN)$names, vertex.frame.color=NA, vertex.label.cex=0.45)
