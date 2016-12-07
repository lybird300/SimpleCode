install.packages("statnet")
library(statnet)
library(ISwR)
library(MASS)

############Naming rules
Pajek network data: net(tick #)_(execution order #)(replication #).net

Repast outputter: tieStrength/knowledgeLevel/strategyChoice_(execution order #)(replication #).txt

Detailed network data:


##############calculate the average degree distribution
###probably useless
net40.degreeDist <- rep(0, 40)
net120.degreeDist <- rep(0, 40)
net200.degreeDist <- rep(0, 40)
net40.num = 0
net120.num = 0
net200.num = 0

for(i in 1:35){
	for(j in 1:3){
		myNet <- read.paj(paste("net1500_", i, j, ".net", sep=""))
		if(length(myNet[,])>0){
			myNet.size <- network.size(myNet)
			if(myNet.size==40){
				net40.num <- net40.num+1
				net40.degreeDist <- net40.degreeDist + degree(myNet, cmode="indegree")
			}
			if(myNet.size==120){
				net120.num <- net120.num+1
				net120.degreeDist <- net120.degreeDist + degree(myNet, cmode="indegree")
			}
			if(myNet.size==200){
				net200.num <- net200.num+1
				net200.degreeDist <- net200.degreeDist + degree(myNet, cmode="indegree")
			}
		}
	}
}

par(mfrow=c(2,1))
hist(net40.degreeDist/net40.num, xlab="Degree", main="Degree Distribution (N = 40)", freq=F, prob=TRUE)
hist(net120.degreeDist/net40.num, xlab="Degree", main="Degree Distribution (N = 120)", freq=F, prob=TRUE)
hist(net200.degreeDist/net40.num, xlab="Degree", main="Degree Distribution (N = 200)", freq=F, prob=TRUE)

###display degree distribution (as histogram and line graph)
par(mfrow=c(3,3), mex=.6)
#for(i in 1:35){
	i=12
	n=0
	for(name in 1:3){
		net.degreeDist <- vector()
		for(j in 1:3){
			myNet <- read.paj(paste("net", 500*name, "_", i, j, ".net", sep=""))
			if(length(myNet[,])>0){
				if (n==0){
					n <- network.size(myNet)
					net.degreeDist <- rep(0, n)
				}
				net.degreeDist <- net.degreeDist + degree(myNet, cmode="indegree")
			}
		}
		if(n>0){
			net.degreeDist <- net.degreeDist/3
			h <- hist(net.degreeDist, breaks = seq(0,n-1,by=1), main=paste("i=", i, ";T=", 500*name, sep=""))
			#h <- hist(net.degreeDist, main=paste("i=", i, ";T=", 500*name, sep=""))
			lines(h, lty=2, lwd=3)
		}
	}
#}
?paste
?hist
?is.isolate


g1 <- read.paj("netStat/net20302/net_P13.net")
g1 <- read.paj(paste("net", 500*3, "_", 7, 3, ".net", sep=""))
gplot(g1)
plot.network(g1, usearrows = F)
g2 <- read.paj(paste("net", 500, "_", 27, 3, ".net", sep=""))
gplot(g2)
dc = centralization(g2, degree, mode="graph")

##################################DOE experiment results(network structure)
result <- result/3
?gden
?geodist
for(i in 1:35){  					###indicate the execution order of experiments
	for(j in 1:3){					###indicate the replication #
		for(name in 1:3){				###indicate the tick number (included in the name of .net files)
			g <- read.paj(paste("net", 500*name, "_", i, j, ".net", sep=""))

			g <- read.paj(paste("net", 1500, "_", 5, 1, ".net", sep=""))
			###length(g1[1,])
			if(length(g[,])>0){		###some of the sna measures below are uncalculable for empty networks (also, some of them cannot be calculated on one-edge network...)
				den = gden(g, mode="digraph") ###even though you define social network as a graph instead of a digraph in Repast; otherwise the result will be greater than 1
				geo = geodist(g)
				###geo = geodist(g, inf.replace=40)
				netDiameter = max(geo$gdist)
				avgPathLength = mean(geo$gdist) ###average geodesic distance
			
				deg = degree(g, cmode="indegree")
				maxDegree = max(deg)
				avgDegree = mean(deg)
				dc = centralization(g, degree, mode="graph")
		
				bet <- betweenness(g, gmode="graph")       # Geographic betweenness
				maxBetween = max(bet)
				avgBetween = mean(bet)
				bc = centralization(g, betweenness, mode="graph")

			
				# computes the clustering coefficient of Watts-Strogatz (local clustering coefficient, the density of ego network)
				total = 0
				for (k in 1:length(g[1,])){               ###length(g[1,]) is the number of vectors in the entire network
					# get the local network (including the ego and the alters)
					if(is.isolate(g,k)){
						next;
					}else{	
						neighbors=ego.extract(g,k)
						total = total + gden(neighbors,mode="digraph")
					}
				}
				cc <- total/length(g[1,]) ##average local clustering coefficient
				result <- c(den, netDiameter, avgPathLength, maxDegree, avgDegree, dc, maxBetween, avgBetween, bc, cc)
			
				###write.table(result,file="netStructure.csv",sep=",",row.names=c("den", "netDiameter", "avgPathLength", "maxDegree", "avgDegree", "dc", "maxBetween", "avgBetween", "bc", "cc"))
				write.table(matrix(result, nrow=1, byrow=T),file="netStructure.csv",sep=",",row.names=F, col.names=F, append=T)
			}else{
				result <- rep(0, 10)
				write.table(matrix(result, nrow=1, byrow=T),file="netStructure.csv",sep=",",row.names=F, col.names=F, append=T)
				next
			}
		}
	}
}
?centralization


###test the small-world feature
?rgraph
?abline
par(mfrow=c(1,2))
smallworld.test<-function(agd.obs, lcc.obs, netSize, density, niter=1000){  #Define a simple test function
	agd.rep = rep(0, niter)
	lcc.rep = rep(0, niter)
	g <- rgraph(n=netSize, m=niter, tprob = density, mode="graph")
	#g <- rgraph(n=50, tprob = dens)
	#g[1,,]
	for(i in 1:niter){
		#if(is.connected(g[i,,])){
		#	geo = geodist(g[i,,])
		#	agd.rep[i] = mean(geo$gdist)
		#}else{
		#	agd.rep[i] = netSize
		#}
		geo = geodist(g[i,,], inf.replace=0)
		agd.rep[i] = mean(geo$gdist)
		
		total = 0
		for (k in 1:netSize){               ###length(g[1,]) is the number of vectors in the entire network g
			if(is.isolate(g[i,,],k)){
				next;
			}else{	
				neighbors=ego.extract(g[i,,],k)
				total = total + gden(neighbors,mode="digraph")
			}
		}
		lcc.rep[i] <- total/netSize

	}
  	cat("Observed average geodesic distance: ",agd.obs,"\n")
  	cat("\t\tPr(rho>=obs):",mean(agd.rep>=agd.obs),"\n")
  	cat("\t\tPr(rho<=obs):",mean(agd.rep<=agd.obs),"\n")
  	cat("\t\tPr(|rho|>=|obs|):",mean(abs(agd.rep)>=abs(agd.obs)),"\n")
  	hist(agd.rep,freq=F,main="Histogram of Difference in the Characteristic Path Length")
  	abline(v=agd.obs, lty=2, col="red")

	cat("Observed average cluster coefficient: ",lcc.obs,"\n")
  	cat("\t\tPr(rho>=obs):",mean(lcc.rep>=lcc.obs),"\n")
  	cat("\t\tPr(rho<=obs):",mean(lcc.rep<=lcc.obs),"\n")
  	cat("\t\tPr(|rho|>=|obs|):",mean(abs(lcc.rep)>=abs(lcc.obs)),"\n")
	hist(lcc.rep,freq=F,main="Histogram of Difference in the Local Clustering Coefficient")
  	abline(v=lcc.obs, lty=2, col="red")
}
?ego.extract
###use the data of 105 samples
###change directory if necessary
myData <- read.csv("structureSample.csv", header=T)
detach(myData)
attach(myData)
agd.mean <- rep(0,105)
lcc.mean <- rep(0,105)
for(j in 1:dim(myData)[1]){
		g <- rgraph(n=f1[j], m=100, tprob = density, mode="graph")
		agd.rep = rep(0,100)
		lcc.rep = rep(0,100)
		for(i in 1:100){
			if(is.connected(g[i,,])){
				geo = geodist(g[i,,])
				agd.rep[i] = mean(geo$gdist)
			}else{
				agd.rep[i] = f1[j]
			}

			total = 0
			for (k in 1:f1[j]){               
				if(is.isolate(g[i,,],k)){
					next;
				}else{	
					neighbors=ego.extract(g[i,,],k)
					total = total + gden(neighbors,mode="digraph")##please notice that you made a mistake here by using the digraph mode
				}
			}
			lcc.rep[i] <- total/f1[j]

		}
		agd.mean[j] = mean(agd.rep)
		lcc.mean[j] = mean(lcc.rep)
}
avgPathLength
clust_coeff
?plot
?legend
?points
plot()
?rgraph

plot(avgPathLength, type="l", lty=1, lwd=1, xlim=c(1,105), xlab="Scenario #", ylab="Average geodesic distance")
lines(agd.mean, lty=2, lwd=1)
legend("topright", legend=c("Simulated networks", "Random networks"), lwd=2, lty=c(1,2))

plot(clust_coeff, type="l", lty=1, lwd=1, xlim=c(1,105), xlab="Scenario #", ylab="Local clustering coefficient")
lines(lcc.mean, lty=2, lwd=1)
legend("topright", legend=c("Simulated networks", "Random networks"), lwd=2, lty=c(1,2))

which((agd.mean > avgPathLength) & (lcc.mean/2 < clust_coeff))


write.table(matrix(agd.mean, nrow=1, byrow=T),file="smallworld.csv",sep=",",row.names=F, col.names=F, append=T)
write.table(matrix(lcc.mean, nrow=1, byrow=T),file="smallworld.csv",sep=",",row.names=F, col.names=F, append=T)


###mark out scenario networks which show "small-world" effect: 
###(1) characteristic path length is relatively small but lower bounded by that of a random network
###(2) local clustering coefficient is larger than that of a random network
sampleData <- read.csv("smallworldTest.csv", header=T)##in this data set, lcc.mean has been divided by 2
attach(sampleData)
detach(sampleData)
par(mfrow=c(1,2))

##x=which((Avg_geo_s >= Avg_geo_r) & (clust_coeff_s > clust_coeff_r))
mark=c(5, 15 ,21, 25, 26, 33)

###plot the average geodesic distance of the sameple networks and random networks (same density)
plot(Scenario, Avg_geo_s, type="p", pch=1, ps=2, xaxt="n", xlab="Scenario #", ylim=c(1,250), ylab="Average geodesic distance") ###plot without automatic x-axis
axis(side=1, at=Scenario, labels=Scenario, las=1) ###draw the x axis by myself
points(Avg_geo_r, pch=16, ps=2)
legend("topright", legend=c("Simulated networks", "Random networks"), pch=c(1,16), cex=0.9)
abline(v=mark, lty=2)
###plot the average local clustering coefficient of the sameple networks and random networks (same density)
plot(Scenario, clust_coeff_s, type="p", pch=1, ps=2, xaxt="n", xlab="Scenario #", ylim=c(0,1), ylab="Average local clustering coefficient") ###plot without automatic x-axis
axis(side=1, at=Scenario, labels=Scenario, las=1) ###draw the x axis by myself
points(clust_coeff_r, pch=16, ps=2) # the clustering coefficient of random networks has been timed 2
legend("topright", legend=c("Simulated networks", "Random networks"), pch=c(1,16), cex=0.9)
abline(v=mark, lty=2)

?abline
?lines

?rep

##examine the small-world effect with small-world Q (the larger than 1, the better)
which(sampleData$clust_coeff_s > sampleData$clust_coeff_r)
x = (clust_coeff_s/clust_coeff_r)/(Avg_geo_s/Avg_geo_r) ###x must be >> 1
order(x)

#compare the "small-world" effect in low-cost and high-cost networks in graph
lowCost <- read.csv("summary_lowcost.csv", header=T)
highCost <- read.csv("summary_highcost.csv", header=T)

###the code is used for low-cost network, change corresponding variable names from "lowCost" to "highCost" for high-cost network
attach(highCost)
nTick = dim(highCost)[1]
agd.eachTick = rep(0, nTick)
lcc.eachTick = rep(0, nTick)
for(j in 1:nTick){ 
	g <- rgraph(n=50, m=100, tprob = Den[j], mode="graph")
	##g <- rgraph(n=50, m=1, tprob = 0.33306122, mode="graph")
	##gplot(g)
	##gden(g)
	agd.rep = rep(0,100) ##store the average geodesic distance of m=100 random networks
	lcc.rep = rep(0,100) ##store the average local clustering-coefficient of m=100 random networks
	for(i in 1:100){
		if(is.connected(g[i,,])){
			geo = geodist(g[i,,])
			agd.rep[i] = mean(geo$gdist)
		}else{
			agd.rep[i] = 50 ##the total number of nodes in both low-cost network and high-cost network is 50
		}
		total = 0
		for (k in 1:50){               
			if(is.isolate(g[i,,],k)){ ###the node has no neighbors
				next;
			}else{	
				neighbors=ego.extract(g[i,,],k)
				total = total + gden(neighbors,mode="graph")
			}
		}
		lcc.rep[i] <- total/50
	}
	agd.eachTick[j] = mean(agd.rep)
	lcc.eachTick[j] = mean(lcc.rep)
}
plot(Avg_geod ~ Tick, type="l", lty=1, lwd=1, ylim = c(0,49), xlab="Tick", ylab="Average geodesic distance")
ramNet = cbind(c(Tick), agd.eachTick, lcc.eachTick)
#points(ramNet[,2] ~ ramNet[,1], pch=16)
?points
lines(ramNet[,2] ~ ramNet[,1], lty=2, lwd=1)
legend("topright", legend=c("Simulated networks", "Random networks"), lwd=2, lty=c(1,2))
plot(Loc_CC, type="l", lty=1, lwd=1, ylim = c(0.5, 3), xlab="Tick", ylab="Average local clustering coefficient")
lines(ramNet[,3] ~ ramNet[,1], lty=2, lwd=1)
?line
?dim
?rgraph

###draw the degree distribution of the sample graph
g <- read.paj(paste("net", 500*name, "_", i, j, ".net", sep=""))

###################consider isolates in the calculation of local cc
for(i in 23:35){  					###indicate the execution order of experiments
	for(j in 1:3){					###indicate the replication #
		for(name in 1:3){				###indicate the tick number (included in the name of .net files)
			g <- read.paj(paste("net", 500*name, "_", i, j, ".net", sep=""))
			###g <- read.paj(paste("net", 500*3, "_", 24, 1, ".net", sep=""))
			if(length(g[,])>0){		###some of the sna measures below are uncalculable for empty networks (also, some of them cannot be calculated on one-edge network...)
				# computes the clustering coefficient of Watts-Strogatz (local clustering coefficient, the density of ego network)
				total = 0
				for (k in 1:length(g[1,])){               ###length(g[1,]) is the number of vectors in the entire network
					# get the local network (including the ego and the alters)
					if(is.isolate(g, k)){
						next;
					}else{	
						neighbors=ego.extract(g,k)
						total = total + gden(neighbors,mode="digraph")
					}
				}
				cc <- total/length(g[1,])
				write.table(matrix(cc, nrow=1, byrow=T),file="cc.csv",sep=",",row.names=F, col.names=F, append=T)
			}else{
				write.table(matrix(c(0), nrow=1, byrow=T),file="cc.csv",sep=",",row.names=F, col.names=F, append=T)
				next
			}
		}
	}
}

?geodist
?read.paj
rm(list=ls())
library(network)
library(sna)
###ERGM
myNet <- read.graph("net150_k.net", format="Pajek")
myNet <- read.paj("net0_s.net")
myNet %e% "weight" <- myNet %e% "net150_k1"
delete.edge.attribute(myNet, "net150_k1")
gplot(myNet)
dyad.census(myNet)
list.edge.attributes(myNet)
mat <- matrix(0, nrow=40, ncol=40)
myNet <- network(mat)
plot.sociomatrix(myNet)
par(mfrow=c(1,1))
	myAttribute = read.csv("nodeAttributes_275.csv", header=F)
	myNet = read.paj("net275_s.net")
	temp = get.edge.value(myNet, "net275_s1")
	set.edge.attribute(myNet, "weight", temp)
	#get.edge.value(myNet, "net275_s1")
	delete.edge.attribute(myNet, "net275_s1")
	myNet %v% "vertex.names" <- myAttribute[,2]
	myNet %v% "knowledge" <- myAttribute[,3]
	myNet %v% "revenue" <- myAttribute[,5]
	myNet %v% "strategy" <- myAttribute[,4]
	myNet %v% "change" <- myAttribute[,6]
	k <- myNet %v% "knowledge"
	w <- myNet %e% "weight"
	##plot(myNet,  vertex.col="strategy", displaylabels = F, label.border = F, label = network.vertex.names(myNet1))
	plot(myNet, vertex.cex=k*0.5, vertex.col="strategy", displaylabels = F, label.border = F, edge.lwd=w*30)
	legend("topright",fill=0:3,legend=c("W&L","S&L","W&M","S&M"), bty="n", trace=T)

as.sociomatrix(myNet, "weight")


# Visualization with gplot-------------------------------------------

#use gplot arguments to refine the plot: add labels, scale and color the vertices, etc

gplot(myNet, gmode="graph", label=network.vertex.names(myNet), label.border = F, vertex.col = myNet%v%"knowledge")

#knowledge distribution
plot(myNet1, label=network.vertex.names(myNet1), label.border = F, label.bg=F, vertex.col = myNet1%v%"knowledge")
vals<-sort(unique(round(myNet1 %v% "knowledge")))
legend(x="topleft",fill=2:6, legend = round(vals), bty="n", trace=T)

#revenue distribution
plot(myNet, label=network.vertex.names(myNet), label.border = F, label.bg=F, vertex.col = (-1)* (myNet%v%"revenue"))
vals<-sort(unique(round(myNet %v% "revenue")))
legend(x="topleft",fill=(-1)*vals, legend = vals, bty="n", trace=T)


?round

# for more information, see

?gplot

# Density-------------------------------

gden(myNet)
?gden

#Centrality and centralization----------------------------

deg = degree(myNet,gmode="graph")			#Degree centrality
hist(deg,xlab="Degree",breaks= seq(0,49,by=1), prob=T)
hist(deg,xlab="Degree",freq=F)
hist((myNet%e%"weight"), xlab="Tie strength", prob=T)
w <- myNet %e% "weight"
k <- myNet %v% "knowledge"
mean(k)
max(k)
k[44]
k[27]
deg[27]
deg[44]
mean(deg)
max(deg)
#curve(dnorm(x, mean=mean(x), sd=sd(x)), add=T)
#lines(summary(myNet~degree(0:11)),lty=1,lwd=3)

betweenness(fauxhigh,gmode="graph")
closeness(fauxhigh,gmode="graph")
evcent(fauxhigh,gmode="graph")			#Eigenvector centrality

centralization(fauxhigh, betweenness)
?centralization

#Brokerage

?brokerage

fh.br<- brokerage(fauxhigh,fauxhigh%v%"Race" )
fh.br$raw.nli

#Transitivity and the triad census--------------------------

tc <- triad.census(myNet,mode="graph")
dim(tc)
?triad.census

#Conditional uniform graph tests---------------------------

is.connected(myNet)
ct1<-cug.test(myNet,centralization, mode="graph", cmode="size",FUN.arg=list(FUN=degree)) 
ct2<-cug.test(myNet,centralization,mode="graph", cmode="size",FUN.arg=list(FUN=betweenness)) 
ct3<-cug.test(myNet,gden,mode="graph", cmode="size",FUN.arg=list(mode="graph"))
# Note that we here must pass not only arguments to cug.test, but also to centralization!
ct2                                       # Print the result
plot(ct2) 

?cug.test

#pertumation test-------

?geodist
?gden
geo = geodist(myNet, inf.replace=0)
agd.obs = mean(geo$gdist)
ego.den = rep(0,50)
total = 0
for (k in 1:50){               ###length(g[1,]) is the number of vectors in the entire network
	if(is.isolate(myNet,k)){
		next;
	}else{	
		neighbors=ego.extract(myNet,k, neighborhood = "combined")
		ego.den[k] = gden(neighbors,mode="graph")*2
		total = total + ego.den[k]
	}
}
lcc.obs <- total/50
dens <- gden(myNet, mode="graph")*2

#see above for the codes of this function
smallworld.test(agd.obs, lcc.obs, 50, dens, niter=1000)

#Positions, roles, blockmodels

?equiv.clust
?blockmodel
ec<-equiv.clust(myNet,mode="graph")          # Complete link SE clustering using Hamming dist
ec                              # Gives not so useful summary info
plot(ec)                        # Shows dendrogram

# Now that we have our SE clustering, let's try a blockmodel
bm<-blockmodel(myNet,ec,h=25)               # Cut the tree at 25
summary(bm)
bm                                          # Examine the blockmodel
# One way to plot the block image:
gplot(bm$block.model,diag=T,edge.lwd=bm$block.model*3) 

# Can plot the role structure
gplot(bm$block.model,diag=T,edge.lwd=bm$block.model*3, vertex.cex=table(bm$block.membership)) #whoa, not good! What happened?

table(bm$block.membership)

gplot(bm$block.model,diag=T,edge.lwd=bm$block.model*3, vertex.cex=log(table(bm$block.membership))) #whoa, not good!
																   #when the size distribution of blocks is not even, use log function

#-Using the blockmodel function with fixed membership---------------------------
#
# The "blockmodel" function in sna can be used for various sorts of
# blockmodeling; here, we will pass it a membership vector, and use it to
# extract some basic information.
memb<-match(fauxhigh%v%"Race",sort(unique(fauxhigh%v%"Race"))) # Race vector
bm<-blockmodel(fauxhigh,memb)                                  # Block by race
bm
names(bm)                                              # Examine components
bmr<-bm$block.model                                    # Get the block densities

#Plot the reduced form blockmodel
gplot(bmr,edge.lwd=bmr*200,diag=TRUE,usearrows=FALSE)

#connectivity and substructure

cd <- component.dist(myNet)
sum(cd$cdist) ##the number of components


#Cohesive subgroups--------------------------------------------------


?kcores
?clique.census


cs<-clique.census(myNet,mode="graph", tabulate.by.vertex = F)
cs$clique.count                                     # Count of cliques by vertex
cs$cliques                                          # Enumerate all cliques

# Can also get co-membership information
cs<-clique.census(myNet,clique.comembership="sum", tabulate.by.vertex = F)    # Add up all
cs$clique.comemb
cs<-clique.census(myNet,clique.comembership="bysize", tabulate.by.vertex = F) # Sort by size
cs$clique.comemb[1,,]                                # 1-cliques
cs$clique.comemb[2,,]                                # 2-cliques
cs$clique.comemb[3,,]     


data(fauxhigh)
kc1<-kcores(myNet,mode="graph")                    # Calculate kcores 

#Examine the core distribution
table(kc1)

?gplot

# Visualize the core structure
gplot(myNet,vertex.col=heat.colors(max(kc1)+1)[kc1+1], gmode="graph", label=network.vertex.names(myNet), label.border = F, label.bg=F, )
gplot(myNet[kc1>1,kc1>1],vertex.col=
    heat.colors(max(kc1[kc1>1])+1)[kc1[kc1>1]+1],gmode="graph")        # 2-core
gplot(myNet[kc1>3,kc1>3],vertex.col=
    heat.colors(max(kc1[kc1>3])+1)[kc1[kc1>3]+1],gmode="graph")        # 4-core

?heat.colors


#Homophily----------------------------

mms<-mixingmatrix(myNet,"strategy")               # Compute with a network object and
mmg<-mixingmatrix(fauxhigh,"Grade")             # attribute name  
mmr<-mixingmatrix(fauxhigh,"Race")




#Network regression-----------------------------

gcor(as.sociomatrix(myNet),as.sociomatrix(myNet1), mode="graph")
gcor(as.sociomatrix(yourNet),as.sociomatrix(yourNet1), mode="graph")
qt1<-qaptest(list(as.sociomatrix(myNet),as.sociomatrix(myNet_pre)),gcor,g1=1,g2=2)
summary(qt2)                                       # Examine the results
plot(qt2, xlim=c(-0.1, 0.7))
qt2<-qaptest(list(as.sociomatrix(yourNet),as.sociomatrix(yourNet1)),gcor,g1=1,g2=2)

mod <- netlm(yourNet1, yourNet1_pre, mode="graph")
summary(mod)
?netlm
###compare "black" and "white" network
###before running this part of script, make sure you've changed to the right working directory
library(network)
library(sna)

####collect the organizational-level attributes at each time run
#avg_k = rep(0,21)
#avg_s = rep(0,21)
#avg_r = rep(0,21)
#per_us = rep(0,21)
#per_strategy = matrix(0, nrow=21, ncol=4)
#den = rep(0,21)
#avgDegree = rep(0,21)
#avgBetween = rep(0,21)
#degC = rep(0,21)
#betC = rep(0,21)
#triad = matrix(0, nrow=21, ncol=4)
agd.obs = rep(0,21)
lcc.obs = rep(0,21)
nComp = rep(0,21)

	myAttribute = read.csv("nodeAttributes_1.csv", header=F)
	myNet = read.paj("net1_s.net")
	temp = get.edge.value(myNet, "net1_s1")
	set.edge.attribute(myNet, "weight", temp)
	#get.edge.value(myNet, "net1_s1")
	delete.edge.attribute(myNet, "net1_s1")
	myNet %v% "vertex.names" <- myAttribute[,2]
	myNet %v% "knowledge" <- myAttribute[,3]
	myNet %v% "revenue" <- myAttribute[,5]
	myNet %v% "strategy" <- myAttribute[,4]
	myNet %v% "change" <- myAttribute[,6]
	#k <- myNet %v% "knowledge"
	#w <- myNet %e% "weight"
	##plot(myNet, vertex.cex=k*0.5, vertex.col="strategy", displaylabels = F, label.border = F, label = network.vertex.names(myNet1), edge.lwd=w)
	#plot(myNet, vertex.col="strategy", displaylabels = F, label.border = F, label = network.vertex.names(myNet), edge.lwd=w)
	#vals<-sort(unique(myNet %v% "strategy"))   # Let's add a legend....
	#legend("topleft",fill=1:4,legend=c("Weak&Less","Strong&Less","Weak&More","Strong&More"), bty="n")
	
	#avg_k[1] = mean(myNet %v% "knowledge")
	#avg_s[1] = mean(myNet %e% "weight")
	#avg_r[1] = mean(myNet %v% "revenue")
	#per_us[1] = length(which(myNet %v% "change" == 1))/50
	#for(k in 1:4){
		#per_strategy[1,k] = length(which(myNet %v% "strategy" == k-1))/50
	#}
	#for(k in 1:4){
	#	tc = triad.census(myNet,mode="graph")
	#	triad[1,k] = tc[1,k]
	#}
	geo = geodist(myNet, inf.replace=length(myNet1[1,]))
	agd.obs[1] = mean(geo$gdist)
	ego.den = rep(0,50)
	total = 0
	for (k in 1:50){               ###length(g[1,]) is the number of vectors in the entire network
		if(is.isolate(myNet,k)){
			next;
		}else{	
			neighbors=ego.extract(myNet,k, neighborhood = "combined")
			ego.den[k] = gden(neighbors,mode="graph")*2
			total = total + ego.den[k]
		}
	}
	lcc.obs[1] <- total/50
	cd <- component.dist(myNet)
	nComp[1] = sum(cd$cdist) ##the number of components
	#den[1] = gden(myNet, mode="graph") ###even though you define social network as a graph instead of a digraph in Repast; otherwise the result will be greater than 1
	#deg = degree(myNet, gmode="graph")
	#avgDegree[1] = mean(deg)
	#degC[1] = centralization(myNet, degree, mode="graph")
			
	#bet = betweenness(myNet, gmode="graph")       # Geographic betweenness
	#avgBetween[1] = mean(bet)
	#betC[1] = centralization(myNet, betweenness, mode="graph")


for(i in 1:20){
	i=19
	myAttribute = read.csv(paste("nodeAttributes_", 25*i, ".csv", sep=""), header=F)
	yourNet_475 = read.paj(paste("net", 25*i, "_s.net",sep=""))
	temp = get.edge.value(yourNet_475, paste("net", 25*i, "_s1",sep=""))
	set.edge.attribute(yourNet_475, "weight", temp)
	#list.edges.value(yourNet_475)
	delete.edge.attribute(yourNet_475, paste("net", 25*i, "_s1",sep=""))
	yourNet_475 %v% "vertex.names" <- myAttribute[,2]
	yourNet_475 %v% "knowledge" <- myAttribute[,3]
	yourNet_475 %v% "revenue" <- myAttribute[,5]
	yourNet_475 %v% "strategy" <- myAttribute[,4]
	yourNet_475 %v% "change" <- myAttribute[,6]

	#k <- yourNet1 %v% "knowledge"
	#w <- yourNet1 %e% "weight"
	##plot(myNet, vertex.cex=k, vertex.col="strategy", displaylabels = F, label.border = F, label = network.vertex.names(myNet), edge.lwd=w)
	#plot(yourNet1, vertex.cex=k*0.5, vertex.col="strategy", displaylabels = F, label.border = F, label = network.vertex.names(myNet2), edge.lwd=w*3)
	#vals<-sort(unique(myNet %v% "strategy"))   # Let's add a legend....
	#legend("topleft",fill=1:4,legend=c("WL","SL","WM","SM"), bty="n")

	#avg_k[i+1] = mean(myNet %v% "knowledge")
	#avg_s[i+1] = mean(myNet %e% "weight")
	#avg_r[i+1] = mean(myNet %v% "revenue")
	#per_us[i+1] = length(which(myNet %v% "change" == 1))/50
	#for(k in 1:4){
		#per_strategy[i+1,k] = length(which(myNet %v% "strategy" == k-1))/50
	#}
	#for(k in 1:4){
	#	tc = triad.census(myNet,mode="graph")
	#	triad[i+1,k] = tc[1,k]
	#}
	geo = geodist(myNet, inf.replace=0)
	agd.obs[i+1] = mean(geo$gdist)
	ego.den = rep(0,50)
	total = 0
	for (k in 1:50){               ###length(g[1,]) is the number of vectors in the entire network
		if(is.isolate(myNet,k)){
			next;
		}else{	
			neighbors=ego.extract(myNet,k, neighborhood = "combined")
			ego.den[k] = gden(neighbors,mode="graph")*2
			total = total + ego.den[k]
		}
	}
	lcc.obs[i+1] <- total/50
	cd <- component.dist(myNet)
	nComp[i+1] = sum(cd$cdist) ##the number of components
	#den[i+1] = gden(myNet, mode="digraph") ###even though you define social network as a graph instead of a digraph in Repast; otherwise the result will be greater than 1
	#deg = degree(myNet, gmode="graph")
	#deg1 = degree(myNet4, gmode="graph")
	#avgDegree[i+1] = mean(deg)
	#degC[i+1] = centralization(myNet, degree, mode="graph")
			
	#bet = betweenness(myNet, gmode="graph")       # Geographic betweenness
	#avgBetween[i+1] = mean(bet)
	#betC[i+1] = centralization(myNet, betweenness, mode="graph")
}

write.table(cbind(agd.obs, lcc.obs, nComp), file="triad_lowcost.csv", sep=",", row.names=F, col.names=F)

performance <- cbind(avg_k, avg_s, avg_r, per_us, per_strategy, den, avgDegree, avgBetween, degC, betC)
write.table(performance, file="summary_highcost.csv", sep=",", row.names=F, col.names=F)

#######################################################
library(ergm)
delete.edge.attribute(myNet,"names")
delete.edge.attribute(myNet_pre,"names")
delete.edge.attribute(yourNet,"names")
delete.edge.attribute(yourNet_pre,"names")
gplot(myNet_500)
gplot(yourNet_500)

###myNet--high direct cost, net500; myNet_pre--net475; myNet1--net275; myNet1_pre--net125
model1.1 <- ergm(myNet_300 ~ edges + absdiff("knowledge") + gwdegree(0.4, fixed=T) + triangle, interval=1000) ###this one is good
#model1.2 <- ergm(myNet_500 ~ edges + edgecov(myNet_475, "weight") + nodecov("knowledge") + absdiff("knowledge") + gwdegree(3) + kstar(2), burnin=15000, MCMCsamplesize=15000, interval=1000)
#model1.3 <- ergm(myNet1 ~ edges + edgecov(myNet1_pre, "weight") + nodecov("knowledge") + absdiff("knowledge") + gwdegree(3))
 


###yourNet--low direct cost, net500; yourNet_pre--net475; yourNet1--net275; yourNet1_pre--net125
model2.1 <- ergm(yourNet_325 ~ edges + nodecov("knowledge") + absdiff("knowledge") + gwdegree(0.4, fixed=T) + triangle, interval=1000) ###this one is good
#model2.2 <- ergm(yourNet_500 ~ edges + edgecov(yourNet_475, "weight") + nodecov("knowledge") + absdiff("knowledge") + gwdegree + triangle)
#model2.3 <- ergm(yourNet ~ edges + edgecov(yourNet_pre, "weight") + nodecov("knowledge") + absdiff("knowledge") + gwdegree(3))

mcmc.diagnostics(model2.1)
summary(model1.1)

mgof <- gof(model1.1, GOF=~model + distance + degree + triadcensus, interval = 1000)
#mgof <- gof(model1.2, GOF=~model + distance + triadcensus)
par(mfrow=c(2,2))
plot(mgof)
print(mgof)
gof.model3.2
?gof

?legend
####copy summary_highcost.csv and summary_lowcost.csv to the upper fold and then change the working directory
highCost <- read.csv("summary_highcost.csv", header=T)
lowCost <- read.csv("summary_lowcost.csv", header=T)
par(mfrow=c(1,2), mar=c(4,4,2,4))
plot(highCost$Avg_geod ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Average geodesic distance (isolates excluded)")
lines(lowCost$Avg_geod ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Loc_CC ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Local Clustering Coefficient")
lines(lowCost$Loc_CC ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$nComponent ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Number of Components")
lines(lowCost$nComponent ~ lowCost$Tick, lty=2, lwd=1)
legend("topright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))

plot(highCost$Avg_K ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Average knowledge level", ylim=c(3,6))
lines(lowCost$Avg_K ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))

lowCost$Den[lowCost$Tick]
plot(highCost$Avg_K ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Average knowledge level")
lines(lowCost$Avg_K ~ lowCost$Tick, lty=2, lwd=1)
legend("topright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Avg_S ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Average tie strength")
lines(lowCost$Avg_S ~ lowCost$Tick, lty=2, lwd=1)
legend("topright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Avg_R ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylim = c(-20,0.1),ylab="Average revenue")
lines(lowCost$Avg_R ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Per_US ~ highCost$Tick, type="l", lty=1, lwd=1, ylim=c(0.65,1), xlab="Tick", ylab="Proportion of Unsatisfied Actors")
lines(lowCost$Per_US ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Den ~ highCost$Tick, type="l", lty=1, lwd=1, ylim=c(0, 0.34), xlab="Tick", ylab="Network Density")
lines(lowCost$Den ~ lowCost$Tick, lty=2, lwd=1)
legend("topleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Avg_deg ~ highCost$Tick, type="l", lty=1, lwd=1, ylim=c(0,8), xlab="Tick", ylab="Average degree centrality")
lines(lowCost$Avg_deg ~ lowCost$Tick, lty=2, lwd=1)
legend("topleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$Avg_bet ~ highCost$Tick, type="l", lty=1, lwd=1, ylim=c(0,71), xlab="Tick", ylab="Average betweenness centrality")
lines(lowCost$Avg_bet ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$degC ~ highCost$Tick, type="l", lty=1, lwd=1, xlab="Tick", ylab="Network Degree Centralization", ylim=c(0,0.14))
lines(lowCost$degC ~ lowCost$Tick, lty=2, lwd=1)
legend("bottomright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))
plot(highCost$betC ~ highCost$Tick, type="l", lty=1, lwd=1, ylim=c(0, 0.12), xlab="Tick", ylab="Average Betweenness Centralization")
lines(lowCost$betC ~ lowCost$Tick, lty=2, lwd=1)
legend("topleft", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))

barplot(as.matrix(rbind(lowCost$Per_WL, lowCost$Per_SL, lowCost$Per_WM, lowCost$Per_SM)), beside = F, xlab="Tick", names.arg=highCost$Tick, ylab="Total", main="Strategy Distribution of Simulation I", legend.text=c("Weak & Less", "Strong & Less", "Weak & More", "Strong & More"), args.legend=list(x="topright",ncol=1, bg="white"))
barplot(as.matrix(rbind(highCost$Per_WL, highCost$Per_SL, highCost$Per_WM, highCost$Per_SM)), beside = F, xlab="Tick", names.arg=highCost$Tick, ylab="Total", main="Strategy Distribution of Simulation II", legend.text=c("Weak & Less", "Strong & Less", "Weak & More", "Strong & More"), args.legend=list(x="bottomright", ncol=1, bg="white"))

barplot(as.matrix(rbind(lowCost$Triad0, lowCost$Triad1, lowCost$Triad2, lowCost$Triad3)), beside = F, xlab="Tick", names.arg=lowCost$Tick, ylab="Total", main="Triad census of Simulation I", legend.text=c("0 edge", "1 edge", "2 edges", "3 edges"), args.legend=list(x="bottomright",ncol=1, bg="white"))
barplot(as.matrix(rbind(highCost$Triad0, highCost$Triad1, highCost$Triad2, highCost$Triad3)), beside = F, xlab="Tick", names.arg=highCost$Tick, ylab="Total", main="Triad census of Simulation II", legend.text=c("0 edge", "1 edge", "2 edges", "3 edges"), args.legend=list(x="bottomright",ncol=1, bg="white"))

###degree distribution of two networks (side-by-side comparison)
#To plot a bar graph of two variables on the same graph, you need to begin with a contingency table or crosstabulation
degdist.low <- summary(yourNet~degree(0:49))
degdist.high<- summary(myNet~degree(0:49))
?with
?table
barplot(rbind(degdist.low, degdist.high),main="",xlab="",ylab="",axis.lty=1,legend=T,beside=T)
?degree


###degree distribution comparison with line graphs
plot(summary(yourNet~degree(0:49)), type="l", lty=2, lwd=1,xlab="Degree", ylim=c(0,17), ylab= "Count")
lines(summary(myNet~degree(0:49)), lty=1, lwd=1)
legend("topright", legend=c("Simulation I", "Simulation II"), lwd=2, lty=c(2,1))


?barplot
?par
?legend
?vector
?network
?read.csv

avg_s_k1 <- rep(0,6)
den_k1 <- rep(0,6)
for(i in 6:11){
	myNet = read.paj(paste("net", 50*i, "_k.net",sep=""))
	temp = get.edge.value(myNet, paste("net", 50*i, "_k1",sep=""))
	set.edge.attribute(myNet, "weight", temp)
	#get.edge.value(myNet, "weight")
	delete.edge.attribute(myNet, paste("net", 50*i, "_k1",sep=""))
	avg_s_k1[i-5] = mean(myNet %e% "weight")
	den_k1[i-5] = gden(myNet, mode="digraph")
}

plot(avg_s_k, type="l",lty=1, ylim=c(-0.4,0.1))
lines(avg_s_k1,lty=2)
plot(den_k, type="l",lty=1, ylim=c(0,0.008))
lines(den_k1, lty=2)

