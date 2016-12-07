##Collect data: lrrk2 carriers, promiscuous controls, etc
lrrk2Carriers <- read.table("lrrk2_sample_ID.txt", header=F)
# 23 carriers

finalIBDClst_PairFoundWithSubNull <- scan()
1: 320 160 330 309 329 328 138 130 285 324 149 335 166 143 306 158 134 163 326 325 1039 738 1053 402 823
# 25 items
setdiff(lrrk2Carriers$V1, finalIBDClst_PairFoundWithSubNull)
# excluded carriers: 396 455 316 319 331
setdiff(finalIBDClst_PairFoundWithSubNull, lrrk2Carriers$V1)
# included non-carriers: 149  143 1039  738 1053  402  823

finalIBDClst_PairFoundWithRawPV0.1 <- scan()
1: 163 160 330 309 329 328 320 306 130 324 285 138 335 143 166 149 158 326 134 1039 1053 738 402
# 23 items; excluded carriers: 396 455 316 319 325 331; included non-carriers: 143  149 1039 1053  738  402

finalIBDClst_completeSearch <- scan()
1: 138 324 131 130 306 335 143 166 149 285 152 160 330 158 309 320 329 163 119 328 605 929 1053 397 736 402 731 622
# 28 items; excluded carriers: 396 455 134 316 319 325 326 331; included non-carriers: 131  143  149  152  119  605  929 1053  397  736  402  731  622

finalIBDClst_PairFoundbyGermline <- scan()
1: 160 328 309 329 330 320 163 166 138 149 130 324 285 306 335 143 158 326 134 1039 1053 738 402
# 23 items; excluded carriers: 396 455 316 319 325 331; included non-carriers: 149  143 1039 1053  738  402

finalIBDClst_PairFoundbyRefinedIBD <- scan()
1: 320 329 328 163 324 330 160 309 326 306 143 130 285 166 138 149 331 158 335 738 580 402 1053 823 688 1035
# 26 items; excluded carriers 396 455 134 316 319 325; included non-carriers: 143  149  738  580  402 1053  823  688 1035

install.packages('gplots', repos='http://cran.us.r-project.org')
library(gplots)
comple <- c("396","455","134","316","319","325","326","331")
subNull <- c("396","455","316","319","331")
rawLsh <- c("396","455","316","319","331","325")
germline <- c("396","455","316","319","331","325")
refineIBD <- c("396","455","134","316","319","325")
excludedKnownCarriersList <- list(COM=comple,SUBNULL=subNull,RAWLSH=rawLsh,GERM=germline,REFINE=refineIBD)
venn(excludedKnownCarriersList)

comple <- c(131,143,149,152,119,605,929,1053,397,736,402,731,622)
subNull <- c(149,143,1039,738,1053,402,823)
rawLsh <- c(143,149,1039,1053,738,402)
germline <- c(143,149,738,580,402,1053,823,688,1035)
refineIBD <- c(143,149,738,580,402,1053,823,688,1035)
includedNonCarriersList <- list(COM=comple,SUBNULL=subNull,RAWLSH=rawLsh,GERM=germline,REFINE=refineIBD)
venn(includedNonCarriersList)
###################################################################################
library(ndtv)
library)tsna)
##read node file
rawVerts<-read.table("SubGraph_node.txt",header=TRUE)
   
##combine edge files
files <- list.files(pattern = "*_edge\\.txt")
tables <- lapply(files, read.table, header = TRUE)
rawEdges <- do.call(rbind , tables)

#We will pass in the vertex data via the vertex.spells argument, but as it is expecting data in the order ‘onset’,‘terminus’,‘id’, we need to reorder the data.frame as we pass it in. 
#Similarly, the edge data must be ordered by ‘onset’,‘terminus’,‘tail’,‘head’ to be processed correctly by the edge.spells argument.
IBDPairFromRawLSH <-networkDynamic(vertex.spells=rawVerts[,c(7,8,1)], edge.spells=rawEdges[,c(3,4,1,2,5)],create.TEAs = TRUE, edge.TEA.names = c('weight'))
list.edge.attributes(IBDPairFromRawLSH)

# add in the unchanging vertex attribute data
set.vertex.attribute(IBDPairFromRawLSH,"dx",as.numeric(rawVerts$dx))
set.vertex.attribute(IBDPairFromRawLSH,"sex",as.numeric(rawVerts$sex))
set.vertex.attribute(IBDPairFromRawLSH,"carrier",as.numeric(rawVerts$carrier))
set.vertex.attribute(IBDPairFromRawLSH,"inClst",as.numeric(rawVerts$includedInClst))


#Now for some additional (and optional) book-keeping: we can also bring in the data_id – which is the subject research id used across classrooms throughout the McFarland dataset – and set it as the “persistent id”. This allows us to easily translate from the vertex ids (which must always range from 1 to the size of the network) even when working with vertex subsets of the original network. (See ?persistent.ids)
set.vertex.attribute(IBDPairFromRawLSH,"data_id",as.numeric(rawVerts$data_id))
set.network.attribute(IBDPairFromRawLSH,'vertex.pid','data_id')

IBDPairFromRawLSH %v% "col" <- c("green", "red")[IBDPairFromRawLSH %v% "carrier"]
IBDPairFromRawLSH %v% "shape" <- IBDPairFromRawLSH %v% "dx" + 2
IBDPairFromRawLSH %v% "borderWidth" <- c(1, 5)[IBDPairFromRawLSH %v% "inClst"]

plot(network.extract(IBDPairFromRawLSH, at=), vertex.cex=(), vertex.col="col")
plot(network.extract(IBDPairFromRawLSH, onset=, terminus=, rule="any"), vertex.cex=(), vertex.col="col")

filmstrip(IBDPairFromRawLSH, displaylabels=F, mfrow=c(1, 5), slice.par=list(start=0, end=49, interval=10, aggregate.dur=10, rule='any'))

head(as.data.frame(IBDPairFromRawLSH))

compute.animation(IBDPairFromRawLSH, animation.mode = "kamadakawai", slice.par=list(start=25005030, end=54933510, interval=300000, aggregate.dur=1, rule='any'))

render.d3movie(IBDPairFromRawLSH, displaylabels=F, usearrows = F, vertex.rot = -30, vertex.lwd = IBDPairFromRawLSH %v% "borderWidth", vertex.border = "#333333",
    vertex.col=IBDPairFromRawLSH %v% "col", vertex.sides = IBDPairFromRawLSH %v% "shape", 
    edge.col='gray', edge.lwd = (IBDPairFromRawLSH %e% "weight")/5, vertex.tooltip = paste("<b>",(IBDPairFromRawLSH %v% 'data_id'), "</b>"),
   edge.tooltip = paste("<b>", (IBDPairFromRawLSH %e% 'weight'), "</b>"), launchBrowser=T, filename="RawLSH0.1-IBDgraph.html" ) 

proximity.timeline(IBDPairFromRawLSH, default.dist=6, mode='sammon',labels.at=25,vertex.cex=4)

plot(tSnaStats(IBDPairFromRawLSH,'gtrans') )

###################################################################################
#!!!!!!!!!make sure the current directory is the one that has all the graph files needed for the following operations

rm(list = ls())
library(igraph)

nodes <- read.csv("SubGraphAtPoint3703node.csv", header=T, as.is=T)
colnames(nodes) = c("ID","DX","Carrier","Promiscuous")


links <- read.csv("SubGraphAtPoint3703edge.csv", header=F, as.is=T)
colnames(links)= c("MainSub", "AssoSub", "AvgScore")
links <- unique(links)

#nrow(nodes); length(unique(nodes$ID))
#nrow(links); nrow(unique(links[,c("MainSub", "AssoSub")]))
net3703 <- graph_from_data_frame(d=links, vertices=nodes, directed=T)

#set color of vertex attribute according to node characteristic
colrs <- c("green", "red")
V(net3703)$color <- colrs[V(net3703)$Carrier+1]

#set shape of vertex attribute according to node characteristic
shps <- c("circle", "square")
V(net3703)$shape <- shps[V(net3703)$DX]


# Compute node degrees (#links) and use that to set node size:
deg <- degree(net3703, mode="out")
V(net3703)$size <- log(deg)*5

# Set edge width based on weight:
E(net3703)$width <- E(net3703)$AvgScore/5
E(net3703)$arrow.size <- 0.5

l <- layout_with_fr(net3703)
plot(net3703, layout=l)
hist(links$AvgScore)
mean(links$AvgScore)
sd(links$AvgScore)
cut.off <- mean(links$AvgScore) 
net3703.sp <- delete_edges(net3703, E(net3703)[AvgScore<cut.off])
plot(net3703.sp, layout=l)
net3703.case <- delete_vertices(net3703, V(net3703)[DX==1])
plot(net3703, layout=layout_in_circle)
net3703.carrier <- delete_vertices(net3703, V(net3703)[Carrier==0])
plot(net3703.carrier, layout=layout_in_circle)

legend(x=-1.5, y=-1.1, c("Non catrrier","Carrier"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

#find all connected components (at least 2 nodes) and transform each of them to a graph
comps <- decompose.graph(net3703)
plot(comps[[2]]) #4 lrrk2 carrier, 1 non-carrier

#color the edges of the graph based on their source node color
edge.start <- ends(net3703, es=E(net3703), names=F)[,1]
edge.col <- V(net3703)$color[edge.start]
plot(comps[[2]], edge.color=edge.col, edge.curved=.3) 


#Clique operations
cliqs = cliques(net3703) # list of cliques       
sapply(cliqs, length) # clique sizes
largest_cliques(cliqs) # cliques with max number of nodes



#Animation/dynamic network visualization
library('animation') 
library('igraph')


//merge subject (node) files
file_list <- list.files(path=".", pattern="*node\\.csv$")
allNodes <- do.call("rbind",lapply(file_list,FUN=function(files){read.table(files, header=TRUE, sep=",")}))
colnames(allNodes) = c("ID","DX","Carrier","Promiscuous")
allNodes <- unique(allNodes)
write.csv(allNodes, file="allNodes.csv",row.names=F,quote=F)

library(ndtv)
#Great tutorial: http://statnet.csde.washington.edu/workshops/SUNBELT/current/ndtv/ndtv_workshop.html
#Render a network animation as an interactive web page
render.d3movie(short.stergm.sim)
#This should bring up a browser window displaying the animation. In addition to playing and pausing the movie, it is possible to click on vertices and edges to show their ids, double-click to highlight neighbors, and zoom into the network at any point in time.

#plotting the active spells of edges and vertices. In this view, only the activity state of the network elements are shown–the structure and network connectivity is not visible. The vertices in this network are always active so they trace out unbroken horizontal lines in the upper portion of the plot, while the edge toggles are drawn in the lower portion.
timeline(short.stergm.sim)

pairwiseDemo <- network.initialize(12)
add.edges.active(pairwiseDemo,tail=1,head=c(2:5),onset=2:5,terminus=c(6,4,5,6))
add.edges.active(pairwiseDemo,tail=2,head=c(6:9),onset=7:10,terminus=10)
plot(pairwiseDemo)
render.d3movie(pairwiseDemo,out.mode='htmlWidget')



