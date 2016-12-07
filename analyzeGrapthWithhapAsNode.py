'''
Created on Mar 26, 2016

@author: linly
'''
import os
import argparse
#import snap
#import the sys package to use one of its methods later
import sys
#append the path of the parent directory of your package
sys.path.append("/home/linly/.local/lib/python2.6/site-packages/")
#import your package
import networkx as nx
import matplotlib.pyplot as P
import gc
import igraph

def main() :
    # parse the command line
    parser=argparse.ArgumentParser()
    parser.add_argument("-m", "--markPos", help="specifies the seed marker position of the haplotypes used to draw the graph", action="store", required=True)
    args=parser.parse_args()
    
    #g = nx.read_weighted_edgelist(os.getcwd() + "lrrkPairsAtSite_" + args.markPos + ".out", delimiter=",", nodetype=int)
    #nx.write_pajek(g, os.getcwd()+args.markPos+".net")
    g = nx.read_weighted_edgelist("/projects/sequence_analysis/vol4/CHAT_simGWAS/Chat1.2Lrrk2Graph/lrrkPairsAtSite_3999.out",comments='h',delimiter=",", nodetype=str)
    nx.write_pajek(g, "/projects/sequence_analysis/vol4/CHAT_simGWAS/Chat1.2Lrrk2Graph/lrrk_ego/lrrkPairsAtSite_3999.net")
    gc.collect()#garbage collection
    del g
    g = nx.read_pajek("/projects/sequence_analysis/vol4/CHAT_simGWAS/Chat1.2Lrrk2Graph/lrrkPairsAtSite_8003.net")
    
    #assign graph attributes
    #g.graph['marker'] = args.markPos
    g.graph['marker'] = '3999'
    
    #use nx.info() to see whether the network has been loaded correctly
    #display basic information of the network (num of nodes, edges, average degree)
    nx.info(g)
    
    #calculate degree distribution
    dh=nx.degree_histogram(g)
    #Plot using same method as http://networkx.lanl.gov/examples/drawing/degree_histogram.html
    #pos=nx.spring_layout(g)
    P.figure(figsize=(8,8))
    P.loglog(dh,'b-',marker='o')
    P.title("Degree distribution (log-log)")
    P.ylabel("Degree")
    P.xlabel("Frequency")
    #Draw graph in inset
    #P.axes([0.45,0.45,0.45,0.45])
    #P.axis('off')
    #nx.draw_networkx_nodes(g,pos,node_size=10)
    #nx.draw_networkx_edges(g,pos,alpha=0.4)
    P.savefig("8003_degDist.png")
    P.clf()
    #igraph version
    
    
    #Generate a sorted list of connected components, largest first.
    [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    
    #Get the largest connected component
    largest_cc = max(nx.connected_components(G), key=len)
    
    from itertools import combinations
    top_overlap = [list(combinations(c, 2)) for c in components if len(c) > 1]
    top_overlap = [item for sublist in top_overlap for item in sublist]
    #top_overlap
    #[(1, 3), (1, 2), (3, 2)]
    
    #Caution: the following step is memory intensive
    c = list(nx.k_clique_communities(g, 3))
    
    
    

if __name__ == '__main__':
    main()