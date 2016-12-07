#
# Generate a report on ChatSetGenomeAdjust results
#  GenomeAdjustChatSet.gz file columns
#			0: chrom ID
#			1: analysis point ID
#     2: analysis point base pair position
#			3: analysis point rs name
#			4: clique ID
#			5: ct of cases in the clique
#			6: ct of controls in the clique
#			7: -log10(P value of Fisher's exact test on the clique)
#			8: -log10(P value after local adjustment)
#			9: -log10(P value after genome-wide adjustment)
#			10: the set of subjects in the clique (subject IDs separated by ",")
#			11: the set of haplotypes in the clique (haplotype hash IDs separated by ",")
#			10: the threshold used when forming the current clique

import os
import fnmatch 
import argparse
import gzip
import numpy as np
import re
from pylab import *

cliqueinfo={}
regionStart = 3703
regionEnd = 3969
#for simulated data
#regionStart = 2405
#regionEnd = 2406


def main() :
    # parse the command line
    parser=argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", help="example: python <absolutePathToScript>/scriptName.py -i workDir/../CliqueGenomeAdjust/PP_0.001/ -o workDir/../CliqueGenomeAdjust/PP_0.001/", action="store", required=True)
    parser.add_argument("-v", '--verbose', help="verbose",action="store",required=False)
    parser.add_argument("-o", '--outpre', help="output file prefix", action="store", required=False)
    args=parser.parse_args()
    
    for dirName, subdirList, fileList in os.walk(args.indir):
        for x in fileList:
            if x.startswith("GenomeAdjustClique"):
                completePath = os.path.join(dirName,x)     
                findSth = processResultfile(completePath)
                # then write out report
                if(findSth):
                    #pass just the last directory name (model name) to writereport function  
                    if args.outpre is None:
                        writereport(args.indir,"")
                    else:
                        writereport(args.indir,args.outpre+"_")
                else:
                    print "Nothing has been found."             

# write report
def writereport(workDir, reportPrefix) :  
    cliquelist=cliqueinfo.keys()
    cliquelist.sort(key=int)   
    pltx=[]
    plty=[]
    caseCt=[]
    for clique in cliquelist :
        pltx.append(cliqueinfo[clique]['analysisPoint'])
        #plty.append(np.mean(cliqueinfo[clique]['pvs']))
        #plty.append(np.mean(cliqueinfo[clique]['localadjpvs']))
        plty.append(cliqueinfo[clique]['genomeadjpvs'])
        caseCt.append(cliqueinfo[clique]['casect'])
    nojitter(pltx,plty,caseCt,workDir,reportPrefix,False)


def nojitter(x,y,caseCt,workDir,reportPrefix,jitter) :
    if jitter is True:
        sc=plt.scatter(rand_jitter(x), rand_jitter(y), s=20, c=caseCt, cmap="gray")#jitter
    else:
        sc=plt.scatter(x, y, s=30, c=caseCt, cmap="gray")#no jitter
    plt.colorbar(sc)
    plt.xlabel('Analysis Point Index')
    #plt.xlabel('Start Marker Position')
    #plt.ylabel('Fisher-exact-test P value')
    #plt.ylabel('Local ajusted -log10 p value')
    plt.ylabel('Global ajusted -log10P')
    #plt.title(reportPrefix + ": Model " + modelName)
    plt.xlim([0,8030]) #included marker index [0, 8026]
    #plt.ylim([0,1])
    plt.axvspan(regionStart, regionEnd, facecolor='red', alpha=0.5)
    plt.axhline(y=1.3, color='blue')
    #plt.savefig(modelName+'_pv.png', format='png')
    #plt.savefig("Ca1k_Co1k_grr10_rep0_"+modelName+"_"+generation+"_localpv.png", format='png')
    plt.savefig(os.path.join(workDir, reportPrefix+"genomepv.png"), format='png')
    #plt.savefig(os.path.join(workDir, modelName+"_localpv.png", format='png')
    plt.clf()
    
def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

# update clique info
def processResultfile(filename) :
    print "Processing file " + filename
    infile=gzip.open(filename,'r')   
    cliqueinfo.clear()
    cliqueID=0
    for line in infile.readlines():
       if(line.startswith("lct")):
           break
       linetoks=line.split()
       if cliqueID not in cliqueinfo:
           cliqueinfo[cliqueID]={}
           cliqueinfo[cliqueID]['analysisPoint']=[]
           cliqueinfo[cliqueID]['pvs']=[]
           cliqueinfo[cliqueID]['localadjpvs']=[]
           cliqueinfo[cliqueID]['genomeadjpvs']=[]
           cliqueinfo[cliqueID]['casect']=[]
       clique=cliqueinfo[cliqueID]
       clique['analysisPoint'].append(float(linetoks[1]))
       clique['casect'].append(float(linetoks[5]))
       clique['pvs'].append(float(linetoks[7]))
       clique['localadjpvs'].append(float(linetoks[8]))
       clique['genomeadjpvs'].append(float(linetoks[9]))
       cliqueID += 1
    infile.close()
    if(cliqueID==0):
        return False
    return True
    
        
if __name__ == "__main__":
    main()