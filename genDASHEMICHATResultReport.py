import os
import fnmatch 
import argparse
import gzip
import numpy as np
import re
from pylab import *

cliqueinfo={}

def main() :
    # parse the command line
    parser=argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", help="absolute path to the directory that stores the result file, which ends with .clstPV", action="store", required=True)
    parser.add_argument("-p", "--prog", help="specify the program that produces this result", action="store", required=True)
    #parser.add_argument("-m", "--marker", help="absolute path to the file that stores the bp of markers; must provided when using chat", action="store", required=False)
    parser.add_argument("-o", '--outpre', help="output file prefix", action="store", required=False)
    parser.add_argument("-v", "--verbose", help="verbose",action="store",required=False)
    args=parser.parse_args() 
    
    #if args.prog.lower().startswith("chat"):
        #markerFile = open(args.marker)
        #markerList=[]
        #for line in markerFile.readlines():
            #markerList.append(long(line))
            
    foundFile=False
    for dirName, subdirList, fileList in os.walk(args.dir):
      for x in fileList:
        if (x.endswith(".clstPV") or x.startswith("GenomeAdjustClique")):
          foundFile=True    
          processResultfile(os.path.join(dirName,x), args.prog)
          # then write out report
          if args.outpre is None:
            writereport(args.dir, args.prog, 0, 23, "IBDCluster")
          else:
            writereport(args.dir, args.prog, 0, 23,args.outpre+"_IBDCluster")
    if foundFile is not True:
       print "No result file."          

# write report
def writereport(workDir, usedProgram, vminColorBar, vmaxColorBar, reportPrefix) :  
    cliquelist=cliqueinfo.keys()
    cliquelist.sort(key=int)
    
    pltx_sPos=[]
    pltx_ePos=[]
    plty=[]
    caseCt=[]
    
    for clique in cliquelist :
        pltx_sPos.append(np.mean(cliqueinfo[clique]['startPos']))
        if not usedProgram.lower().startswith("chat"):
            pltx_ePos.append(np.mean(cliqueinfo[clique]['endPos']))
        plty.append(np.mean(cliqueinfo[clique]['pv']))
        caseCt.append(np.mean(cliqueinfo[clique]['casect']))
    if usedProgram.lower().startswith("chat"):
        plotResult(pltx_sPos,plty,caseCt,usedProgram,"Fisher's exact test PVs",vminColorBar,vmaxColorBar,38673387,39399147, reportPrefix)#the dataset chat uses follows hg18/build36
    else:
        lineSeg(pltx_sPos,pltx_ePos,plty,caseCt,usedProgram,"Fisher's exact test PVs",True,vminColorBar,vmaxColorBar,40387120,41112880, reportPrefix)#the dataset dash/emi uses follows hg19/build37

    
# Draw scatter plot with jitter (each clique is represented by the start marker of the locating bin)
def plotResult(x,y,caseCt,usedProgram,plotTitle,vminColorBar,vmaxColorBar,regionStart,regionEnd,reportPrefix) :
    # s needs to be big enough to cover a region as each node indicates the start marker of a bin; the line plot function below (i.e., lineSeg) provides more accuracy
    #sc=plt.scatter(rand_jitter(x), rand_jitter(y), s=50, c=colorType, cmap="gray")
    cm = plt.cm.get_cmap('jet')
    sc = plt.scatter(x, y, c=caseCt, vmin=vminColorBar, vmax=vmaxColorBar, s=35, cmap=cm)
    plt.colorbar(sc, ticks=range(0,24,4)) 
    #plt.scatter(x, y, s=30, c=colorType, cmap="gray")#no jitter
    #plt.xlabel('Start Marker Index')
    plt.xlabel('Marker Position')
    plt.ylabel('Fisher-exact-test P value')
    plt.ylim(0,12)
    plt.title(usedProgram + ": " + plotTitle)
    plt.axvspan(regionStart, regionEnd, facecolor='grey', alpha=0.5)
    #plt.xlim([0,20000000])
    #plt.xlim([4525000,11350000])
    #plt.axvline(x=cv, color='red', linewidth=2)
    #plt.axhline(y=1.3, color='blue', linewidth=2)
    #plt.savefig("IBDCluster_" + usedProgram + ".png", format='png')
    plt.savefig(reportPrefix + ".png", format='png')
    plt.clf()

def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev   
    
# Draw line graph (each clique is represented by a horizontal line that connects the start and end markers of the current bin)
def lineSeg(xmin,xmax,y,caseCt,usedProgram,plotTitle,isBPPos,vminColorBar,vmaxColorBar,regionStart,regionEnd, reportPrefix) :
    # create a ScalarMappable and initialize a data structure (cmap -- color map, norm is a class which, when called, can normalize data into the [0.0, 1.0] interval)
    # norm=plt.Normalize(vmin=np.min(caseCt), vmax=np.max(caseCt))
    norm=plt.Normalize(vmin=vminColorBar, vmax=vmaxColorBar)
    sm = plt.cm.ScalarMappable(cmap="jet", norm=norm)
    sm.set_array([])
    
    # plotting line segments; the colors are chosen by calling the ScalarMappable that was initialised with cmap and norm
    #plt.hlines(rand_jitter(y),xmin,xmax,colors=sm.to_rgba(caseCt))
    plt.hlines(y,xmin,xmax,colors=sm.to_rgba(caseCt))
    # having plotted the 11 curves we plot the colorbar, using again our ScalarMappable
    plt.colorbar(sm, ticks=range(0,13,2))
    plt.margins(0.01)# So the lines aren't at the plot boundaries.
    if isBPPos :
        plt.xlabel('Marker Position')
        #plt.axvline(x=regionStart, color='grey', linewidth=1)
        #plt.axvline(x=regionEnd, color='grey', linewidth=1)
        plt.axvspan(regionStart, regionEnd, facecolor='grey', alpha=0.5)
    else :
        plt.xlabel('Marker Index')
        plt.axvspan(2570, 2571, color='red', alpha=0.5)
    #plt.xlim([plotXmin,plotXmax])    
    #plt.axhline(y=1.3, color='blue', linewidth=1)
    plt.ylabel('Fisher-exact-test P value')
    plt.ylim(range(0,13,2))
    plt.title(usedProgram + ": " + plotTitle)
    if isBPPos :
      plt.savefig(reportPrefix + "_" + usedProgram + ".png", format='png')
    else :
      plt.savefig(reportPrefix + "_" + repID+"_" + modelName+"_l_xid_" + str(plotXmin) + "_" + str(plotXmax) + ".png", format='png')
    plt.clf()   

# collect clique info
def processResultfile(filename, usedProgram) :
    print "Processing file " + filename 
    cliqueinfo.clear()
    cliqueID=0
    if usedProgram.lower().startswith("chat"):
        infile=gzip.open(filename,'r')
    else:
        infile = open(filename)    
    for line in infile.readlines():
       if(line.startswith("lct")):
          break
       linetoks=line.split('\t')
       #print linetoks[2]
       if cliqueID not in cliqueinfo:
           cliqueinfo[cliqueID]={}
           cliqueinfo[cliqueID]['startPos']=[]
           cliqueinfo[cliqueID]['endPos']=[]
           cliqueinfo[cliqueID]['pv']=[]
           cliqueinfo[cliqueID]['casect']=[]
       clique=cliqueinfo[cliqueID]
       if not usedProgram.lower().startswith("chat"):
           clique['startPos'].append(float(linetoks[1]))
           clique['endPos'].append(float(linetoks[2]))
           clique['casect'].append(float(linetoks[3]))
           clique['pv'].append(float(linetoks[5]))
       else:
           clique['startPos'].append(float(linetoks[2]))
           clique['casect'].append(float(linetoks[5]))
           clique['pv'].append(float(linetoks[7]))
       cliqueID += 1
    infile.close()
      
if __name__ == "__main__":
    main()
