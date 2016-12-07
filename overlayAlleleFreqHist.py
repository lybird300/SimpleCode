import os
import fnmatch 
import argparse
import gzip
import numpy as np
import re
import pandas as pd
from pylab import *

freqFileList={}
#for ctrlSize in ['100', '500', '1000', '2000']:
#    freqFileList[ctrlSize]={}
#cv = 7734610

common_params = dict(bins=20, range=(0,0.6), histtype='step', color='grey')

def main() :
    # parse the command line
    parser=argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", help="absolute path of the cv directory; sth like workDir/SNP_<cv position>_<cv maf range>", action="store", required=True)
    parser.add_argument("-v", '--verbose', help="verbose",action="store",required=False)
    args=parser.parse_args()
    
    # now walk through the directory tree looking for freq files (i.e., Chr08majorAlleleFreq.txt)
    ctrlCt=''
    chromKey=''
    for dirName, subdirList, fileList in os.walk(args.indir):
        subdirList[:] = [x1 for x1 in subdirList if x1.startswith(('Case', 'ChatInput', 'Rep', 'ChromosomeSpecific'))]
        for x in subdirList:
            if x.startswith('Case'):
                ctrlCt = x.split('_')[1]
                #print ctrlCt
                if ctrlCt not in freqFileList:
                    freqFileList[ctrlCt]={}
                    #print freqFileList.keys()
            if x.startswith("Rep"):
                for repSubDirName in os.listdir(os.path.join(dirName,x)):
                    #print repSubDirName
                    if repSubDirName.startswith('ChromosomeSpecificData'):
                        collectFreqFilePaths(os.path.join(dirName,x,repSubDirName), ctrlCt)
                        break

    ctrlCt=''
    chromKey=''
    printAlleleFreqSpectrum(ctrlCt, chromKey)

def collectFreqFilePaths(chatInputDirComplete, numOfCtrls) :
    for x2 in os.listdir(chatInputDirComplete):
        if x2.endswith('majorAlleleFreq.txt'):
            chromKey = x2[3:5]
            #print "chromKey for File ", x2, " is ", chromKey
            if chromKey not in freqFileList[numOfCtrls]:
                freqFileList[numOfCtrls][chromKey]=[]
                #print freqFileList[numOfCtrls].keys()
            freqFileList[numOfCtrls][chromKey].append(os.path.join(chatInputDirComplete, x2))

        
def printAlleleFreqSpectrum(numOfCtrls, chromID) :
    fig,ax=plt.subplots()
    for ctrlCt in freqFileList:
        if(numOfCtrls and ctrlCt!=numOfCtrls):
            continue
        for chromKey in freqFileList[ctrlCt]:
            if(chromID and chromKey!=chromID):
                continue          
            for f in freqFileList[ctrlCt][chromKey]:
                data = loadtxt(f, delimiter='\t', skiprows=1, usecols=(1,))
                data[:] = [1-majAlleFreq for majAlleFreq in data]
                #df = pd.DataFrame(data)
                #df.plot(kind='density',ax=ax,legend=False,xlim=(0,0.6), color='grey')    
                plt.hist(data, **common_params)
    popData = loadtxt("/projects/sequence_analysis/vol4/CHAT_simGWAS/newSimGWASData/SNP_7734610_1.0E-5_1.0E-4/tagSNPMAF.txt", delimiter='\t', skiprows=1, usecols=(1,))
    #df = pd.DataFrame(popData)
    #df.plot(kind='density',ax=ax,legend=False,xlim=(0,0.6), linewidth=3, color='black')
    plt.hist(data, bins=20, range=(0,0.6), histtype='step', color='black')
    plt.title('Tag SNPs minor allele frequency distribution')
    plt.savefig('majAlleFreqOverlay_hist.png')
    #plt.show()

if __name__ == "__main__":
    main()