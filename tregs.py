# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:12:35 2020

@author: Dale
"""

## Some Notes:
    #To run this script, you should have three files in your directory:
        #tregs.py
        #GSE76598_family.soft
        #probe2gene571.txt
    #This script is split into three cells. Run the first cell (Ctrl+Enter) to
    #establish functions, the second cell to process data from GEO file, and the
    #third cell to choose proteins, plot, and save graph to current directory.
    
## IMPORTANT NOTE: You need to install Biopython to run this script. To install Biopython,
#run the following line (use F9 to run one line):
    
pip install biopython

#%% Run this cell first
def makeDataframe(data):
    import pandas as pd
    output = pd.DataFrame.from_dict(data[0],orient = 'index')
    for i in range(len(data)-1):
        df = pd.DataFrame.from_dict(data[i+1],orient = 'index')
        output[i+1] = df
        
    return output

def expPlotter(values,proteins,labels,title,std):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    
    x = np.arange(len(proteins))  # the label locations
    width = 0.5  # the width of the bars
    diff = width/len(labels)
    num = len(labels)
    if (num % 2) == 0:
        shift = num/2
        partition = np.arange(num)
        partition = partition - shift
    else:
        shift = num/2-0.5
        partition = np.arange(num)
        partition = partition - shift
    fig, ax = plt.subplots(dpi=300)
    
    for i in range(len(labels)):
        rects1 = ax.bar(x + partition[i]*diff,values[:,i],diff,label = labels[i],yerr=std[:,i],
       align='center',
       alpha=1,
       ecolor='black',
       capsize=4)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Raw RNA Measure',fontname="Arial", fontsize=12)
    ax.set_xlabel('Genes',fontname="Arial", fontsize=12)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(proteins)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    ax.legend(prop={'size': 8})
    fig.tight_layout()
    saveName = title+'-fig.tiff'
    fig.savefig(saveName)
    
    plt.show()
#%% Run this cell next
    import pandas as pd
    import numpy as np
    from Bio import Geo

    sample_characteristics = {}
    handle = open('GSE76598_family.soft')
    records = Geo.parse(handle)
    records2 = []
    for record in records:
        records2.append(record)
    x = records2[3].entity_attributes['Sample_characteristics_ch1']
    for i in x:
        q = i.split(':')
        sample_characteristics[q[0]] = []
    for i in records2[3:]:
        x = i.entity_attributes['Sample_characteristics_ch1']
        for j in x:
           q = j.split(':')
           sample_characteristics[q[0]].append(q[1])
    counter = 0
    naiveTregUnstim = []
    naiveTregStim = []
    memoryTregUnstim = []
    memoryTregStim = []
    naiveTconvUnstim = []
    naiveTconvStim = []
    memoryTconvUnstim = []
    memoryTconvStim = []
    x = sample_characteristics['cell subset']
    x1 = sample_characteristics['treatment']
    x = [i.strip() for i in x]
    x1 = [i.strip() for i in x1]
    
    for i in range(len(x)):
        if x[i] == 'naiveTreg' and x1[i] == 'unstimulated':
            naiveTregUnstim.append(i)
        elif x[i] == 'naiveTreg' and x1[i] == 'stimulated':
            naiveTregStim.append(i)
        elif x[i] == 'memoryTreg' and x1[i] == 'unstimulated':
            memoryTregUnstim.append(i)
        elif x[i] == 'memoryTreg' and x1[i] == 'stimulated':
            memoryTregStim.append(i)
        elif x[i] == 'naiveTconv' and x1[i] == 'unstimulated':
            naiveTconvUnstim.append(i)
        elif x[i] == 'naiveTconv' and x1[i] == 'stimulated':
            naiveTconvStim.append(i)
        elif x[i] == 'memoryTconv' and x1[i] == 'unstimulated':
            memoryTconvUnstim.append(i)
        elif x[i] == 'memoryTconv' and x1[i] == 'stimulated':
            memoryTconvStim.append(i)

    probe2gene = open('probe2gene571.txt')
    probe2gene = probe2gene.read()
    probe2gene = probe2gene.split(',')
    probe2geneDict = {}
    for i in probe2gene:
        i = i.split(':')
        i[0] = i[0].replace(' ','')
        i[1] = i[1].replace(' ','')
        i[0] = i[0].replace("'","")
        i[1] = i[1].replace("'","")
        probe2geneDict[i[0]] = i[1]
    Data = []
    for i in records2[3:]:
        x = i.table_rows
        ii = []
        for j in x[1:]:
            j[0] = probe2geneDict[j[0]]
            ii.append(j)
        Data.append(ii)
    naiveTregUnstim = [Data[i] for i in naiveTregUnstim] 
    naiveTregStim = [Data[i] for i in naiveTregStim]
    memoryTregUnstim = [Data[i] for i in memoryTregUnstim] 
    memoryTregStim = [Data[i] for i in memoryTregStim]
    naiveTconvUnstim = [Data[i] for i in naiveTconvUnstim] 
    naiveTconvStim = [Data[i] for i in naiveTconvStim]
    memoryTconvUnstim = [Data[i] for i in memoryTconvUnstim] 
    memoryTconvStim = [Data[i] for i in memoryTconvStim]
    
    naiveTregUnstimDict = []
    for i in naiveTregUnstim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        naiveTregUnstimDict.append(x)
    naiveTregStimDict = []
    for i in naiveTregStim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        naiveTregStimDict.append(x)
    memoryTregUnstimDict = []
    for i in memoryTregUnstim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        memoryTregUnstimDict.append(x)
    memoryTregStimDict = []
    for i in memoryTregStim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        memoryTregStimDict.append(x)
    naiveTconvUnstimDict = []
    for i in naiveTconvUnstim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        naiveTconvUnstimDict.append(x)
    naiveTconvStimDict = []
    for i in naiveTconvStim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        naiveTconvStimDict.append(x)
    memoryTconvUnstimDict = []
    for i in memoryTconvUnstim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        memoryTconvUnstimDict.append(x)
    memoryTconvStimDict = []
    for i in memoryTconvStim:
        x = {}
        for j in i:
            x[j[0]] = j[1]
        memoryTconvStimDict.append(x)
    NaiveTregUnstim = makeDataframe(naiveTregUnstimDict)
    NaiveTregStim = makeDataframe(naiveTregStimDict)
    MemoryTregUnstim = makeDataframe(memoryTregUnstimDict)
    MemoryTregStim = makeDataframe(memoryTregStimDict)
    NaiveTconvUnstim = makeDataframe(naiveTconvUnstimDict)
    NaiveTconvStim = makeDataframe(naiveTconvStimDict)
    MemoryTconvUnstim = makeDataframe(memoryTconvUnstimDict)
    MemoryTconvStim = makeDataframe(memoryTconvStimDict)
    
#%% Change the proteins to look at (use list of strings as shown below) and change save name to change file name.
proteins = ['CUL5','SOCS1','SOCS2','SOCS3','CISH','ASB1','ASB6']
saveName = 'Cullin 5 Substrate Receptors'
values = np.zeros(shape = (len(proteins),4))
std = np.zeros(shape=(len(proteins),4))
for i in range(len(proteins)):
    values[i,0] = pd.Series.mean(pd.to_numeric(NaiveTconvUnstim.loc[proteins[i]]))
    std[i,0] = pd.Series.std(pd.to_numeric(NaiveTconvUnstim.loc[proteins[i]]))
    values[i,1] = pd.Series.mean(pd.to_numeric(NaiveTconvStim.loc[proteins[i]]))
    std[i,1] = pd.Series.std(pd.to_numeric(NaiveTconvStim.loc[proteins[i]]))
    values[i,2] = pd.Series.mean(pd.to_numeric(MemoryTconvUnstim.loc[proteins[i]]))
    std[i,2] = pd.Series.std(pd.to_numeric(MemoryTconvUnstim.loc[proteins[i]]))
    values[i,3] = pd.Series.mean(pd.to_numeric(MemoryTconvStim.loc[proteins[i]]))
    std[i,3] = pd.Series.std(pd.to_numeric(MemoryTconvStim.loc[proteins[i]]))
    
labels = ['Conv_Unstim','Conv_40hr_Stim','Memory_Unstim','Memory_40hr_Stim']
expPlotter(values,proteins,labels,saveName,std)
