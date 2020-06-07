# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:02:03 2020

@author: Dale
"""
#%% Run this function first.
def tmtcorrFilter(pop,threshold,raw):
    import numpy as np
    import pandas as pd
    import networkx as nx
    data = raw['cell population']==pop
    data = raw[data]
    p_val = data['spearman p-value']<=0.05
    data_pval = data[p_val]
    data_temp = data_pval.iloc[:,8].abs()
    data_corr = data_temp > threshold
    data_corr = data_pval[data_corr]
    output = data_corr.iloc[:,[1,2,8,9]]
    return output

def lfqFind(proteins,LFQ):
    import numpy as np
    import pandas as pd
    ind = np.zeros(shape=(len(proteins),6))
    for i in range(len(proteins)):
        temp = LFQ.iloc[:,0]==proteins[i]
        temp = LFQ[temp]
        ind[i,:] = temp.iloc[:,1:7]
    values = ind
    return values

def LFQplotter(values,proteins,labels,title):
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
        rects1 = ax.bar(x + partition[i]*diff,values[:,i],diff,label = labels[i])
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Summed LFQ',fontname="Arial", fontsize=12)
    ax.set_xlabel('Genes',fontname="Arial", fontsize=12)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(proteins)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    ax.legend()
    fig.tight_layout()
    saveName = title+'-fig.tiff'
    fig.savefig(saveName)
    
    plt.show()
    
def LFQplotter1(values,proteins,labels,title):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    
    x = np.arange(len(labels))  # the label locations
    width = 0.5  # the width of the bars
    diff = width/len(proteins)
    num = len(proteins)
    if (num % 2) == 0:
        shift = num/2
        partition = np.arange(num)
        partition = partition - shift
    else:
        shift = num/2-0.5
        partition = np.arange(num)
        partition = partition - shift
    fig, ax = plt.subplots(dpi=300)
    
    for i in range(len(proteins)):
        rects1 = ax.bar(x + partition[i]*diff,values[i,:],diff,label = proteins[i])
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Summed LFQ',fontname="Arial", fontsize=12)
    ax.set_xlabel('Genes',fontname="Arial", fontsize=12)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    ax.legend()
    fig.tight_layout()
    saveName = title+'-fig.tiff'
    fig.savefig(saveName)
    
    plt.show()
#%% Script to generate outputs for different cell populations.
import numpy as np
import pandas as pd

rawTMT = pd.read_csv('tmt.csv')
rawLFQ = pd.read_csv('lfq.csv')
LFQ = rawLFQ.iloc[:,[0,11,12,13,14,15,16]]

#%% For TMT data correlations
thresh = 0.5
ERP_output = tmtcorrFilter('ERP',thresh,raw)
GRA_output = tmtcorrFilter('GRA',thresh,raw)
HSC_output = tmtcorrFilter('HSC',thresh,raw)
LYM_output = tmtcorrFilter('LYM',thresh,raw)
MON_output = tmtcorrFilter('MON',thresh,raw)
MSC_output = tmtcorrFilter('MSC',thresh,raw)

s1 = set(ERP_output.iloc[:,0])
s2 = set(GRA_output.iloc[:,0])
s3 = set(HSC_output.iloc[:,0])
s4 = set(LYM_output.iloc[:,0])
s5 = set(MON_output.iloc[:,0])
s6 = set(MSC_output.iloc[:,0])

s12 = s1.intersection(s2)
s13 = s1.intersection(s3)
s14 = s1.intersection(s4)
s15 = s1.intersection(s5)
s16 = s1.intersection(s6)
s23 = s2.intersection(s3)
s24 = s2.intersection(s4)
s25 = s2.intersection(s5)
s26 = s2.intersection(s6)
s34 = s3.intersection(s4)
s35 = s3.intersection(s5)
s36 = s3.intersection(s6)
s45 = s4.intersection(s5)
s46 = s4.intersection(s6)
s56 = s5.intersection(s6)

#%%
proteins = ['DDB2','VPRBP','DCAF5','DCAF6','DCAF7','DCAF8','DCAF11','DCAF12','DCAF13']
labels = ['HPC','LYM','MON','GRA','ERP','MSC']
values = lfqFind(proteins,LFQ)
LFQplotter(values,proteins,labels,'CUL4B Substrate Receptors')
LFQplotter1(values,proteins,labels,'CUL4B Sub. Recept._1')

proteins = ['CUL5','ARIH2']
labels = ['HPC','LYM','MON','GRA','ERP','MSC']
values = lfqFind(proteins,LFQ)
LFQplotter(values,proteins,labels,'ARIH2')
LFQplotter1(values,proteins,labels,'ARIH2_1')
