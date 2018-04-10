from Bio.SeqUtils import GC
from Bio import SeqIO
from sys import argv
import os
from Bio.SeqFeature import SeqFeature,FeatureLocation
import numpy as np

script, genome_folder = argv

folder = os.listdir(genome_folder)

oupGC = open("GCcontent.txt","w")
oupSize = open("genome_sizes.txt","w")

oupGC.write("DATASET_BOXPLOT\nSEPARATOR COMMA\nDATASET_LABEL,GC content\nCOLOR,#ff0000\nDATA\n")
oupSize.write("DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Size\nCOLOR,#ff1000\nDATA\n")
    
for inl in folder:
    print inl
    GClist = []
    size = 0
    genome = SeqIO.parse("%s/%s"%(genome_folder,inl),"fasta")
    
    #dividing genome sequence in 1kb windows
    for contig in genome:
        size = size + len(contig.seq)
        i=0
        
        while (i +1000)< len(contig.seq):
            start = i
            end = i + 1000
            seqGC = contig.seq[start:end]
            if float(GC(seqGC))>0:
                GClist.append(float(GC(seqGC)))
            else:
                pass
            i= i + 1001
    
    GClist = sorted(GClist)
    GClist = GClist[50:-50] #to remove extremes
    
    if GClist != []:
        Q1 = np.percentile(GClist,25)
        Q3 = np.percentile(GClist,75)
        median = np.percentile(GClist,50)
        maximum = np.percentile(GClist,99)
        minimum = np.percentile(GClist,1)
        
        #name = ""
        #for feat in genome.features:
            #if feat.type == "source":
                #try:
                    #name = feat.qualifiers["organism"][0]
                    #name = name.replace(" ","_")
                #except KeyError:
                    #pass
            #else:
                #pass
        
        oupGC.write("%s,%.2f,%.2f,%.2f,%.2f,%.2f\n"%(inl.split(".")[0],minimum,Q1,median,Q3,maximum))
        oupSize.write("%s,%.2f\n"%(inl.split(".")[0],float(size)/float(1e6)))

    else:
        pass
             
oupGC.close()
oupSize.close()
