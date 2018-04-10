from Bio import AlignIO
from itertools import groupby
from operator import itemgetter
from sys import argv
import os
from Bio.SeqUtils import GC

script, align_folder = argv

folder = os.listdir(align_folder)


#oup_ins = open("ins_statistics.txt","w")
oup_GC = open("GCstatistics.txt","w")

del_list = []
ins_list = []
num_indels = []
GClist_pseudo = []
GClist_funct = []
identities = []
main_dict = {}
ortho_groups= []

def percentID(seq1,seq2):
    x = 0
    for i,j in enumerate(seq1):
        if seq2[i] == j:
            x = x+1
        else:
            pass
    percID = 100*float(x)/len(seq1)
    
    return percID



species_list = []

for inl2 in folder:
    align2 = AlignIO.read("%s/%s"%(align_folder,inl2),"clustal")
    for seq2 in align2:
        species_list.append(seq2.id)

species_list = set(species_list)

oup_GC.write("group")
for sp in species_list:
    main_dict[sp] = {}
    #oup_GC.write("\t%s"%sp)
    
oup_GC.write("group\tspecies\tGC\tGC1\tGC2\tGC3\n")

for inl in folder:
    print inl
    group = inl.split(".aln")[0]
    ortho_groups.append(group)
    os.system("trimal -in %s/%s -out trimmed/%s.aln -gt 0.9"%(align_folder,inl,inl)) #generate a trimmed alignment
    # to calculate GC content on conserved sections
    
    align = AlignIO.read("trimmed/%s.aln"%inl,"clustal")
    
   
    for seq in align:
        species = seq.id
        seqstring = str(seq.seq).replace("-","")
        GCcont = GC(seqstring)
        GC1 = GC(seqstring[0::3])
        GC2 = GC(seqstring[1::3])
        GC3 = GC(seqstring[2::3])
        

        main_dict[species][group] = (GCcont,GC1,GC2,GC3)
        print species, group, GCcont
    

    

ortho_groups = set(ortho_groups)

for i in ortho_groups:
    for spec in species_list:
        oup_GC.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\n"%(i,spec,main_dict[spec][i][0],main_dict[spec][i][1],main_dict[spec][i][2],main_dict[spec][i][3]))

oup_GC.close()
        
        
