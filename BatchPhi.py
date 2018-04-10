import os
import subprocess

folder = os.listdir("alignments_fasta")
oup = open("Phi_output.txt","w")

for i in folder:
    process = subprocess.Popen("Phi -p -f alignments_fasta/%s"%i, shell = True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    
    lines = process.stdout.readlines()
    pvalue = "NA"
    
    for line in lines:
        if "(Permutation)" in line:
            pvalue = line.split(":")[1].split("(")[0]
            pvalue = pvalue.replace(" ","")
            pvalue = float(pvalue)
            print i,pvalue
        else:
            pass
        
    oup.write("%s\t%s\n"%(i,pvalue))
    
    
oup.close()

    
    
