from sys import argv
import os

script,folder, control_template = argv

os.system("mkdir results")
files = os.listdir(folder)

for inl in files:
    line = "      seqfile = %s/%s\n"%(folder,inl)
    ctrl = open(control_template,"rU")
    temp = open("yn00.temp.ctl","w")
    
    for i in ctrl.readlines():
        if "seqfile" in i:
            temp.write(line)
        elif "outfile" in i:
            temp.write("      outfile = results/yn_%s\n"%inl)
        else:
            temp.write(i)
    temp.close()
    os.system("yn00 yn00.temp.ctl")
    
