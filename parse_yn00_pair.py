import os
from sys import argv
from Bio.Phylo.PAML import yn00

script, folder = argv

folder_in = os.listdir(folder)

oup = open("pairwise_dNdS.txt","w")
oup.write("file\tpair\tdNdS\tdN\tdS\n")

for inl in folder_in:
    print inl
    result = yn00.read("%s/%s"%(folder,inl))
    
    for k in result.keys():
        if "BSIN" in k:
            dNx = dSx = Wx = 0
            for x in result[k].keys():
                if "WS68" in x:
                    dNx = result[k][x]["YN00"]["dN"]
                    dSx = result[k][x]["YN00"]["dS"]
                    Wx = result[k][x]["YN00"]["omega"]
                    oup.write("%s\tBSIN_WS68\t%r\t%r\t%r\n"%(inl,Wx,dNx,dSx))
                else:
                    pass
                
        elif "AQ15" in k:
            dNy = dSy = Wy = 0
            for x in result[k].keys():
                if "BG90" in x:
                    dNy = result[k][x]["YN00"]["dN"]
                    dSy = result[k][x]["YN00"]["dS"]
                    Wy = result[k][x]["YN00"]["omega"]
                    oup.write("%s\tAQ15_BG90\t%r\t%r\t%r\n"%(inl,Wy,dNy,dSy))
                else:
                    pass
        elif "BTN" in k:
            dNy = dSy = Wy = 0
            for x in result[k].keys():
                if "BG90" in x:
                    dNz = result[k][x]["YN00"]["dN"]
                    dSz = result[k][x]["YN00"]["dS"]
                    Wz = result[k][x]["YN00"]["omega"]
                    oup.write("%s\tBTN_BG90\t%r\t%r\t%r\n"%(inl,Wz,dNz,dSz))
                else:
                    pass
        else:
            pass
oup.close()
                
              
