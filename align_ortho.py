from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from sys import argv
from Bio.Align.Applications import MuscleCommandline
import subprocess
import sys, os
import random
from Bio import AlignIO

#usage: python whole_genome_align.py orthomcl_out.txt prots.fa nucl.fa 
#orthomcl_out.txt is a set of one_to_one orthologous genes following the same format as the output of orthomcl v1.2
#prots.fa is a fasta file of protein sequences for all the genomes used in the orthomcl run. The sequence identifiers have to match the identifiers in orthomcl_out.txt
#nucl.fa is as prots.fa but with nucleotide sequences. This file is used as the reference for the back translation of the alignments.

current_dir = os.getcwd()

try:
	script, orthomcl_out, fasta_prot, fasta_nucl = argv
except ValueError:
	script, orthomcl_out, fasta_prot = argv

if len(argv) < 4:
	run_tcoffee = 0
else:
	run_tcoffee = 1

fasta_prot_dict = SeqIO.index(fasta_prot,"fasta")

if run_tcoffee == 1:
	fasta_nucl_dict = SeqIO.index(fasta_nucl,"fasta")
else:
	pass


def fasta_extract(fasta_dict,IDlist):
	#where IDlist is a list of tuples ("species","ID")
	
	protlist = []
	#sorting list of IDs per species name
	#IDlist = sorted(IDlist, key= lambda species: species[0])

	
	for x in IDlist:
		species = x[0]
		seq = fasta_dict[x[1]].seq
		seqrec = SeqRecord(seq, name=species, id=species)
		protlist.append(seqrec)
	
	return protlist
	


ortho_out = open(orthomcl_out, "rU")

line_in = ortho_out.readline()
taxa_numbers = []

while line_in != "":
	num_taxa = int(line_in.split(",")[1].split(" ")[0])
	taxa_numbers.append(num_taxa)
	line_in = ortho_out.readline()

max_taxa = max(taxa_numbers)


ortho_out = open(orthomcl_out, "rU")
inl = ortho_out.readline()


align_protlist = []
align_nucllist = []
	
counter = 0
while inl != "":
	
	orthogroup = inl.split("(")[0]
	num_taxa = int(inl.split(",")[1].split(" ")[0])
	geneslist = inl.split(" ")[3:]
	
	
	species_list = []

	
	if (num_taxa == max_taxa) and (num_taxa == len(geneslist)):

		for gene in geneslist:
                        species = gene.split("(")[1].split(")")[0]
			ID = gene.split("(")[0]
			species_list.append((species,ID))
			
		protlist = fasta_extract(fasta_prot_dict,species_list)
		SeqIO.write(protlist,"temp.fa","fasta")
		
		if run_tcoffee == 1:
			nucllist = fasta_extract(fasta_nucl_dict,species_list)
			SeqIO.write(nucllist, "tempnucl.fa","fasta")
		else:
			pass
		
		os.system("muscle -in temp.fa -out temp.clw -clwstrict -quiet")
		os.system("cp temp.clw prot_alignments/%s"%orthogroup)
		if run_tcoffee == 1:
			os.system("t_coffee -other_pg seq_reformat -in tempnucl.fa -in2 temp.clw -action +thread_dna_on_prot_aln -output clustalw > temp.clustal")
			alignnucl= AlignIO.read("temp.clustal","clustal")
			os.system("mv temp.clustal alignments/%s.aln"%orthogroup)
			alignnucl.sort()
			align_nucllist.append(alignnucl)
		else:
			pass
		alignprot= AlignIO.read("temp.clw","clustal")
		alignprot.sort()
		align_protlist.append(alignprot)
		
		
		
				
	else:
		pass
	inl = ortho_out.readline()
	
	if ((counter*100)/len(taxa_numbers)) in range(1,100):
		print ("progress: ", counter*100/len(taxa_numbers), "%", end="\r") 
		
	elif ((counter*100)/len(taxa_numbers)) >= 99:
		print ("complete!")
	else:
		pass
	counter = counter +1
		
#create directory for raxml output
if os.path.exists("raxml_out") == False:
	os.system("mkdir raxml_out")
else:
	pass

