import matplotlib  # must import first
matplotlib.use('Agg')  # allows you to not have an x-server running
import os
import sys, os.path
import matplotlib.pyplot as plt
import pylab
import numpy as np

#####################################
# written by Reed Stein
# 10/2017 - 4/2018
#
#######################################


def convert_to_perc(property_list):

	y,binEdges=np.histogram(property_list,bins=10)

	perc_prop_list = []
	for term in y:
		new_term = float(term) / float(len(property_list))
		perc_prop_list.append(new_term)

	return(perc_prop_list, binEdges)


def my_plot(property_list_for_lig, property_list_for_dec, system_name, property_name, property_type):

	lig_perc_list, lig_binEdges = convert_to_perc(property_list_for_lig)
	dec_perc_list, dec_binEdges = convert_to_perc(property_list_for_dec)
	
	bincenters1 = 0.5*(lig_binEdges[1:]+lig_binEdges[:-1])
	bincenters2 = 0.5*(dec_binEdges[1:]+dec_binEdges[:-1])
	plt.margins(0.2)
	plt.plot(bincenters1,lig_perc_list,'b--', label='Ligands', lw=2)
	plt.plot(bincenters2,dec_perc_list,'r:', label='Decoys', lw=2)
	plt.ylabel("Percentage")
	plt.legend(loc='upper left')
	pylab.title(system_name+" "+property_name)
	plt.xlabel(property_name)
	#plt.show()
	#sys.exit()
	
	pwd = os.getcwd()+"/"
	pylab.savefig(pwd+system_name+"_new_decoys_"+property_type+".png", dpi=600)
	plt.clf()

def collect_properties(file_list, property_type):

	#SMILES ID MW LogP RotB HBD HBA Q
	#SMILES ID MW LogP RotB HBD HBA Q PROT_ID

	if property_type.lower()[0] == "m":
		lig_num = 3
		dec_num = 4
		property_name = "MOLECULAR WEIGHT (Da)"
	elif property_type.lower()[0] == "l":
		lig_num = 4
		dec_num = 5
		property_name = "LogP"
	elif property_type.lower()[0] == "r":
		lig_num = 5
		dec_num = 6
		property_name = "# Rotatable Bonds"
	elif property_type.lower() == "hbd":
		lig_num = 6
		dec_num = 7
		property_name = "# Hydrogen Bond Donors"
	elif property_type.lower() == "hba":
		lig_num = 7
		dec_num = 8
		property_name = "# Hydrogen Bond Acceptors"
	elif property_type.lower()[0] == "c" or property_type.lower()[0] == "q":
		lig_num = 8 
		dec_num = 9
		property_name = "Net Charge"
	

	prop_lig_list = []
	prop_dec_list = []

	for txt in file_list:
		decoy_count = 0
		open_txt = open(txt, 'r')
		read_txt = open_txt.readlines()
		open_txt.close()
		
		for line in read_txt:
			splitline = line.strip().split()
			if len(splitline) > 0:
				if splitline[0] == "LIGAND:":
					lig_ID = splitline[2]
					#mw = int(round(float(splitline[3])))
					lig_prop = int(round(float(splitline[lig_num])))
					prop_lig_list.append(lig_prop)
					

				elif splitline[0] == "DECOY":
					decoy_count += 1
					dec_ID = splitline[3]
					dec_prop = int(round(float(splitline[dec_num])))
					prop_dec_list.append(dec_prop)

	return(prop_lig_list, prop_dec_list, property_name)

def main():

	pwd = os.getcwd()+"/"
	direc_name = sys.argv[1]
	if direc_name[-1] == "/":
		direc_name = direc_name[0:-1]

	smiles_dir = pwd+direc_name+"/"

	os.chdir(smiles_dir)

	suffix_len = len("_final_property_matched_decoys.txt")
	
	decoy_file_list = [name for name in os.listdir(".") if (os.path.isfile(name) and name[-suffix_len:] == "_final_property_matched_decoys.txt")]

	full_txt_list = []
	for decoy in decoy_file_list:
		full_txt_list.append(smiles_dir+decoy)

	#os.chdir(pwd)

	property_list = ["mw", "logp", "rb", "hbd", "hba", "q"]
	for prop in property_list:
		prop_list_for_lig, prop_list_for_dec, property_name = collect_properties(full_txt_list, prop)
		my_plot(prop_list_for_lig, prop_list_for_dec, direc_name, property_name, prop)


main()
