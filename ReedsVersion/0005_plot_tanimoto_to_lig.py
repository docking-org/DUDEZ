import matplotlib  # must import first
matplotlib.use('Agg')  # allows you to not have an x-server running
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pylab

###################################
# written by Reed Stein
# 10/2017 - 4/2018
#
#
####################################


def boxplot(data, namelist, system):

	#data = [es_list1, ld_list1, vdw_list1, es_list3, ld_list3, vdw_list3]
	#lbl = ["ES_up", "LD_up", "vdW_up", "ES_down", "LD_down", "vdW_down"]
	#plt.xticks(list('123456'), ["ES_up", "LD_up", "vdW_up", "ES_down", "LD_down", "vdW_down"])
	
	
	position_list = []
	xtick_list = []
	for i in range(1, len(namelist)+1):
		xtick_list.append(i)
		position_list.append(float(i))

	plt.xticks(xtick_list, namelist, rotation='vertical', fontsize=8)
	fig = plt.figure(1, figsize=(25,10))
	ax = fig.add_subplot(111)
	bp = ax.boxplot(data, patch_artist=True)
	#bp = ax.boxplot(data, showfliers=False)
	
	for box in bp['boxes']:
		box.set(color='#ffa07a', linewidth=2)

	for flier in bp['fliers']:
		#flier.set(marker='-', color='#e7298a', alpha=0.2)
		flier.set(marker='o', color='w',alpha=0.2)

	for cap in bp['caps']:
		cap.set(color='#7570b3', linewidth=4)


	#pylab.legend([bp["boxes"][0], bp["boxes"][1]], ["ES", "LD"], loc='upper right', prop={'size':6})
	pylab.ylabel("Maximum Matched Decoy Tc to All Ligands",fontsize=18)
	pylab.xlabel("Ligand",fontsize=18)
	pylab.title(system+" Matched Decoy Similarity Comparison" ,fontsize=15)
	#pylab.figure(figsize=(20,20))
	pylab.subplots_adjust(bottom=0.15)
	
	#plt.figtext(0.80, 0.08, "vdW", backgroundcolor='#ffe4e1')
	
	#plt.show()
	pwd = os.getcwd()+"/"
	pylab.savefig(pwd+"/"+system+"_tanimoto_similarity.png", dpi=600)


def collect_IDs(smiles_dir):

	decoy_file_list = [name for name in os.listdir(smiles_dir) if (os.path.isfile(smiles_dir+name) and name.endswith("_final_property_matched_decoys.txt"))]
	
	lig_dict = {}
	decoy_dict = {}
	name_list = []
	tanimoto_dict = {}
	for decoy_file in decoy_file_list:
		print("COLLECTING DECOYS OF "+decoy_file.split("_final")[0])
		open_decoy = open(decoy_file, 'r')
		read_decoy = open_decoy.readlines()
		open_decoy.close()
		
		for line in read_decoy:
			splitline = line.strip().split()
			if len(splitline) > 0:
				if splitline[0] == "LIGAND:":
					lig_ID = splitline[2]
					if lig_ID not in lig_dict:	
						lig_dict[lig_ID] = []
						name_list.append(lig_ID)
					if lig_ID not in tanimoto_dict:
						tanimoto_dict[lig_ID] = []

				if splitline[0] == "DECOY":
					dec_ID = splitline[3]
					tc_to_any_lig = float(splitline[11]) 
					if dec_ID not in decoy_dict:
						decoy_dict[dec_ID] = lig_ID
					if dec_ID not in lig_dict[lig_ID]:
						lig_dict[lig_ID].append(dec_ID)
						tanimoto_dict[lig_ID].append(tc_to_any_lig)

	name_list = sorted(name_list)
	data = []
	for name in name_list:
		print(name, tanimoto_dict[name])
		print(np.mean(tanimoto_dict[name]))
		data.append(tanimoto_dict[name])

	return(name_list, data)


def main():

	pwd = os.getcwd()+"/"
	
	smiles_dir = pwd+sys.argv[1]+"/"
	os.chdir(smiles_dir)
	
	name_list, data = collect_IDs(smiles_dir)
	
	if sys.argv[1].endswith("/"):
		plot_name = sys.argv[1].split("/")[0]
	boxplot(data, name_list, plot_name)

main()
