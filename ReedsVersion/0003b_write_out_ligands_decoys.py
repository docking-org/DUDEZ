import os, sys

###################################
# written by Reed Stein
# 10/2017 - 4/2018
#
####################################

def write_out_lig_decoys(lig_list, dec_list):

	output = open("ligands.smi", 'w')
	for lig in lig_list:
		name = lig[0]
		smiles = lig[1]
		if "_" in name:
			name = name.split("_")[0]
		output.write(smiles+" "+name+"\n")

	output.close()
	
	output1 = open("decoys.smi", 'w')
	output2 = open("decoy_protomers.smi", 'w')
	for dec in dec_list:
		name = dec[0]
		smiles = dec[1]
		prot_smiles = dec[2]
		if "_" in name:
			name = name.split("_")[0]
		output1.write(smiles+" "+name+"\n")
		output2.write(prot_smiles+" "+name+"\n")

	output1.close()
	output2.close()
	
	pwd = os.getcwd()+"/"
	decoy_dir = pwd+"decoys/"
	if not os.path.isdir(decoy_dir):
		os.system("mkdir "+decoy_dir)

	os.chdir(decoy_dir)
	os.system("cp ../decoys.smi .")
	os.system("cp ../decoy_protomers.smi .")


def check_int(prot_id):

	try:
		int(prot_id)
		return(True)

	except ValueError:
		return(False)
		
def main():

	pwd = os.getcwd()+"/"
	if len(sys.argv) != 3:
		print("Syntax: python 0004b_write_ligands_decoys.py smiles_dir new_dir_name")
		sys.exit()

	system_dir = sys.argv[1]+"/"
	new_dir_name = sys.argv[2]

	if not os.path.isdir(system_dir):
		print(system_dir+" does not exist")
		sys.exit()

	decoy_files = [name for name in os.listdir(system_dir) if (os.path.isfile(system_dir+name) and name.endswith("_final_property_matched_decoys.txt"))]

	lig_list = []
	dec_list = []
	repeat_list = []
	for decoy_file in decoy_files:
		decoy_file_name = system_dir+decoy_file
		if os.path.isfile(decoy_file_name):

			open_fin = open(decoy_file_name, 'r')
			read_fin = open_fin.readlines()
			open_fin.close()
			
			for line in read_fin:
				splitline = line.strip().split()
				if len(splitline) > 0:
					if splitline[0] == "LIGAND:":
						lig_smiles = splitline[1]
						lig_name = splitline[2]
						
						if lig_name not in repeat_list:
							repeat_list.append(lig_name)
							lig_list.append([lig_name, lig_smiles])

					if splitline[0] == "DECOY":
						dec_smiles = splitline[2]
						dec_name = splitline[3]
						dec_prot_id = splitline[10]
						
						int_or_not = check_int(dec_prot_id)
						
						if int_or_not == False: ### make sure it does not have a protomer ID
							if dec_name not in repeat_list:
								repeat_list.append(dec_name)
								dec_list.append([dec_name, dec_smiles, dec_prot_id])

	full_new_dir_path = pwd+new_dir_name+"/"
	if not os.path.isdir(full_new_dir_path):
		os.system("mkdir "+full_new_dir_path)
	os.chdir(full_new_dir_path)
						
	write_out_lig_decoys(sorted(lig_list), sorted(dec_list))


main()
