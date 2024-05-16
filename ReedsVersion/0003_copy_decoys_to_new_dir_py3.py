import os, sys
import gzip
import subprocess

###################################
# written by Reed Stein
# 10/2017 - 4/2018
#
####################################

def count_gz_lines(gz_file):

	protomer_open = subprocess.Popen(['gzip', '-cdfq', gz_file], stdout=subprocess.PIPE)
	count = 0
	for line in protomer_open.stdout:
		splt = line.strip().split()
		if len(splt) > 0:
			if splt[0].decode('utf-8') == "M": ## has db2.gz file format
				count += 1
				break
		count += 1

	return(count)

def find_path(prot_id):

        #print(prot_id)
	lp = len(prot_id)
	if lp < 7:
		print("Error")
		exit()
	dir1 = prot_id[lp-2:lp]
	dir2 = prot_id[lp-4:lp-2]
	dir3 = prot_id[lp-6:lp-4]
	#print(dir1, dir2, dir3)
	
	# /nfs/dbraw/zinc/49/76/33/419497633.db2.gz
	db2file = "/nfs/dbraw/zinc/"+dir3+"/"+dir2+"/"+dir1+"/"+prot_id+".db2.gz"
	#mol2file = "/nfs/dbraw/zinc/"+dir3+"/"+dir2+"/"+dir1+"/"+prot_id+".mol2.gz"
	if not os.path.isfile(db2file):
		print(db2file+" does not exist")

	if os.path.isfile(db2file):
		try:
			num_db2_lines = count_gz_lines(db2file)

		except IOError:
			print(prot_id+" is empty")

		if num_db2_lines > 0:
			return(db2file)
			
	return("0")


        #if not os.path.isfile(mol2file):
        #        print(mol2file, prot_id+" has not been built")
        #        return("0")
        #if os.path.isfile(mol2file):
        #        try:
        #                num_mol2_lines = count_gz_lines(mol2file)
        #        except IOError:
        #                return("0")
        #        if num_mol2_lines > 0:
        #                return(mol2file)

 #       return("0")

def write_fil(outlist, outname):

	output = open(outname, 'w')
	for l in outlist:
		name = l[0]
		smi = l[1]
		if "_" in name:
			name = name.split("_")[0]
		output.write(smi+" "+name+"\n")
	output.close()

def write_out_lig_decoys(lig_list, dec_list, decoy_build_list):

	write_fil(lig_list, "ligands.smi")

	write_fil(dec_list, "decoys.smi")

	if len(decoy_build_list) > 0:
		write_fil(decoy_build_list, "decoys_to_build.smi")

	pwd = os.getcwd()+"/"
	decoy_dir = pwd+"decoys/"
	if not os.path.isdir(decoy_dir):
		os.system("mkdir "+decoy_dir)

	os.system("cp decoys.smi "+decoy_dir)
	os.chdir(decoy_dir)

	for dec in dec_list:
		dec_id_name = dec[0]
		prot_id_path = dec[2]
		os.system("cp "+prot_id_path+" "+dec_id_name+".db2.gz")


def check_int(prot_id):

	try:
		int(prot_id)
		return(True)

	except ValueError:
		return(False)
		
def main():

	pwd = os.getcwd()+"/"
	if len(sys.argv) != 3:
		print("Syntax: python 0004_copy_decoys_to_new_dir.py smiles_dir new_dir_name")
		sys.exit()

	system_dir = sys.argv[1]+"/"
	new_dir_name = sys.argv[2]

	if not os.path.isdir(system_dir):
		print(system_dir+" does not exist")
		sys.exit()

	decoy_files = [name for name in os.listdir(system_dir) if (os.path.isfile(system_dir+name) and name.endswith("_final_property_matched_decoys.txt"))]
	
	lig_list = []
	dec_list = []
	decoy_build_list = []
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
						if int_or_not == False:
							if dec_name not in repeat_list:
								repeat_list.append(dec_name)
								decoy_build_list.append([dec_name, dec_prot_id])
							#print(dec_name+" does not have a protomer ID. You will need to build this yourself.")

						if dec_name not in repeat_list:
							repeat_list.append(dec_name)
							path = find_path(dec_prot_id)
							#print(dec_prot_id, path)
							if path == "0":
								print(dec_prot_id+" does not have a db2 file")
							else:
								dec_list.append([dec_name, dec_smiles, path])

	full_new_dir_path = pwd+new_dir_name+"/"
	if not os.path.isdir(full_new_dir_path):
		os.system("mkdir "+full_new_dir_path)
	os.chdir(full_new_dir_path)
						
	write_out_lig_decoys(sorted(lig_list), sorted(dec_list), sorted(decoy_build_list))

		

main()
