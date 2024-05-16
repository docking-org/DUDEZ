from __future__ import print_function, absolute_import
import os, sys
from rdkit import Chem as C
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD

def get_stuff(smiles):

	mol = C.MolFromSmiles(smiles)
	#hac = D.HeavyAtomCount(mol)
	if len(str(mol)) == 0:
		return(0, 0, 0, 0, 0, 0)

	else:
		mw = CD.CalcExactMolWt(mol)
		logp = C.Crippen.MolLogP(mol)
		rotB = D.NumRotatableBonds(mol)
		HBA = CD.CalcNumHBA(mol)
		HBD = CD.CalcNumHBD(mol)
		q = C.GetFormalCharge(mol)
	return(int(mw), int(logp), rotB, HBA, HBD, q)


def compare(list1, list2):

	true_count = 0
	for i in range(len(list1)):
		val1 = list1[i]
		val2 = list2[i]
		if val1 == val2:
			true_count += 1

	if true_count == 6:
		return(True)
	else:
		return(False)
			


def make_uniq_prot(prot_dict):

	repeat_list = []
	uniq_prot_dict = {}
	for p in prot_dict:
		if len(prot_dict[p]) == 1:
			uniq_prot_dict[p+"_0"] = [prot_dict[p][0][6]] # set SMILES
		else:
			for i in range(len(prot_dict[p])-1):
				f_mw = prot_dict[p][i][0]
				f_logp = prot_dict[p][i][1]
				f_rotB = prot_dict[p][i][2]
				f_HBA = prot_dict[p][i][3]
				f_HBD = prot_dict[p][i][4]
				f_q = prot_dict[p][i][5]
				f_smiles = prot_dict[p][i][6]
				true1 = prot_dict[p][i][7]
				list1 = [f_mw, f_logp, f_rotB, f_HBA, f_HBD, f_q]
				if true1 == True:
					lig_ID_count = repeat_list.count(p)
					dict_entry = p+"_"+str(lig_ID_count)
					repeat_list.append(p)				
					uniq_prot_dict[dict_entry] = [f_smiles]
	
					for j in range(i+1, len(prot_dict[p])):
						s_mw = prot_dict[p][j][0]
						s_logp = prot_dict[p][j][1]
						s_rotB = prot_dict[p][j][2]
						s_HBA = prot_dict[p][j][3]
						s_HBD = prot_dict[p][j][4]
						s_q = prot_dict[p][j][5]
						s_smiles = prot_dict[p][j][6]
						true2 = prot_dict[p][j][7]
						if true2 == True:
							list2 = [s_mw, s_logp, s_rotB, s_HBA, s_HBD, s_q]
							same_or_not = compare(list1, list2)
							if same_or_not == True:
								prot_dict[p][j][7] = False # set to False if it is identical
							else:
								lig_ID_count = repeat_list.count(p)
								dict_entry = p+"_"+str(lig_ID_count)
								repeat_list.append(p)

								uniq_prot_dict[dict_entry] = [s_smiles] 


	print("LENGTH OF UNIQUE PROTOMERS = ",len(uniq_prot_dict))
	return(uniq_prot_dict)


def protonate_ligands(path, prot_dict):

	prot_dir = path+"test_ligand_protonation/"
	if not os.path.isdir(prot_dir):
		os.system("mkdir "+prot_dir)
	os.chdir(prot_dir)

	to_prot_file = prot_dir+"test_ligand_ids.smi"
	output = open(to_prot_file, 'w')
	for zinc_id in prot_dict:
		smiles = prot_dict[zinc_id][0]
		output.write(smiles+" "+zinc_id+"\n")
	output.close()

	os.system("source /nfs/soft/python/current/env.csh")
	os.system("source /nfs/home/rstein/.cshrc_dbgen_corina")
	os.system("source /nfs/soft/jchem/current/env.csh")
	os.system("bash /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/generate_protomers_from_smiles.sh -H 7.4 "+to_prot_file)

	built_dir = prot_dir+"/working/protonate/"
        #built_file = built_dir+"test_ligand_ids-protomers.ism"
	built_file = built_dir+"test_ligand_ids-protomers-expanded.ism"

	new_prot_dict = {}
	repeat_list = []
	if os.path.isfile(built_file):
		open_built = open(built_file, 'r')
		read_built = open_built.readlines()
		open_built.close()

		for line in read_built:
			splitline = line.strip().split()
			smiles = splitline[0]
			lig_ID = splitline[-2]
			
			mw, logp, rotB, HBA, HBD,q = get_stuff(smiles)
			if lig_ID not in new_prot_dict:
				new_prot_dict[lig_ID] = [[mw, logp, rotB, HBA, HBD, q, smiles, True]]
			elif lig_ID in new_prot_dict:
				new_prot_dict[lig_ID].append([mw, logp, rotB, HBA, HBD, q, smiles, True])

		uniq_prot_dict = make_uniq_prot(new_prot_dict)
				

	os.chdir(path)
	os.system("rm -rf "+prot_dir)

	return(uniq_prot_dict)

def gather_ligands(ligand_file):

	open_lig = open(ligand_file,'r')
	read_lig = open_lig.readlines()
	open_lig.close()

	#lig_list = []
	lig_dict = {}
	for line in read_lig:
		splitline = line.strip().split()
		lig_ID = splitline[-1]
		smiles = splitline[0]
		lig_dict[lig_ID] = [smiles]
		#lig_num = lig_list.count(lig_ID)
		#key_entry = lig_ID+"_"+str(lig_num)
		#lig_list.append(lig_ID)
		#lig_dict[key_entry] = [smiles]

	return(lig_dict)

def create_lig_dirs(smiles_dir, prot_dict, infile, smiles_file):
	
	### first calculate fingerprints for all ligands to be copied into all directories
	insmiles = "all_ligands.smi"
	os.system("cp "+smiles_file+" "+insmiles)
	os.system("python /mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/teb_chemaxon_cheminf_tools/generate_chemaxon_fingerprints.py "+insmiles+" all_ligands")
	os.system("/mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/convert_fp_2_fp_in_16unit/convert_fp_2_fp_in_uint16 all_ligands.fp "+insmiles+" all_ligands")

	lig_uint16_fp = "all_ligands_uint16.fp"
	lig_uint16_count = "all_ligands_uint16.count"

	if not (os.path.isfile(lig_uint16_fp) and os.path.isfile(lig_uint16_count)):
		print("fingerprints for ligands failed")
		sys.exit()

	lig_count = 0
	map_output = open(smiles_dir+"LIGAND_MAP.txt", 'w') ### write out a map of ligand names
	for prot in prot_dict:
		smiles = prot_dict[prot][0]
		lig_count += 1
		
		print(smiles, prot)
		lig_dir = smiles_dir+"/ligand_"+str(lig_count)+"/"
		if not os.path.isdir(lig_dir):
			os.system("mkdir "+lig_dir)
			
			os.chdir(lig_dir)
			os.system("cp "+infile+" .")
			output = open("ligand_"+str(lig_count)+".smi", 'w')
			output.write(smiles+" "+prot+"\n")
			output.close()
			map_output.write("ligand_"+str(lig_count)+" "+smiles+" "+prot+"\n")
			#lig_set_name = "lig_set_"+str(prot)+".smi"
			#output = open(lig_set_name, 'w')
			#output.write(smiles+" "+prot+"\n")
			#output.close()
			
			#os.system("cp "+lig_set_name+" ligand_"+str(lig_count)+".smi")
			
			os.chdir(smiles_dir)

	map_output.close()

def read_infile(infile):

	open_in = open(infile, 'r')
	read_in = open_in.readlines()
	open_in.close()

	for line in read_in:
		splitline = line.strip().split()
		if len(splitline) > 0:
			if splitline[0].lower() == "protonate":
				tf = splitline[-1]

	return(tf)

def main():

	
	pwd = os.getcwd()+"/"
	if len(sys.argv) != 3:
		print("Syntax: python 0000_protonate_setup_dirs.py smiles_file.smi directory_name")
		sys.exit()

	infile = pwd+"decoy_generation.in"
	protonate_ligs = False
	if not os.path.isfile(infile):
		print("decoy_generation.in does not exist")
		sys.exit()
	else:
		tf = read_infile(infile)
		if tf.lower()[0] == "y":
			protonate_ligs = True
		
	smiles_file = pwd+sys.argv[1]
	smiles_dir = pwd+sys.argv[2]+"/"

	if not os.path.isdir(smiles_dir):
		os.system("mkdir "+smiles_dir)

	os.chdir(smiles_dir)
	os.system("cp "+infile+" .")
	os.system("cp "+smiles_file+" all_ligands.smi")
		
	lig_dict = gather_ligands(smiles_file)
	if protonate_ligs == True:
		prot_dict = protonate_ligands(smiles_dir, lig_dict) 
		if len(prot_dict) == 0:
			print("No protomers generated. Make sure to run the following commands:\n")
			print("source /nfs/soft/python/current/env.csh")
			print("source /nfs/home/rstein/.cshrc_dbgen_corina")
			print("source /nfs/soft/jchem/current/env.csh")
			sys.exit()

		create_lig_dirs(smiles_dir, prot_dict, infile, smiles_file)
	else:
		create_lig_dirs(smiles_dir, lig_dict, infile, smiles_file)


main()
