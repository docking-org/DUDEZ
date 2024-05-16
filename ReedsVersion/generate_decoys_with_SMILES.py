from __future__ import print_function, absolute_import
import os, sys, random
from rdkit import Chem as C
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD
from datetime import datetime

def get_stuff(smiles, q_or_not):

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
		if q_or_not[0] == "y":
			q = C.GetFormalCharge(mol)
			return(mw, logp, rotB, HBD, HBA, q)
		else:
			return(mw, logp, rotB, HBD, HBA)


def compare(mw, logp, rotB, HBD, HBA, q, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, prop_range_dict):

	#CHG_RANGES =  [  2]
	#HBD_RANGES =  [  3]
	#HBA_RANGES =  [  4]
	#RB_RANGES =   [  5]
	#MWT_RANGES =  [125]
	#LOGP_RANGES = [3.6]
	
	CHG_RANGES = prop_range_dict["CHG_RANGES"]
	HBD_RANGES = prop_range_dict["HBD_RANGES"]
	HBA_RANGES = prop_range_dict["HBA_RANGES"]
	RB_RANGES = prop_range_dict["RB_RANGES"]
	MWT_RANGES = prop_range_dict["MWT_RANGES"]
	LOGP_RANGES = prop_range_dict["LOGP_RANGES"]
	
	#print(prop_range_dict)
	
	mwt_round = int(round(mw))
	mwt_low = mwt_round-MWT_RANGES[0]
	mwt_high = mwt_round+MWT_RANGES[0]
	
	logp_round = round(logp, 2)
	logp_low = logp_round-LOGP_RANGES[0]
	logp_high = logp_round+LOGP_RANGES[0]
	
	rotB_low = rotB-RB_RANGES[0]
	rotB_high = rotB+RB_RANGES[0]
	
	HBD_low = HBD-HBD_RANGES[0]
	HBD_high = HBD+HBD_RANGES[0]
	
	HBA_low = HBA-HBA_RANGES[0]
	HBA_high = HBA+HBA_RANGES[0]
	
	CHG_low = int(q)-CHG_RANGES[0]
	CHG_high = int(q)+CHG_RANGES[0]
	
	#print(mw, logp, rotB, HBD, HBA, q)
	#print(dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q)
	if dec_mw >= mwt_low and dec_mw <= mwt_high:

		if dec_logp >= logp_low and dec_logp <= logp_high:

			if dec_rotB >= rotB_low and dec_rotB <= rotB_high:

				if dec_HBD >= HBD_low and dec_HBD <= HBD_high:

					if dec_HBA >= HBA_low and dec_HBA <= HBA_high:

						if dec_q != "NA": ### only compare charge when the SMILES have been protonated
							if dec_q >= CHG_low and dec_q <= CHG_high: 
								return(True)
						else:
							return(True)

	return(False) 

def map_tranche(mwt, logp):

	if logp <= -1:
		logp_tranche = 'A'
	elif logp > -1 and logp <= 0:
		logp_tranche = 'B'
	elif logp > 0 and logp <= 1:
		logp_tranche = 'C'
	elif logp > 1 and logp <= 2:
		logp_tranche = 'D'
	elif logp > 2 and logp <= 2.5:
		logp_tranche = 'E'
	elif logp > 2.5 and logp <= 3:
		logp_tranche = 'F'
	elif logp > 3 and logp <= 3.5:
		logp_tranche = 'G'
	elif logp > 3.5 and logp <= 4:
		logp_tranche = 'H'
	elif logp > 4 and logp <= 4.5:
		logp_tranche = 'I'
	elif logp > 4.5 and logp <= 5:
		logp_tranche = 'J'
	elif logp > 5:
		logp_tranche = 'K'
	
	
	if mwt <= 200:
		mwt_tranche = 'A'
	elif mwt > 200 and mwt <= 250:
		mwt_tranche = 'B'
	elif mwt > 250 and mwt <= 300:
		mwt_tranche = 'C'
	elif mwt > 300 and mwt <= 325:
		mwt_tranche = 'D'
	elif mwt > 325 and mwt <= 350:
		mwt_tranche = 'E'
	elif mwt > 350 and mwt <= 375:
		mwt_tranche = 'F'
	elif mwt > 375 and mwt <= 400:
		mwt_tranche = 'G'
	elif mwt > 400 and mwt <= 425:
		mwt_tranche = 'H'
	elif mwt > 425 and mwt <= 450:
		mwt_tranche = 'I'
	elif mwt > 450 and mwt <= 500:
		mwt_tranche = 'J'
	elif mwt > 500:
		mwt_tranche = 'K'
	
	tranche = mwt_tranche+logp_tranche

	return(tranche)
	

def extend_tranche(mwt_or_logp):

	if mwt_or_logp == 'A':
		ex_tranche = '[AB]' 
	if mwt_or_logp == 'B':
		ex_tranche = '[ABC]' 
	if mwt_or_logp == 'C':
		ex_tranche = '[BCD]' 
	if mwt_or_logp == 'D':
		ex_tranche = '[CDE]' 
	if mwt_or_logp == 'E':
		ex_tranche = '[DEF]' 
	if mwt_or_logp == 'F':
		ex_tranche = '[EFG]' 
	if mwt_or_logp == 'G':
		ex_tranche = '[FGH]' 
	if mwt_or_logp == 'H':
		ex_tranche = '[GHI]' 
	if mwt_or_logp == 'I':
		ex_tranche = '[HIJ]' 
	if mwt_or_logp == 'J':
		ex_tranche = '[IJK]' 
	if mwt_or_logp == 'K':
		ex_tranche = '[JK]' 

	return(ex_tranche)
	
def check_len(smi_file):

#	time.sleep(10)
	num_lines = sum(1 for line in open(smi_file, 'r'))
	return(num_lines)

def zinc_subfunc(smi_file, decoy_count, decoy_dict, smiles_list, compound_list, mwt, logp, rotB, hbd, hba, q, prop_range_dict, generate_num):

	print(smi_file, check_len(smi_file), decoy_count)
	open_smi = open(smi_file, 'r')
	read_smi = open_smi.readlines()
	open_smi.close()
	line_count = 0
	#for line in open_smi:
	for line in read_smi:
		line_count += 1
		if line_count == 1: ### skip the header
			continue
		pick_or_not = random.randint(0,1)
		if pick_or_not == 0:
			continue
		else:
			splt = line.strip().split()
			if len(splt) > 1:
				dec_smiles = splt[0]
				dec_zinc_id = splt[1]
				if "ZINC" in dec_zinc_id:
					dec_mwt, dec_logp, dec_rotB, dec_HBD, dec_HBA = get_stuff(dec_smiles, 'n')
					boolliee = compare(mwt, logp, rotB, hbd, hba, q, dec_mwt, dec_logp, dec_rotB, dec_HBD, dec_HBA, "NA", prop_range_dict)
					if boolliee == True:
						if dec_zinc_id not in decoy_dict:
							#print(dec_smiles, dec_zinc_id)
							decoy_count += 1
							decoy_dict[dec_zinc_id] = [dec_smiles, dec_mwt, dec_logp, dec_rotB, dec_HBD, dec_HBA, "NA"] ### no prot_ids from SMILES so put "NA"
							smiles_list.append(dec_smiles)
							compound_list.append(dec_zinc_id)
							if decoy_count >= generate_num:
								print("DECOYS FOUND: ",decoy_count)
								break
	print("DECOYS FOUND: ",decoy_count)
	
	#open_smi.close()

	return(decoy_count, decoy_dict, smiles_list, compound_list)

def check_time(start_time):

	new_time = float(str(datetime.now() - start_time).split(":")[-1])
	return(new_time)

def query_zinc(lig_id, smiles, pwd, prop_range_dict, ligand_file):

	generate_num = prop_range_dict["GENERATE_NUM"]
	LIG_TC_RANGE = prop_range_dict["LIGAND_TC_RANGE"]
	
	mwt, logp, rotB, hbd, hba, q = get_stuff(smiles, 'y')
	lig_list = [smiles, lig_id, mwt, logp, rotB, hbd, hba, q]
	
	#MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
	#LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]
	
	MWT_RANGES = prop_range_dict["MWT_RANGES"]
	LOGP_RANGES = prop_range_dict["LOGP_RANGES"]
	
	mwt_round = int(round(mwt))
	logp_round = round(logp, 2)
	
	### only pick the farthest range
	mwt_low = mwt-MWT_RANGES[0]
	mwt_high = mwt+MWT_RANGES[0]
	
	logp_low = logp_round-LOGP_RANGES[0]
	logp_high = logp_round+LOGP_RANGES[0]
	
	print("MWT RANGE for "+str(mwt)+":", mwt_low, mwt_high)
	print("LogP RANGE for "+str(logp)+":", logp_low, logp_high)
	
	tranche = map_tranche(mwt, logp)
	
	### only consider most specific tranche first
	twoD_dir = "/mnt/nfs/ex3/published/2D/"+tranche+"/"
	smi_list = [name for name in os.listdir(twoD_dir) if (os.path.isfile(twoD_dir+name) and name[-4:] == ".smi")]
	random.shuffle(smi_list) ### randomly shuffle the SMILES file list to minimize redundancy of decoys retrieved
	
	decoy_count = 0
	decoy_dict = {}
	smiles_list = [] ### not including ligand SMILES anymore
	compound_list = [] ### not including ligand ID anymore

	#start_time = datetime.now()
	for smi in smi_list:
		if decoy_count >= generate_num:
			break
		else:
			decoy_count, decoy_dict, smiles_list, compound_list = zinc_subfunc(twoD_dir+smi, decoy_count, decoy_dict, smiles_list, compound_list, mwt, logp, rotB, hbd, hba, q, prop_range_dict, generate_num)

	if decoy_count <= 50:
		### This should not happen when querying SMILES
		print("NOT ENOUGH DECOYS. CONSIDER WIDENING PROPERTY RANGE")
		sys.exit()

	print("Running Tanimoto")
	print("LENGTH OF SMILES_LIST = ",len(smiles_list))
	print("LENGTH OF COMPOUND_LIST = ",len(compound_list))
	ori_lig_props = [mwt, logp, rotB, hbd, hba, q]
	map_dict = run_tanimoto(smiles_list, compound_list, ligand_file, LIG_TC_RANGE) 
	
	filtered_decoy_list = remove_wrong_protomers(ori_lig_props, map_dict, pwd, prop_range_dict)
	
	if len(filtered_decoy_list) > 0:
		write_out_approved_decoys(lig_list, filtered_decoy_list)

def run_tanimoto(smiles_list, compound_list, ligand_file, lig_tc_list):

        ### This new run_tanimoto function was added 8/14/2019.
        ### Should speed up the calculation as only the ligand
        ### is being compared to the decoys, not everything vs everything
        ### Also now running Tanimoto calculation before protomer generation
        ### step to avoid protonating molecules that are too similar to ligands

	print("RUNNING TANIMOTO CALCULATION")
	
	map_dict = {}
	output = open("test_decoy_smiles.smi", 'w')
	for i in range(len(smiles_list)):
		comp_name = compound_list[i]
		comp_smiles = smiles_list[i]
		output.write(comp_smiles+" "+comp_name+"\n")
		map_dict[comp_name] = comp_smiles
	output.close()

	### these were already calculated while setting up the directories
	tot_ligand_file = "../all_ligands.smi"
	lig_uint16_fp = "../all_ligands_uint16.fp"
	lig_uint16_count = "../all_ligands_uint16.count"
	
	### calculate fingerprints of property matched decoys
	dec_file_pref = "test_decoy_smiles"
	decoy_file = "test_decoy_smiles.smi"
	#os.system("python /mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/teb_chemaxon_cheminf_tools/generate_chemaxon_fingerprints.py "+decoy_file+" "+dec_file_pref)
	os.system("python /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/generate_chemaxon_fingerprints.py "+decoy_file+" "+dec_file_pref)
	os.system("/mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/convert_fp_2_fp_in_16unit/convert_fp_2_fp_in_uint16 test_decoy_smiles.fp "+decoy_file+" "+dec_file_pref)
	
	dec_uint16_fp = dec_file_pref+"_uint16.fp"
	dec_uint16_count = dec_file_pref+"_uint16.count"
	
	if not (os.path.isfile(dec_uint16_fp) and os.path.isfile(dec_uint16_count) and os.path.isfile(lig_uint16_fp) and os.path.isfile(lig_uint16_count)):
		print("ERROR! Something went wrong with fingerprinting ligands and decoys")
		sys.exit()
	else:
		os.system("/mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/cal_Tc_matrix_uint16/cal_Tc_matrix_uint16 "+dec_uint16_fp+" "+decoy_file+" "+dec_uint16_count+" "+lig_uint16_fp+" "+tot_ligand_file+" "+lig_uint16_count+" Max_Tc_col > log")

	if not os.path.isfile("Max_Tc_col_max_TC.col"):
		print("Tanimoto calculation failed")
		sys.exit()


	final_map_dict = {} ### this will only include molecules within your Tc range to ligand
	lig_tc_lower_bound = lig_tc_list[0]
	lig_tc_upper_bound = lig_tc_list[1]
	if os.path.isfile("Max_Tc_col_max_TC.col"):

		output = open("test_ligand_ids.smi", 'w')
		filtered_decoy_list = []
		open_max = open("Max_Tc_col_max_TC.col", 'r')
		read_max = open_max.readlines()
		open_max.close()
		
		for line in read_max:
			splt = line.strip().split(",")
			decoy_name = splt[0]
			tc = float(splt[1])
			if (tc >= lig_tc_lower_bound) and (tc < lig_tc_upper_bound):
				final_map_dict[decoy_name] = [tc, map_dict[decoy_name]]
				output.write(map_dict[decoy_name]+"\t"+decoy_name+"\n") ### write out the new SMILES for protonation
			else:
				print("FILTERING OUT "+decoy_name+" because of tanimoto of:",float(tc))

		output.close()
	

	return(final_map_dict)


#def remove_wrong_protomers(lig_props, ori_smiles_list, ori_comp_list, decoy_dict, pwd, prop_range_dict):
def remove_wrong_protomers(lig_props, map_dict, pwd, prop_range_dict):

	mwt = lig_props[0]
	logp = lig_props[1]
	rotB = lig_props[2]
	hbd = lig_props[3]
	hba = lig_props[4]
	q = lig_props[5]
	
	smi_file = pwd+"test_ligand_ids.smi" ### this file was written in run_tanimoto
	prot_dir = pwd+"test_ligand_protonation/"
	if not os.path.isdir(prot_dir):
		os.system("mkdir "+prot_dir)
	os.chdir(prot_dir)
	
	to_prot_file = prot_dir+"test_ligand_ids.smi"
	os.system("cp "+smi_file+" "+to_prot_file)
	
	#os.system("source /mnt/nfs/home/tbalius/.cshrc_dbgen_corina")
	os.system("source /mnt/nfs/home/rstein/.cshrc_dbgen_corina")
	#os.system("bash /mnt/nfs/export/rstein/DUDE_Z/new_DUDE_pipeline_4_2_2018/generate_protomers_from_smiles.sh -H 7.4 "+to_prot_file)
	os.system("bash /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/generate_protomers_from_smiles.sh -H 7.4 "+to_prot_file)
	
	built_dir = prot_dir+"/working/protonate/"
	if not os.path.isdir(built_dir):
		print("SMILES were not protonated. Quitting...Check that your ligand building pipeline is working.")
		sys.exit()

  #      built_file = built_dir+"test_ligand_ids-protomers.ism"
	built_file = built_dir+"test_ligand_ids-protomers-expanded.ism"
	
	new_decoy_dict = {}
	repeat_list = []
	filtered_decoy_list = []
	if os.path.isfile(built_file):
		open_built = open(built_file, 'r')
		read_built = open_built.readlines()
		open_built.close()
		
		for line in read_built:
			splitline = line.strip().split()
			dec_smiles = splitline[0]  ### we want to compare protomer SMILES properties to ligand properties, but we keep the original input SMILES
			dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q = get_stuff(dec_smiles, 'y')
			boolliee = compare(mwt, logp, rotB, hbd, hba, q, dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q, prop_range_dict)
			
			if boolliee == True:
				dec_zinc_ID = splitline[-2]
				lig_ID_count = repeat_list.count(dec_zinc_ID)
				dict_entry = dec_zinc_ID+"_"+str(lig_ID_count)
				repeat_list.append(dec_zinc_ID)
				
				if dec_zinc_ID not in new_decoy_dict: ## the first matching protomer is appended for this ZINC ID, although others may match
					tc_to_ligand = map_dict[dec_zinc_ID][0]
					mapped_dec_smiles = map_dict[dec_zinc_ID][1] ### make sure we append the input SMILES so we get same protomer SMILES during building
					list_entry = [mapped_dec_smiles, dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q, dec_smiles] ### added protomer SMILES 8/14/2019
					new_decoy_dict[dec_zinc_ID] = list_entry ### added protomer SMILES 8/14/2019
					filtered_decoy_list.append([tc_to_ligand, dec_zinc_ID, mapped_dec_smiles, dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q, dec_smiles])

	sorted_filtered_decoy_list = sorted(filtered_decoy_list)
	
	os.chdir(pwd)
	os.system("rm -rf "+prot_dir)
	
	return(sorted_filtered_decoy_list)


def write_out_approved_decoys(lig_list, approved_decoy_list):

	lig_smiles = lig_list[0]
	lig = lig_list[1]
	mw = lig_list[2]
	logp = lig_list[3]
	rotB = lig_list[4]
	HBD = lig_list[5]
	HBA = lig_list[6]
	q = lig_list[7]
	
	output = open(lig+"_property_matched_decoys.txt", 'w')
	output.write("%s %s %s %s %s %s %s %s %s %s\n" % ("SMILES", "ID", "MW", "LogP", "RotB", "HBD", "HBA", "Q", "PROT_ID", "Tc_to_Lig"))
	output.write("%s %s %s %f %f %d %d %d %d\n\n" % ("LIGAND: ", lig_smiles, lig, round(mw, 2), logp, rotB, HBD, HBA, int(q)))
	
	#[float(tc), second_comp]
	#[dec_smiles, dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q, dec_prot_id]
	
	app_dec_count = 1
	for app_dec in approved_decoy_list:
		print(app_dec)
		app_dec_tc = app_dec[0]
		app_dec_name = app_dec[1]
		app_dec_smiles = app_dec[2]
		app_dec_mw = app_dec[3]
		app_dec_logp = app_dec[4]
		app_dec_rotB = app_dec[5]
		app_dec_HBD = app_dec[6]
		app_dec_HBA = app_dec[7]
		app_dec_q = app_dec[8]
		app_dec_prot_id = app_dec[9]
		output.write("%s %s %s %f %f %d %d %d %d %s %f\n" % ("DECOY "+str(app_dec_count)+": ", app_dec_smiles, app_dec_name, app_dec_mw, app_dec_logp, app_dec_rotB, app_dec_HBD, app_dec_HBA, app_dec_q, app_dec_prot_id, app_dec_tc))
		app_dec_count += 1

	output.close()


def gather_ligands(ligand_file):

	open_lig = open(ligand_file,'r')
	read_lig = open_lig.readlines()
	open_lig.close()
	
	lig_dict = {line.strip().split()[1]:line.strip().split()[0] for line in read_lig if len(line.strip().split()) == 2}
	
	return(lig_dict)

def create_range_dict(infile):

	open_in = open(infile, 'r')
	read_in = open_in.readlines()
	open_in.close()

	### DEFAULTS FROM MYSINGER - these are used if property not specified
	#range_dict = {"CHG_RANGES": [  0,   0,   0,   0,   0,   1,   2],
	#"HBD_RANGES": [  0,   0,   1,   1,   2,   2,   3],
	#"HBA_RANGES": [  0,   1,   2,   2,   3,   3,   4],
	#"RB_RANGES": [  1,   2,   2,   3,   3,   4,   5],
	#"MWT_RANGES": [ 20,  35,  50,  65,  80, 100, 125],
	#"LOGP_RANGES": [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]}
	
# Example decoy_generation.in file
##############################
#   SMILES YES
#   PROTONATE YES
#   MWT 20 125
#   LOGP 0.4 3.6
#   RB 1 5
#   HBA 0 4
#   HBD 0 3
#   CHARGE 0 2
#   MINIMUM DECOYS PER LIGAND 20
#   DECOYS PER LIGAND 50
#   GENERATE DECOYS 750
#   MAXIMUM TC BETWEEN DECOYS 0.8
###################################

	range_dict = {"CHG_RANGES": [ 2],
	"HBD_RANGES":[ 3],
	"HBA_RANGES":[ 4],
	"RB_RANGES":[ 5],
	"MWT_RANGES":[ 125],
	"LOGP_RANGES":[ 3.6],
	"GENERATE_NUM":750,
	"TC":True,
	"min_num_decoys":20,
	"pref_num_decoys":50,
	"max_tc_between_decoys":0.8,
	"LIGAND_TC_RANGE": [0.0, 0.35]}
	
	for line in read_in:
		splitline = line.strip().split()
		if len(splitline) == 5:
			prop_type = splitline[0]
			if prop_type.lower() == "ligand":
				tc_lower_bound = float(splitline[3])
				tc_higher_bound = float(splitline[4])
				
				range_dict["LIGAND_TC_RANGE"] = [tc_lower_bound, tc_higher_bound]

		if len(splitline) == 3:
			prop_type = splitline[0]
			if prop_type.lower() == "generate":
				generate_num = int(splitline[-1])

			if prop_type.lower()[0] == "c":
				chg_high = int(splitline[-1])
				range_dict["CHG_RANGES"] = [chg_high]

			elif prop_type.lower() == "hbd":
				hbd_high = int(splitline[-1])
				range_dict["HBD_RANGES"] = [hbd_high]

			elif prop_type.lower() == "hba":
				hba_high = int(splitline[-1])
				range_dict["HBA_RANGES"] = [hba_high]

			elif prop_type.lower() == "rb":
				rb_high = int(splitline[-1])
				range_dict["RB_RANGES"] = [rb_high]

			elif prop_type.lower()[0] == "m":
				mwt_high = float(splitline[-1])
				range_dict["MWT_RANGES"] = [mwt_high]

			elif prop_type.lower()[0] == "l":
				logp_high = float(splitline[-1])
				range_dict["LOGP_RANGES"] = [logp_high]

	return(range_dict)


def main():

	pwd = os.getcwd()+"/"
	
	ligand_file = sys.argv[1]
	infile = pwd+"decoy_generation.in"
	
	if not os.path.isfile(infile):
		print("decoy_generation.in does not exist")
		sys.exit()
	else:
		prop_range_dict = create_range_dict(infile)

	print(prop_range_dict)

	lig_dict = gather_ligands(ligand_file)
	print(lig_dict)

	for lig in lig_dict:
		smiles = lig_dict[lig]
		query_zinc(lig, smiles, pwd, prop_range_dict, ligand_file)
	

main()
