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

                #print(mw, logp, rotB, HBA, HBD)
                #print("MW is ",mw)
                #print("logP is ",logp)
                #print("Rotatable Bonds is ",rotB)
                #print("HB Donors is ",HBD)
                #print("HB Acceptors is ",HBA)
                return(mw, logp, rotB, HBD, HBA, q)



def compare(mw, logp, rotB, HBD, HBA, q, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, step_count, prop_range_dict):

        #CHG_RANGES =  [  0,   0,   0,   0,   0,   1,   2]
        #HBD_RANGES =  [  0,   0,   1,   1,   2,   2,   3]
        #HBA_RANGES =  [  0,   1,   2,   2,   3,   3,   4]
        #RB_RANGES =   [  1,   2,   2,   3,   3,   4,   5]
        #MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
        #LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]

	CHG_RANGES = prop_range_dict["CHG_RANGES"]
	HBD_RANGES = prop_range_dict["HBD_RANGES"]
	HBA_RANGES = prop_range_dict["HBA_RANGES"]
	RB_RANGES = prop_range_dict["RB_RANGES"]
	MWT_RANGES = prop_range_dict["MWT_RANGES"]
	LOGP_RANGES = prop_range_dict["LOGP_RANGES"]

	mwt_round = int(round(mw))
	mwt_low = mwt_round-MWT_RANGES[step_count]
	mwt_high = mwt_round+MWT_RANGES[step_count]

	logp_round = round(logp, 2)
	logp_low = logp_round-LOGP_RANGES[step_count]
	logp_high = logp_round+LOGP_RANGES[step_count]

	rotB_low = rotB-RB_RANGES[step_count]
	rotB_high = rotB+RB_RANGES[step_count]
	
	HBD_low = HBD-HBD_RANGES[step_count]
	HBD_high = HBD+HBD_RANGES[step_count]
	
	HBA_low = HBA-HBA_RANGES[step_count]
	HBA_high = HBA+HBA_RANGES[step_count]
	
	CHG_low = int(q)-CHG_RANGES[step_count]
	CHG_high = int(q)+CHG_RANGES[step_count]
	
	if dec_mw >= mwt_low and dec_mw <= mwt_high:

		if dec_logp >= logp_low and dec_logp <= logp_high:

			if dec_rotB >= rotB_low and dec_rotB <= rotB_high:

				if dec_HBD >= HBD_low and dec_HBD <= HBD_high:

					if dec_HBA >= HBA_low and dec_HBA <= HBA_high:

						if dec_q >= CHG_low and dec_q <= CHG_high: 
							return(True)

	return(False) 

def run_tanimoto(smiles_list, compound_list, decoy_dict, ligand_file, lig_tc_list):

	### This new run_tanimoto function was added 8/14/2019.
	### Should speed up the calculation as only the ligand
	### is being compared to the decoys, not everything vs everything
	print("RUNNING TANIMOTO CALCULATION!")

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
		print("ERROR! The uint16.fp or uint16.count files do not exist for the ligands or decoys")
		sys.exit()
	else:
		os.system("/mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/cal_Tc_matrix_uint16/cal_Tc_matrix_uint16 "+dec_uint16_fp+" "+decoy_file+" "+dec_uint16_count+" "+lig_uint16_fp+" "+tot_ligand_file+" "+lig_uint16_count+" Max_Tc_col > log")

	if not os.path.isfile("Max_Tc_col_max_TC.col"):
		print("Tanimoto calculation failed!")
		sys.exit()

	lig_tc_lower_bound = lig_tc_list[0]
	lig_tc_upper_bound = lig_tc_list[1]
	filtered_decoy_list = []
	if os.path.isfile("Max_Tc_col_max_TC.col"):

		filtered_decoy_list = []
		open_max = open("Max_Tc_col_max_TC.col", 'r')
		read_max = open_max.readlines()
		open_max.close()

		for line in read_max:
			splt = line.strip().split(",")
			decoy_name = splt[0]
			tc = float(splt[1])
			if (tc >= lig_tc_lower_bound) and (tc < lig_tc_upper_bound): ## changed 3/6/2020
				list_entry = [float(tc), decoy_name]
				for i in range(len(decoy_dict[decoy_name])): ### add all decoy properties to list entry but sort by tanimoto to ligand
					list_entry.append(decoy_dict[decoy_name][i])
				filtered_decoy_list.append(list_entry)
			else:
				print("FILTERING OUT "+decoy_name+" because of tanimoto of:",float(tc))

	sorted_filtered_decoy_list = sorted(filtered_decoy_list)
	
	return(sorted_filtered_decoy_list)	



#def query_zinc(lig_id, smiles, mwt, logp, rotB, hbd, hba, q, prop_range_dict, ligand_file):
def query_zinc(lig_id, smiles, prop_range_dict, ligand_file):

	#MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
        #LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]

	MWT_RANGES = prop_range_dict["MWT_RANGES"]
	LOGP_RANGES = prop_range_dict["LOGP_RANGES"]
	CHG_RANGES = prop_range_dict["CHG_RANGES"]
	LIG_TC_RANGE = prop_range_dict["LIGAND_TC_RANGE"]

	mwt, logp, rotB, hbd, hba, q = get_stuff(smiles) ### maybe query from ZINC, don't manually calculate
	lig_list = [smiles, lig_id, mwt, logp, rotB, hbd, hba, q]
	
	mwt_round = int(round(mwt))
	logp_round = round(logp, 2)

	query_num = 1000
	decoy_count = 0
	decoy_dict = {}
        #while decoy_count < 3000:

        #smiles_list = [smiles]
        #compound_list = [lig_id]
	smiles_list = [] ### don't need to include ligand SMILES anymore because the Tanimoto calculation reads in separate files now (8/14/2019)
	compound_list = []

	filtered_decoy_list = []
	for i in range(7):
		if decoy_count > 1500: ## another arbitrary cutoff
			break

		mwt_low = mwt - MWT_RANGES[i]
		mwt_high = mwt + MWT_RANGES[i]
		
		logp_low = logp_round - LOGP_RANGES[i]
		logp_high = logp_round + LOGP_RANGES[i]
		
		chg_low = int(q) - CHG_RANGES[i]
		chg_high = int(q) + CHG_RANGES[i]
		
		#command3 = 'curl "http://zinc15/protomers.txt:smiles+zinc_id+mwt+logp+rb+hbd+hba+net_charge+tpsa+prot_id+model_type_name" -F "mwt-between=%d %d" -F "logp-between=%f %f" -F "sort=no" -F "net_charge=%d" -F "count=%d"' % (mwt_low, mwt_high, logp_low, logp_high, int(q), query_num)
		command3 = 'curl -k "https://zinc15.docking.org/protomers/subsets/usual.txt:smiles+zinc_id+mwt+logp+rb+hbd+hba+net_charge+prot_id+model_type_name" -F "mwt-between=%d %d" -F "logp-between=%f %f" -F "sort=no" -F "net_charge-between=%d %d" -F "count=%d"' % (mwt_low, mwt_high, logp_low, logp_high, int(chg_low), int(chg_high), query_num)

		print(command3)
		
		g = os.popen(command3)
		
		for line in g:
			splitline = line.strip().split()
			if len(splitline) > 6:
				dec_model_name = splitline[-1]
				if dec_model_name == "ref" or dec_model_name == "mid": ### want decoys built at reference pH
					dec_smiles = splitline[0]
					dec_zinc_id = splitline[1]
					dec_mwt = float(splitline[2])
					dec_logp = float(splitline[3])
					dec_rotB = int(splitline[4])
					dec_hbd = int(splitline[5])
					dec_hba = int(splitline[6])
					dec_q = float(splitline[7])
					dec_prot_id = int(splitline[8])
					
					boolliee = compare(mwt, logp, rotB, hbd, hba, q, dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q, i, prop_range_dict)
					if boolliee == True:
						if dec_zinc_id not in decoy_dict:
							decoy_count += 1
							decoy_dict[dec_zinc_id] = [dec_smiles, dec_mwt, dec_logp, dec_rotB, dec_hbd, dec_hba, dec_q, dec_prot_id]
							smiles_list.append(dec_smiles)
							compound_list.append(dec_zinc_id)
							print("DECOYS FOUND: ",decoy_count)


	print("FILTERING FOR TANIMOTO")
	if len(smiles_list) > 0:
		filtered_decoy_list = run_tanimoto(smiles_list, compound_list, decoy_dict, ligand_file, LIG_TC_RANGE)
	else:
		print("NO PROPERTY MATCHED DECOYS FOUND")
		sys.exit()

	if len(filtered_decoy_list) > 0:
		write_out_approved_decoys(lig_list, filtered_decoy_list)


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
		output.write("%s %s %s %f %f %d %d %d %d %d %f\n" % ("DECOY "+str(app_dec_count)+": ", app_dec_smiles, app_dec_name, app_dec_mw, app_dec_logp, app_dec_rotB, app_dec_HBD, app_dec_HBA, app_dec_q, app_dec_prot_id, app_dec_tc))
		app_dec_count += 1

	output.close()


def gather_ligands(ligand_file):

	open_lig = open(ligand_file,'r')
	read_lig = open_lig.readlines()
	open_lig.close()
	
	lig_dict = {line.strip().split()[1]: line.strip().split()[0] for line in read_lig if len(line.strip().split()) == 2}
	
	return(lig_dict)

def make_seven(low, high, prop):

	print(prop)
	seven_list = [low, high]
	
	if low == high:
		for i in range(5):
			seven_list.append(low)

	elif low != high:
		diff = (high - low)
		seven_diff = (float(diff) / float(6))
		new_term = low
		
		for i in range(5):
			new_term += seven_diff
			#print(new_term)
			
			if prop.lower()[0] == "l":
				seven_list.append(round(new_term, 1))
			elif prop.lower()[0] == "m":
				seven_list.append(round(new_term))
			else:
				seven_list.append(int(round(new_term)))

	return(sorted(seven_list))

def create_range_dict(infile):

##########################
## the format of the input file is shown below
#
# CHARGE 0 2
# HBD 0 3
# HBA 0 4
# RB 1 5
# MWT 20 125
# LOGP 0.4 3.6
# Tanimoto YES 
# MINIMUM DECOYS PER LIGAND 20 
# DECOYS PER LIGAND 50
# MAXIMUM TC BETWEEN DECOYS 0.8
#
#########################

	open_in = open(infile, 'r')
	read_in = open_in.readlines()
	open_in.close()
	
	### DEFAULTS FROM MYSINGER - these are used if property not specified
	### altered 10/10/2019
	range_dict = {"CHG_RANGES": [  0,   0,   0,   0,   0,   1,   2],
	"HBD_RANGES": [  0,   0,   1,   1,   2,   2,   3],
	"HBA_RANGES": [  0,   1,   2,   2,   3,   3,   4],
	"RB_RANGES": [  0,   1,   2,   3,   3,   4,   5],
	"MWT_RANGES": [ 0,  20,  40,  60,  80, 100, 125],
	"LOGP_RANGES": [0.0, 0.6, 1.2, 1.8, 2.4, 3.0, 3.6],
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
			if prop_type.lower()[0] == "c":
				chg_low = int(splitline[1])
				chg_high = int(splitline[2])
				
				chg_range = make_seven(chg_low, chg_high, "CHG")
				range_dict["CHG_RANGES"] = chg_range

			elif prop_type.lower() == "hbd":
				hbd_low = int(splitline[1])
				hbd_high = int(splitline[2])
				
				hbd_range = make_seven(hbd_low, hbd_high, "HBD")
				range_dict["HBD_RANGES"] = hbd_range

			elif prop_type.lower() == "hba":
				hba_low = int(splitline[1])
				hba_high = int(splitline[2])
				
				hba_range = make_seven(hba_low, hba_high, "HBA")
				range_dict["HBA_RANGES"] = hba_range

			elif prop_type.lower() == "rb":
				rb_low = int(splitline[1])
				rb_high = int(splitline[2])
				
				rb_range = make_seven(rb_low, rb_high, "RB")
				range_dict["RB_RANGES"] = rb_range

			elif prop_type.lower() == "mwt":
				mwt_low = float(splitline[1])
				mwt_high = float(splitline[2])
				
				mwt_range = make_seven(mwt_low, mwt_high, "MWT")
				range_dict["MWT_RANGES"] = mwt_range

			elif prop_type.lower() == "logp":
				logp_low = float(splitline[1])
				logp_high = float(splitline[2])
				
				logp_range = make_seven(logp_low, logp_high, "LOGP")
				range_dict["LOGP_RANGES"] = logp_range

				
	print(range_dict)
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

	lig_dict = gather_ligands(ligand_file)

	lig_count = 0
	for lig in lig_dict:
		smiles = lig_dict[lig]
		query_zinc(lig, smiles, prop_range_dict, ligand_file)
		lig_count += 1
	

main()
