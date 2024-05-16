import os, sys, gzip
#import numpy as np
#import time

############################################################################
# 3/5/2020 (AGAIN?!) Reed M Stein
#
# New version of filter_decoys.py used for protomers and SMILES
# Has the ability to NOT run the TC calculation if unnecessary
# Also sorts decoys by max Tc to any ligand, then assigns them 
# to ligands based on which ligand has the least # of decoys assigned
# Also allows modifying the maximum Tc allowed between decoys
#
##########################################################################

############################## UNUSED FUNCTIONS######################

def better_assign_decoys(lig_property_dict, decoy_property_dict, log_output):

        ### this function is not used in this iteration 5/3/2019

	current_dict = {}  ### now starting with unsorted dictionary
	for lig in lig_property_dict:
                current_dict[lig] = [] ### initialize dictionary

	print("LENGTH OF CURRENT_LIG_DICT = ",len(current_dict))
	print("LENGTH OF DECOY_PROPERTY_DICT = ",len(decoy_property_dict))
	
	for decoy in decoy_property_dict:
		prot_id = decoy_property_dict[decoy][8]
		if prot_id == "NA":
			mol2_file_present = True
		else:
			mol2file = find_path(prot_id)
			if mol2file == "0":
				print(decoy+" does not have a mol2file")
				continue
			else:
				mol2_file_present = True

		if mol2_file_present == True:
			mapped_lig = decoy_property_dict[decoy][10]
			closest_lig_tc_to_dec = decoy_property_dict[decoy][11]
			if mapped_lig in current_dict:
				current_dict[mapped_lig].append([closest_lig_tc_to_dec, decoy])

	print("FINISHED ASSIGNING ALL DECOYS TO LIGANDS")
	return(current_dict)

def randomly_choose_decoys(current_dict, lig_property_dict, decoy_property_dict, num_decoys, log_output):

        ### This function is not used in this iteration 5/3/2019

	final_dict = {}
	for lig_dir in current_dict:
		### sort so it's in order of highest Tc to known ligands
		sorted_decoy_id_list = [name[1] for name in sorted(current_dict[lig_dir])]
		#print(lig_dir, current_dict[lig_dir])
		if len(sorted_decoy_id_list) < num_decoys:
			print("There are not enough decoys for "+lig_dir, len(sorted_decoy_id_list))
			log_output.write("There are not enough decoys for "+lig_dir+"\n")
			continue

		elif len(sorted_decoy_id_list) == num_decoys:
			final_dict[lig_dir] = sorted_decoy_id_list

		else:

			final_dict[lig_dir] = []
			current_lig_dir_list = []
			
			#for i in range(num_decoys):
			#        final_dict[lig_dir].append(sorted_decoy_id_list[i])
			
			#output = open(lig_dir+"_replacements.txt", 'w')
			#for j in range(num_decoys, len(sorted_decoy_id_list)):
			#        rep_decoy = sorted_decoy_id_list[j]
			#        rep_decoy_properties = decoy_property_dict[rep_decoy]
			#        rep_dec_smiles = rep_decoy_properties[0]
			#        rep_dec_name = rep_decoy_properties[1]
			#        rep_prot_id = rep_decoy_properties[8]
			#        output.write(rep_dec_smiles+" "+rep_dec_name+" "+rep_prot_id+"\n")
			#output.close()
			
			#for i in range(len(current_dict[lig_dir])):
			#        rand_int = int(np.random.random(1)[0] * 1000000) ### assign a random number to each entry
			#       current_lig_dir_list.append([rand_int, i]) ### keep index with random number
			#       print(current_lig_dir_list)
			#       
			#sys.exit()
			print
			current_lig_dir_list = sorted(current_lig_dir_list) ### sort by this random number so we can randomly choose approved decoys for each ligand
			print(current_lig_dir_list)
			print
			
			for i in range(num_decoys):
				decoy_index_to_append = current_lig_dir_list[i][1] ### then add the decoy of that index to the final decoy list
				final_dict[lig_dir].append(current_dict[lig_dir][decoy_index_to_append]) ### add decoy properties to final dictionary
				print(final_dict[lig_dir])

			output = open(lig_dir+"_replacements.txt", 'w')
			for j in range(num_decoys,len(current_lig_dir_list)): ### save the rest as replacements
				j_decoy_index_to_append = current_lig_dir_list[j][1]
				rep_decoy = current_dict[lig_dir][j_decoy_index_to_append]
				rep_decoy_properties = decoy_property_dict[rep_decoy]
				rep_dec_smiles = rep_decoy_properties[0]
				rep_dec_name = rep_decoy_properties[1]
				rep_prot_id = rep_decoy_properties[8]
				output.write(rep_dec_smiles+" "+rep_dec_name+" "+rep_prot_id+"\n")
			output.close()

		lig_list = lig_property_dict[lig_dir]
		write_out_approved_decoys(lig_list, decoy_property_dict, final_dict)

	print("FINISHED ASSIGNING "+str(num_decoys)+" DECOYS TO ALL LIGANDS AND WRITING OUT REPLACEMENTS")

	return(final_dict)

#################### UNUSED FUNCTIONS ########################################################################################


def create_property_dict(infile):

##########################
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
#########################

	#### READ IN A BUNCH OF PARAMETERS FOR CALCULATION

	open_in = open(infile, 'r')
	read_in = open_in.readlines()
	open_in.close()
	
	### set the dictionary as defaults
	range_dict = {"CHG_RANGES": [ 2],
	"HBD_RANGES":[ 3],
	"HBA_RANGES":[ 4],
	"RB_RANGES":[ 5],
	"MWT_RANGES":[ 125],
	"LOGP_RANGES":[ 3.6],
	"TC":True,
	"min_num_decoys":20,
	"pref_num_decoys":50,
	"max_tc_between_decoys":0.8}
	
	for line in read_in:
		splitline = line.strip().split()

		if len(splitline) == 2:
			if splitline[0].lower() == "tanimoto":
				if splitline[1][0].lower() == "n":
					range_dict["TC"] = False

		if len(splitline) == 3:
			prop_type = splitline[0]
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


		if len(splitline) == 4:
			if splitline[0].lower() == "decoys":
				pref_num_decoys = int(splitline[-1])
				range_dict["pref_num_decoys"] = pref_num_decoys

		if len(splitline) == 5:
			if splitline[0].lower() == "minimum":
				min_num_decoys = int(splitline[-1])
				range_dict["min_num_decoys"] = min_num_decoys

			if splitline[0].lower() == "maximum":
				max_tc_bn_decs = float(splitline[-1])
				range_dict["max_tc_between_decoys"] = max_tc_bn_decs

	return(range_dict)

def collect_decoys(lig_property_dict, decoy_property_dict, decoy_file, dec_list, keep_lig_or_not):

	decoy_file_name = decoy_file.split("_property_")[0]
	open_decoy = open(decoy_file, 'r')
	read_decoy = open_decoy.readlines()
	open_decoy.close()
	
	for line in read_decoy:
		splitline = line.strip().split()
		if len(splitline) > 1:
			if splitline[0] == "LIGAND:":
				lig_smiles = splitline[1]
				lig_name = splitline[2]
				lig_mw = splitline[3]
				lig_logp = splitline[4]
				lig_rotB = splitline[5]
				lig_HBD = splitline[6]
				lig_HBA = splitline[7]
				lig_Q = splitline[8]
				
				if keep_lig_or_not[0] == "y":
					lig_property_dict[lig_name] = [lig_smiles, lig_name, lig_mw, lig_logp, lig_rotB, lig_HBD, lig_HBA, lig_Q]

                                #LIGAND:  CC(=O)NC1=C(C=CS1)C([O-])=O 322_0 184.010000 0.070000 2 1 4 -1
                                # "SMILES", "ID", "MW", "LogP", "RotB", "HBD", "HBA", "Q"

			elif splitline[0] == "DECOY":
				dec_smiles = splitline[2]
				dec_name = splitline[3]
				dec_mw = float(splitline[4])
				dec_logp = splitline[5]
				dec_rotB = splitline[6]
				dec_HBD = splitline[7]
				dec_HBA = splitline[8]
				dec_q = splitline[9]
				dec_prot_id = splitline[10] ### may be an actual number or the protomer SMILES depending on which script is used
				dec_tc_to_lig = float(splitline[11])
				
				dict_entry = [dec_smiles, dec_name, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, dec_prot_id, dec_tc_to_lig, lig_name]
				
				if dec_name not in decoy_property_dict:
					decoy_property_dict[dec_name] = dict_entry
					dec_list.append([dec_tc_to_lig, dec_smiles, dec_name]) # sorting by Tc to any ligand

	return(lig_property_dict, decoy_property_dict, dec_list)

def collect_decoys_calc_tanimoto(system_dir, lig_dir_list, prop_range_dict, log_output):

        #txt_files = [[system_dir+"/"+name+"/"+f for f in os.listdir(name) if (os.path.isfile(os.path.join(name,f)) and f.endswith("_property_matched_decoys.txt")] for name in lig_dir_list]
        #txt_files = [[system_dir+"/"+name+"/"+f for f in os.listdir(name) if os.path.isfile(os.path.join(name,f)) and f.endswith("_property_matched_decoys.txt")] for name in lig_dir_list]

	min_num_decoys = prop_range_dict["min_num_decoys"]
	calc_tanimoto = prop_range_dict["TC"]
	
	txt_files = []
	for l in lig_dir_list:
		fil_list = [name for name in os.listdir(system_dir+"/"+l+"/") if os.path.isfile(system_dir+"/"+l+"/"+name)]
		for f in fil_list:
			if f.endswith("_property_matched_decoys.txt"):
				txt_files.append(system_dir+"/"+l+"/"+f)

	lig_dir_list = []
	for t in txt_files:
		sum_lines = sum(1 for line in open(t))
		lig_dir_list.append([sum_lines, t])
 #       sorted_lig_dir_list = sorted([[sum(1 for line in open(t[0])),t[0]] for t in txt_files])
	sorted_lig_dir_list = sorted(lig_dir_list)

	lig_property_dict = {} ## contains ligands as keys and ligand properties
	decoy_property_dict = {} ## contains decoys as keys and decoy properties
	assigned_decoy_dict = {}
	dec_list = []
	repeat_list = []
	lig_count = 0
	for lig in sorted_lig_dir_list:
		sum_lines = lig[0]
		decoy_file_path = lig[1]
		lig_count += 1
		lig_property_dict, decoy_property_dict, dec_list = collect_decoys(lig_property_dict, decoy_property_dict, decoy_file_path, dec_list, "yes")


	print("TOTAL DECOYS:",len(dec_list))
	if len(dec_list) == 0:
		print("NO DECOYS FOUND")
		sys.exit()

	if calc_tanimoto == False:
		if os.path.isfile("test_decoy_smiles.smi") and os.path.isfile("cluster_head.list"):
			print("skipping TC calc")
			log_output.write("Skipping Tanimoto calculation\n")
			open_test = open("test_decoy_smiles.smi", 'r')
			read_test = open_test.readlines()
			open_test.close()
			
			smiles_list = [line.strip().split()[0] for line in read_test]
			zinc_id_list = [line.strip().split()[1] for line in read_test]
		else:
			print("This step requires 'test_decoy_smiles.smi' and 'cluster_head.list'. Need to run TC calculation...\n")
			log_output.write("This step requires 'test_decoy_smiles.smi' and 'cluster_head.list'. Need to run TC calculation...\n")
			sys.exit()

	else:


		dec_list = sorted(dec_list) ### sorted by Tc to any ligand so can maintain most dissimilar decoys at later step
		
		smiles_list = [] ### only includes decoys
		zinc_id_list = [] ### only includes decoys
		final_output = open("test_decoy_smiles.smi", 'w')
		for tc, smiles, zinc_id in dec_list:
			final_output.write(smiles+" "+zinc_id+"\n")
			smiles_list.append(smiles)
			zinc_id_list.append(zinc_id)

		final_output.close()

		tc_cutoff = prop_range_dict["max_tc_between_decoys"]
		#os.system("python ~jklyu/zzz.github/ChemInfTools/utils/teb_chemaxon_cheminf_tools/generate_chemaxon_fingerprints.py test_decoy_smiles.smi test_decoy_smiles")
		os.system("python /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/generate_chemaxon_fingerprints.py test_decoy_smiles.smi test_decoy_smiles")
		os.system("/mnt/nfs/home/jklyu/zzz.github/ChemInfTools/utils/best_first_clustering/best_first_clustering test_decoy_smiles.fp test_decoy_smiles.smi "+str(tc_cutoff)+" "+str(len(smiles_list)))
		#sys.exit()

	return(lig_property_dict, decoy_property_dict, lig_count, smiles_list, zinc_id_list)

def discard_similar_decoys(decoy_property_dict, clust_head_file):
	
	open_c = open(clust_head_file, 'r')
	read_c = open_c.readlines()
	open_c.close()
	
	clustheads = {line.strip().split(",")[2]:[] for line in read_c}
	new_decoy_prop_dict = {d:decoy_property_dict[d] for d in decoy_property_dict if d in clustheads}
	
	print("LENGTH OF DISCARD IS",len(decoy_property_dict) - len(new_decoy_prop_dict))
	print("LENGTH OF NEW DECOY PROP DICT = ",len(new_decoy_prop_dict))
	
	return(new_decoy_prop_dict)


def compare_properties(mw, logp, rotB, HBD, HBA, q, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, prop_range_dict):

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
	
	CHG_range = CHG_RANGES[-1]
	HBD_range = HBD_RANGES[-1]
	HBA_range = HBA_RANGES[-1]
	RB_range = RB_RANGES[-1]
	MWT_range = MWT_RANGES[-1]
	LOGP_range = LOGP_RANGES[-1]
	
	mwt_round = int(round(mw))
	mwt_low = mwt_round-MWT_range
	mwt_high = mwt_round+MWT_range
	
	logp_round = round(logp, 2)
	logp_low = logp_round-LOGP_range
	logp_high = logp_round+LOGP_range
	
	rotB_low = rotB-RB_range
	rotB_high = rotB+RB_range
	
	HBD_low = HBD-HBD_range
	HBD_high = HBD+HBD_range
	
	HBA_low = HBA-HBA_range
	HBA_high = HBA+HBA_range
	
	CHG_low = int(q)-CHG_range
	CHG_high = int(q)+CHG_range
	
	if dec_mw >= mwt_low and dec_mw <= mwt_high:

		if dec_logp >= logp_low and dec_logp <= logp_high:

			if dec_rotB >= rotB_low and dec_rotB <= rotB_high:

				if dec_HBD >= HBD_low and dec_HBD <= HBD_high:
                        
					if dec_HBA >= HBA_low and dec_HBA <= HBA_high:

						if dec_q >= CHG_low and dec_q <= CHG_high:

							return(True)

	return(False)

def check_int(prot_id):

	try:
		int(prot_id)
		return(True)

	except ValueError:
		return(False)


def assign_decoys_deux(lig_property_dict, decoy_property_dict, prop_range_dict, log_output):

	### this function assigns matching decoys to ligands based on which has the least # of decoys already assigned

	min_num_decoys = prop_range_dict["min_num_decoys"]
	pref_num_decoys = prop_range_dict["pref_num_decoys"]
	
	current_dict = {lig:[] for lig in lig_property_dict}
	
	print("LENGTH OF CURRENT_LIG_DICT = ",len(current_dict))
	print("LENGTH OF DECOY_PROPERTY_DICT = ",len(decoy_property_dict))
	
	sorted_decoy_tc_list = sorted([[float(decoy_property_dict[decoy][9]), decoy] for decoy in decoy_property_dict]) ### sort decoys by lowest TC to any ligand
	#### [dec_smiles, dec_name, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, dec_prot_id, dec_tc_to_lig, lig_name]
	#### ['COC(=O)c1coc(CN[C@@H](C(=O)OC)c2ccc(OC)c(O)c2)c1', 'ZINC000274498699', 349.116152, '1.784300', '7', '2', '8', '0', 'COC(=O)[C@H](NCC1=CC(=CO1)C(=O)OC)C1=CC=C(OC)C(O)=C1', 0.242857, '118_0']
	
	for tc, decoy in sorted_decoy_tc_list:
		#decoy = pair[1]
		decoy_map = decoy_property_dict[decoy]
		#prot_id = decoy_property_dict[decoy][8]
		prot_id = decoy_map[8]
		int_or_not = check_int(prot_id)
		if int_or_not == False:
			mol2_file_present = True
		else:
			mol2file = find_path(prot_id)
			if mol2file == "0":
				print(decoy+" does not have a mol2file")
				continue
			else:
				mol2_file_present = True

		if mol2_file_present == True:
	#		mapped_lig = decoy_property_dict[decoy][10]
			#closest_lig_tc_to_dec = decoy_property_dict[decoy][11]
			closest_lig_tc_to_dec = decoy_map[9]
			#if mapped_lig in current_dict:
			#	current_dict[mapped_lig].append([closest_lig_tc_to_dec, decoy])
			
			dec_mw = float(decoy_map[2])
			dec_logp = float(decoy_map[3])
			dec_rotB = int(decoy_map[4])
			dec_HBD = int(decoy_map[5])
			dec_HBA = int(decoy_map[6])
			dec_Q = int(decoy_map[7])
			
			true_list = []
			for lig in lig_property_dict:
				lig_mw = float(lig_property_dict[lig][2])
				lig_logp = float(lig_property_dict[lig][3])
				lig_rotB = int(lig_property_dict[lig][4])
				lig_HBD = int(lig_property_dict[lig][5])
				lig_HBA = int(lig_property_dict[lig][6])
				lig_Q = int(lig_property_dict[lig][7])
				
				#time.sleep(15)
				boolliee = compare_properties(lig_mw, lig_logp, lig_rotB, lig_HBD, lig_HBA, lig_Q, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_Q, prop_range_dict)
				if boolliee == True:
					true_list.append(lig)

			if len(true_list) == 0:
				print(decoy+" does not have a match")
				continue

			sort_len_list = sorted([[len(current_dict[lig]),lig] for lig in true_list])
			#for lig in true_list:
			#        current_dict_len = len(current_dict[lig])
			#        len_list.append([current_dict_len, lig])
			
			#sort_len_list = sorted(len_list)
			
			lig_assigned = sort_len_list[0][1]
			current_dict[lig_assigned].append(decoy)

	num_decoys = calc_max_decoys_available(min_num_decoys, pref_num_decoys, current_dict) 
	#print("RETRIEVING "+str(num_decoys)+" for "+str(len(current_dict))+" ligands")
	#log_output.write("RETRIEVING "+str(num_decoys)+" for "+str(len(current_dict))+" ligands\n")
	
	final_dict = {}
	for lig in current_dict:
		#print(lig, len(current_dict[lig]))
		log_output.write(lig+" ACCEPTABLE DECOYS:"+str(len(current_dict[lig]))+"\n")
		if len(current_dict[lig]) < num_decoys:
			print("THERE ARE NOT ENOUGH DECOYS FOR "+lig+"\n")
			log_output.write("THERE ARE NOT ENOUGH DECOYS FOR "+lig+"\n")
			continue

		elif len(current_dict[lig]) == num_decoys:
			final_dict[lig] = current_dict[lig]
		
		else:
			final_dict[lig] = []
			current_lig_dir_list = []
			
			for i in range(num_decoys):
				final_dict[lig].append(current_dict[lig][i])
			
			output = open(lig+"_replacements.txt", 'w')
			for j in range(num_decoys, len(current_dict[lig])):
				rep_decoy = current_dict[lig][j]
				rep_decoy_properties = decoy_property_dict[rep_decoy]
				rep_dec_smiles = rep_decoy_properties[0]
				rep_dec_name = rep_decoy_properties[1]
				rep_prot_id = rep_decoy_properties[8]
				output.write(rep_dec_smiles+" "+rep_dec_name+" "+rep_prot_id+"\n")
			output.close()

		lig_list = lig_property_dict[lig]
		write_out_approved_decoys(lig_list, decoy_property_dict, final_dict)
			

	print("FINISHED ASSIGNING ALL DECOYS TO LIGANDS")
	return(current_dict)

def count_gz_lines(gz_file):

	protomer_open = gzip.open(gz_file, 'rb')
	count = 0
	for line in protomer_open:
		splt = line.strip().split()
		if len(splt) > 0:
			if splt[0] == "M":
				count += 1
				break
		count += 1
	return(count)



def find_path(prot_id):

	#print(prot_id)
	lp = len(prot_id)
	if lp < 7:
		print("Error")
		sys.exit()
	dir1 = prot_id[lp-2:lp]
	dir2 = prot_id[lp-4:lp-2]
	dir3 = prot_id[lp-6:lp-4]
	#print(dir1, dir2, dir3)
	
	# /nfs/dbraw/zinc/49/76/33/419497633.db2.gz
	db2file = "/nfs/dbraw/zinc/"+dir3+"/"+dir2+"/"+dir1+"/"+prot_id+".db2.gz"
	#mol2file = "/nfs/dbraw/zinc/"+dir3+"/"+dir2+"/"+dir1+"/"+prot_id+".mol2.gz"
	if not os.path.isfile(db2file):
		print(db2file, prot_id+" has not been built")
		return("0")
	if os.path.isfile(db2file):
		try:
			num_mol2_lines = count_gz_lines(db2file)
		except IOError:
			return("0")
		if num_mol2_lines > 0:
			return(db2file)

	return("0")

def calc_max_decoys_available(min_num_decoys, pref_num_decoys, current_dict):

	min_num_decoy_list = [[len(current_dict[lig]), lig] for lig in current_dict]
	
	#min_num_decoys_to_write = min_num_decoys ### assign the user-specified minimum number of decoys for all ligands
	max_num_decoys_list = []
	for i in range(pref_num_decoys, min_num_decoys-1, -1): ## counting down from the preferred number of decoys to the minimum number of decoys
		len_min_num_list_above = len([lig[1] for lig in min_num_decoy_list if lig[0] >= i]) ### count how many ligands have more or equivalent decoys assigned than i
		max_num_decoys_list.append([float(len_min_num_list_above) / float(len(current_dict)), i]) ### add the percentage of total ligands that have more/equivalent decoys than i,
	
	max_num_decoys = sorted(max_num_decoys_list,reverse=True) ## sort the number of decoys assigned by this percentage
	decs_assigned = False
	for perc, num_decs_assigned in max_num_decoys:
		if perc >= 0.6: ### 60% (a totally arbitrary value) of the ligands should have at least ${num_decs_assigned} decoys assigned
			decs_assigned = True
			print("DECOY NUMBER = ", num_decs_assigned)
			min_num_decoys_to_write = num_decs_assigned
			break ### break because we found the maximum 

	if decs_assigned == False:
		min_num_decoys_to_write = max_num_decoys[0][1] ### if <60% of the ligands have ${min_num_decoys} assigned, then just set the number of decoys to write to be
							       ### the number of decoys that the most ligands have assigned

	return(min_num_decoys_to_write)


def write_out_approved_decoys(lig_list, decoy_property_dict, final_decoy_dict):

	print("WRITING OUT DECOYS FOR "+lig_list[1])
	
	lig = lig_list[1]
	
	#[dec_smiles, dec_name, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, dec_prot_id]
	output = open(lig+"_final_property_matched_decoys.txt", 'w')
	output.write("%s %s %s %s %s %s %s %s %s %s\n" % ("SMILES", "ID", "MW", "LogP", "RotB", "HBD", "HBA", "Q", "PROT_ID", "Max_TC_to_LIGS"))
	output.write("%s %s %s %f %f %d %d %d %f\n\n" % ("LIGAND: ", lig_list[0], lig, round(float(lig_list[2]), 2), float(lig_list[3]), int(lig_list[4]), int(lig_list[5]), int(lig_list[6]), float(lig_list[7])))
	
	#[lig_smiles, lig_name, lig_mw, lig_logp, lig_rotB, lig_HBD, lig_HBA, lig_Q]
	
	app_dec_count = 1
	for i in range(len(final_decoy_dict[lig])):
		decoy = final_decoy_dict[lig][i]
		app_dec = decoy_property_dict[decoy]
		
		#dec_smiles, dec_name, dec_mw, dec_logp, dec_rotB, dec_HBD, dec_HBA, dec_q, dec_tpsa, dec_prot_id
		app_dec_smiles = app_dec[0]
		app_dec_name = app_dec[1]
		app_dec_mw = float(app_dec[2])
		app_dec_logp = float(app_dec[3])
		app_dec_rotB = int(app_dec[4])
		app_dec_HBD = int(app_dec[5])
		app_dec_HBA = int(app_dec[6])
		app_dec_q = int(app_dec[7])
		app_dec_prot_id = str(app_dec[8])
		app_dec_tc_to_lig = str(app_dec[9])
		#print("TC = ",str(app_dec[9]))
		output.write("%s %s %s %f %f %d %d %d %d %s %s\n" % ("DECOY "+str(app_dec_count)+": ", app_dec_smiles, app_dec_name, app_dec_mw, app_dec_logp, app_dec_rotB, app_dec_HBD, app_dec_HBA, app_dec_q, app_dec_prot_id, app_dec_tc_to_lig))
		app_dec_count += 1

	output.close()


def main():

	pwd = os.getcwd()+"/"
	
	smiles_dir = pwd
	infile = pwd+"decoy_generation.in"
	
	if not os.path.isfile(infile):
		print("decoy_generation.in does not exist")
		sys.exit()
	else:
		prop_range_dict = create_property_dict(infile) 

	lig_dir_list = [name for name in os.listdir(smiles_dir) if os.path.isdir(name)]
	
	log_output = open(smiles_dir+"FILTER_DECOYS.log", 'w')
	
	lig_property_dict, decoy_property_dict, lig_num, smiles_list, zinc_id_list = collect_decoys_calc_tanimoto(smiles_dir, lig_dir_list, prop_range_dict, log_output)
	#print(decoy_property_dict)
	
	print("LIG NUM = ",lig_num)
	
	clusterhead_file = pwd+"cluster_head.list"
	trim_decoy_prop_dict = discard_similar_decoys(decoy_property_dict, clusterhead_file)
	
	print("NUMBER OF DECOYS:",len(trim_decoy_prop_dict),"from",len(lig_property_dict),"PROTOMERS")
	log_output.write("NUMBER OF DECOYS: "+str(len(trim_decoy_prop_dict))+" from "+str(len(lig_property_dict))+" PROTOMERS\n")
	
	target_number_decoys = len(trim_decoy_prop_dict) / len(lig_property_dict) ### for those with X or more decoys, if there are not enough decoys in total to assign to the ligands, quit
	min_num_decoys = prop_range_dict["min_num_decoys"]
	if target_number_decoys < min_num_decoys:
		print("There are not enough decoys for all ligands.")
		print("There are only "+str(len(trim_decoy_prop_dict))+" for "+str(len(lig_property_dict))+" and your minimum specification is "+str(min_num_decoys))
		print("Either run 'generate_decoys.py' again, or lower your minimum count for decoys assigned.")
		print(target_number_decoys, len(trim_decoy_prop_dict), len(lig_property_dict))
		log_output.write("There are not enough decoys for all ligands\n")
		log_output.write("There are only "+str(len(trim_decoy_prop_dict))+" for "+str(len(lig_property_dict))+" and your minimum specification is "+str(min_num_decoys)+"\n")
		log_output.write("Either run 'generate_decoys.py' again, or lower your minimum count for decoys assigned\n")
		log_output.close()
		sys.exit()
	
	### smiles_list and zinc_id_list contain ligands and decoys - I am only comparing decoys to filter out ones that are too similar
        ### smiles_list and zinc_id_list are sorted by molecular weight in the previous script for best first clustering

	ori_full_dict = assign_decoys_deux(lig_property_dict, trim_decoy_prop_dict, prop_range_dict, log_output)

main()
