import os, sys, time

###################################
# written by Reed Stein
# 10/2017 - 4/2018
# modified 5/3/2019
#####################################

def write_qsub(smiles_directory, curltype):

	output = open("dude_submit_ligands.csh",'w')
	output.write("#$ -S /bin/csh\n")
	output.write("#$ -cwd\n")
	output.write("#$ -q all.q\n")
	output.write("#$ -o dude_submit_stdout\n")
	output.write("#$ -e dude_submit_stderr\n\n")
	if curltype.lower() == "smiles":
		output.write("source /nfs/home/rstein/.cshrc_dbgen_corina\n")
	output.write("source /nfs/soft/jchem/current/env.csh\n")
	output.write("source /nfs/soft/python/envs/complete/current/env.csh\n")
	output.write("cd "+smiles_directory+"/ligand_${SGE_TASK_ID}/\n")
	output.write("echo `pwd`\n")
	if curltype.lower() == "smiles":
		output.write("python /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/generate_decoys_with_SMILES.py "+smiles_directory+"ligand_${SGE_TASK_ID}/ligand_${SGE_TASK_ID}.smi \n")
	elif curltype.lower() == "protomers":
		output.write("python /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/generate_decoys.py "+smiles_directory+"ligand_${SGE_TASK_ID}/ligand_${SGE_TASK_ID}.smi \n")
	output.close()


def read_infile(infile):

	open_in = open(infile, 'r')
	read_in = open_in.readlines()
	open_in.close()
	
	smiles_found = False
	for line in read_in:
		splitline = line.strip().split()
		if splitline[0].lower() == "smiles":
			smiles_found = True
			tf = splitline[-1].lower()
			if tf[0] == 'y':
				return(tf)
			if tf[0] == 'n':
				return(tf)

	if smiles_found == False:
		tf = 'no'
		return(tf)


def main():

	pwd = os.getcwd()+"/"
	
	smiles_dir = pwd+sys.argv[1]+"/"
	
	os.chdir(smiles_dir)
	
	infile = smiles_dir+"decoy_generation.in"
	if not os.path.isfile(infile):
		print("decoy_generation.in does not exist")
		sys.exit()
	
	prot_smi = read_infile(infile)
	if prot_smi.lower()[0] == "y":
		write_qsub(smiles_dir, "SMILES")
	else:
		write_qsub(smiles_dir, "PROTOMERS")

	lig_dirs = [name for name in os.listdir(".") if (os.path.isdir(name) and name.startswith("ligand_"))]
	lig_count = len(lig_dirs)
	print("TOTAL LIG DIRS:",lig_count)

	os.system("qsub -l h='n-9-23' -t 1-"+str(lig_count)+" dude_submit_ligands.csh")


main()
