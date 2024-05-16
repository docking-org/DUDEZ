import os, sys

###################################
# written by Reed Stein
# 10/2017 - 4/2018
# modified 5/3/2019
##################################

def write_qsub(smiles_directory):

	output = open("filter_decoys.csh",'w')
	output.write("#$ -S /bin/csh\n")
	output.write("#$ -cwd\n")
	output.write("#$ -q all.q\n")
	#output.write("#$ -q matt.q\n")
	output.write("#$ -o filter_decoys_submit_stdout\n")
	output.write("#$ -e filter_decoys_submit_stderr\n\n")
	#output.write("source ~tbalius/.cshrc_dbgen_corina\n")
	output.write("source /nfs/soft/jchem/current/env.csh\n")
	output.write("source /nfs/soft/python/envs/complete/current/env.csh\n")
	output.write("cd "+smiles_directory+"/\n")
	output.write("echo `pwd`\n")
	output.write("python /mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/filter_decoys.py\n")
	output.close()



def main():

	pwd = os.getcwd()+"/"
	
	smiles_dir = pwd+sys.argv[1]+"/"
	
	os.chdir(smiles_dir)
	
	write_qsub(smiles_dir)
	
	os.system("qsub filter_decoys.csh")


main()
