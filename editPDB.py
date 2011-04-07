#############################################
#                                           #
#	PDB files manipulation                  #
#                    @author:TRELLET Mikael #
#                    @date: January 2011    #
#############################################


from os import popen, _exit, environ, pathsep
from sys import argv, stderr, stdout
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import PPBuilder

# GLOBAL DEFINITIONS

pdb_name=''
HELP = 0

def printUsage(name, exitCode) :
    if (not HELP) :
        stderr.write("\n\t --- BAD ARGUMENTS (see Usage) ! ---\n\n")
        stderr.write(" To print help on this program, use :   " + name + " -h   (or --help)\n")        
    else :
        stdout.write("\n" + name + " : Perform several modifications in a pdb file, adding specific suffixe to the name of the input file\n")
        stdout.write("To use it : python getFreeFromPDB.py {pdb_file_name} \n\n")
        
    _exit(exitCode)
    return

def testCmd(argv) :
    global HELP
    try:                            # ---------- do you need some help ?
        argv.index("-h")
        HELP = 1
        printUsage(argv[0], 0)
    except:
        try:
            argv.index("--help")
            HELP= 1
            printUsage(argv[0], 0)
        except:
            HELP = 0
	return

def removeDoubleAtoms():# Remove all double atoms defined in a pdb and save the new structure in pdbname_noDouble.pdb
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()
	
	structure.remove_disordered_atoms()

	w = PDBIO()
	w.set_structure(structure)
	w.save(nameStruct+'_noDouble.pdb')

def deleteChain():# Delete a complete chain from a pdb and save the new structure in pdbname_free.pdb
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()
	seq=''
	
	nb_chain=input('How many chain do you want to delete : ')
	for i in range(nb_chain):
		rm_chain=raw_input('What chain you want to delete : ')
		for model in structure:
			for chain in model:
				if(chain.id==rm_chain):
					model.detach_child(chain.id)
	pept = raw_input('Do you want to get a pdb with the sequence in its name : ')
	if(pept == 'y'):
		ppb=PPBuilder()
		for pp in ppb.build_peptides(structure):
			seq = seq + pp.get_sequence()
		seq=seq.lower()
		seq=str(seq)
		w = PDBIO()
		w.set_structure(structure)
		w.save(seq+'_bound.pdb')
	else:
		w = PDBIO()
		w.set_structure(structure)
		w.save(nameStruct+'_without'+rm_chain+'.pdb')

def deleteResidue():# Delete a residue from a pdb and save the new structure in pdbname_noResidue.pdb
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()

	rm_residue=raw_input('What residue you want to delete : ')
	for model in structure:
		for chain in model:
			for residue in chain:
				print residue.id
				if(residue.id[1]==rm_residue):
					print 'HELLO'
					chain.detach_child(residue.id)

	w = PDBIO()
	w.set_structure(structure)
	w.save(nameStruct+'_noResidue.pdb')
	
def removeHetero():# Remove all heteroatoms from a pdb and save the new structure in pdbname_noHetero.pdb
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()
	for model in structure:
		for chain in model:
			for residue in chain:
				id = residue.id				
				if id[0] != ' ':
					chain.detach_child(residue.id)
			if len(chain) == 0:
				model.detach_child(chain.id)

	w = PDBIO()
	w.set_structure(structure)
	w.save(nameStruct+'_noHetero.pdb')

def getSequence(): # Get the sequence of a specific chain
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()
	seq=''
	
	what_chain=raw_input('For what chain do you want the sequence : ')

	for model in structure:
		for chain in model:
			if chain.id != what_chain:
				model.detach_child(chain.id)

	ppb=PPBuilder()
	for pp in ppb.build_peptides(structure):
		seq = seq + pp.get_sequence()
	seq=seq.lower()
	print seq

def renumberChain(): # Allow to renumber from what you want a specific chain
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()
	
	what_chain=raw_input('What is the chain you want to renumber : ')
	number=input('What is the first number of the chain : ')

	for model in structure:
		for chain in model:
			if chain.id == what_chain:
				for residue in chain:
					if residue.id[0] == ' ':
						residue.id=(' ', number, ' ')
						number=number+1
					else:
						chain.detach_child(residue.id)

	w = PDBIO()
	w.set_structure(structure)
	w.save(nameStruct+'_ren.pdb')

def assembleChain(): # Allow to assemble 2 chains together
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()

	what_chain=raw_input('What is the 1st chain you want to assemble : ')
	what_chain2=raw_input('What is the 2nd chain you want to assemble : ')

	for model in structure:
		for chain in model:
			if chain.id == what_chain:
				parent=chain;
			elif chain.id == what_chain2:
				for residue in chain:
					residue.get_parent().id=what_chain

	w = PDBIO()
	w.set_structure(structure)
	w.save(nameStruct+'_assemble.pdb')
	
def renameChain():
	parser = PDBParser()
	nameStruct=pdb_name.partition('.')[0]
	structure = parser.get_structure(nameStruct, pdb_name)
	header = parser.get_header()
	trailer = parser.get_trailer()
	
	what_chain=raw_input('What is the chain you want to rename : ')
	what_chain2=raw_input('What is the new name of this chain : ')
	
	for model in structure:
		for chain in model:
			if chain.id == what_chain:
				chain.id = what_chain2
				
	w = PDBIO()
	w.set_structure(structure)
	w.save(nameStruct+'_rename.pdb')

def main(argv):
	global pdb_name
	
	testCmd(argv)
	
	try:
        	pdb_name = argv[1] # file.pdb
    	except:
        	stderr.write("\n [ERROR] Impossible to open the input file <file.pdb>\n\n")
        	printUsage(argv[0], 0)

	choix = raw_input('Delete a chain (d), delete a residue (dr), remove heteroatoms (rh), remove double atomes (rda), assemble 2 chains (ac), renumber a chain (rc), rename a chain (rnc) or get the sequence (gs) : ')
	if choix == 'd':
		deleteChain()
	elif choix == 'rh':
		removeHetero()
	elif choix == 'gs':
		getSequence()
	elif choix == 'dr':
		deleteResidue()
	elif choix == 'rc':
		renumberChain()
	elif choix == 'rnc':
		renameChain()
	elif choix == 'ac':
		assembleChain()
	elif choix == 'rda':
		removeDoubleAtoms()

end = main(argv)

#nomF = raw_input('Nom du fichier a traiter : ')


	
