import os
from Bio.PDB.PDBParser import PDBParser
from distUtils import get_distance
from Bio.PDB.Polypeptide import three_to_one

# gap_len = 1
# found_dist = (13.422670354004389, 12.01411067192551, 15.601742709167148, 9.833812337568627)

def getScore(dist1, dist2):
	return abs(dist1[0]-dist2[0]) + abs(dist1[1]-dist2[1]) + abs(dist1[2]-dist2[2] + abs(dist1[3]-dist2[3]))

def get_first_last_atoms(res, i):
	start_index = 0
	end_index = 0
	first_found = False
	index = 0
	for atom in res:
		if not first_found and atom.get_id() == "N":
			start_index=index
			first_found = True
		if atom.get_id() == "CB":
			end_index=index
		index+=1
	return (start_index + i, end_index + i)

def get_best_res(target_dist, label=None):
	minimal_score=1000000000000000000
	best_part=None
	ava=0
	for subdir, dirs, files in os.walk("top500H"):
	    for file in files:
	    	try:
		        ava+=1
		        filepath = subdir + os.sep + file
		        p = PDBParser(PERMISSIVE=1, QUIET=True)
		        structure = p.get_structure('file', filepath)
		        residues = [residue for model in structure for chain in model for residue in chain]
		        atoms = [atom for residue in residues for atom in residue]
		        print("Processing ", file, ava/5, "%")
		        for i in range(0,len(residues)):
		        	res = residues[i]
		        	if label == None or label == three_to_one(res.get_resname()):
		        		first, last = get_first_last_atoms(res, i)
			        	dist = get_distance(atoms, first, last)
			        	if getScore(dist, target_dist) < minimal_score:
			        		minimal_score = getScore(dist, target_dist)
			        		best_part = res
	    	except Exception:
    			pass
	return best_part


def get_gaps_parts(gaps_dict, gap_label=None):
	res = {}
	for k in gaps_dict.keys():
		v = gaps_dict[k]
		if v == "#NA":
			continue
		label = None
		if gap_label != None:
			label = gap_label[k]
		best_res = get_best_res(v, label)
		res[k] = best_res
	return res

def get_insertion_type(al0, al1):
	res = {}
	for i in range(0,len(al1)):
		if al1[i] == '-':
			res[i] = al0[i]
	return res

# for subdir, dirs, files in os.walk("top500H"):
#     for file in files:
#     	try:
# 	        ava+=1
# 	        filepath = subdir + os.sep + file
# 	        p = PDBParser(PERMISSIVE=1, QUIET=True)
# 	        structure = p.get_structure('file', filepath)
# 	        residues = [residue for model in structure for chain in model for residue in chain]
# 	        atoms = [atom for residue in residues for atom in residue]
# 	        print("Processing ", file, ava/5, "%")
# 	        for i in range(0,len(residues)):
# 	        	res = residues[i]
# 	        	first, last = get_first_last_atoms(res, i)
# 	        	dist = get_distance(atoms, first, last)
# 	        	if getScore(dist, found_dist) < minimal_score:
# 	        		minimal_score = getScore(dist, found_dist)
# 	        		best_part = res
#     	except Exception:
#     		pass
#     print(best_part, minimal_score)