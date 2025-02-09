import pytraj as pt
import numpy as np
from itertools import combinations
import argparse
from Bio.PDB import *
import glob
import os
import json
from Bio.PDB.PDBParser import PDBParser

# VARIABLES
BACKBONE = ["CA", "C", "H", "O", "N"]
RESIDUES = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HIE","HID","HIP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
ARO_RES = ["TYR", "PHE", "TRP"]
parameter_aro_distance = 7

# AROMATIC INTERACTIONS
P_dist_min = 3.0
P_dist_max = 5.0
P_angle_max = 30
TS_dist_min = 4.5
TS_dist_max = 7.0

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f

def TRAJ_loader(traj, parm):
	trajout = pt.load(traj, parm)
	return trajout

def SINAPs_aro_aro(traj_aro, traj):

    ########################################~

    summary = {"Frames":0, "Pi-Stacking":0, "T-Shape":0, "L-Shape":0}
    aro_per_time_step = {}

    ########################################~

    for combination_aro in combinations(traj_aro, 2):

        # Distances
        temp_distance = pt.distance(traj, "{0} {1}".format(combination_aro[0], combination_aro[1]))
        #print(temp_distance[0])
        
        if temp_distance.any() == None:
            continue

        temp_distance_min = min(temp_distance)

        if temp_distance_min > float(TS_dist_max):
            continue
        else:
            temp_corrplane1 = pt.vector.corrplane(traj, combination_aro[0])
            temp_corrplane2 = pt.vector.corrplane(traj, combination_aro[1])
            temp_centers_vector = pt.vector.vector(traj, "{0} {1}".format(combination_aro[0], combination_aro[1]))

            # Planar angle
            temp_planar_angle = []
            for i,j in zip(temp_corrplane1, temp_corrplane2):
                temp_cos = np.dot(i, j) / (np.linalg.norm(i) * np.linalg.norm(j))
                temp_angle = np.rad2deg(np.arccos(np.clip(temp_cos, -1, 1)))
                if temp_angle > 90:
                    temp_planar_angle.append(180 - temp_angle)
                else:
                    temp_planar_angle.append(temp_angle)

            # Orientation angle
            temp_orientation_angle = []
            for i, j in zip(temp_corrplane1, temp_centers_vector):
                temp_cos = np.dot(i, j) / (np.linalg.norm(i) * np.linalg.norm(j))
                temp_angle = np.rad2deg(np.arccos(np.clip(temp_cos, -1, 1)))
                if temp_angle > 90:
                    temp_orientation_angle.append(180 - temp_angle)
                else:
                    temp_orientation_angle.append(temp_angle)

            # Detecting aromatic interactions type and write into dictionary
            # Creating 
            for frame in range(len(temp_distance)):
                aro_per_time_step.update({frame: ""})          
            
            for frame, distance_frame, planar_frame, orientation_frame in zip(range(len(temp_distance)), temp_distance, temp_planar_angle, temp_orientation_angle):            
                # Pi-Stacking
                if 0 <= planar_frame < 20:
                    if P_dist_min < distance_frame < P_dist_max:
                        if 0 <= orientation_frame < P_angle_max:
                            summary["Frames"] += 1
                            summary["Pi-Stacking"] += 1
                            aro_per_time_step[frame] = "Pi-Stacking"
                        else:
                            summary["Frames"] += 1
                            aro_per_time_step[frame] = None
                    else:
                        summary["Frames"] += 1
                        aro_per_time_step[frame] = None

                elif 60 < planar_frame <= 90:
                    if TS_dist_min < distance_frame < TS_dist_max:
                        if 60 <= orientation_frame < 90:
                            summary["Frames"] += 1
                            summary["T-Shape"] += 1
                            aro_per_time_step[frame] = "T-Shape"
                        elif 30 <= orientation_frame < 60:
                            summary["Frames"] += 1
                            summary["L-Shape"] += 1
                            aro_per_time_step[frame] = "L-Shape"
                        elif 0 <= orientation_frame < 30:
                            summary["Frames"] += 1
                            summary["T-Shape"] += 1
                            aro_per_time_step[frame] = "T-Shape"
                        else:
                            summary["Frames"] += 1
                            aro_per_time_step[frame] = None
                    else:
                        summary["Frames"] += 1
                        aro_per_time_step[frame] = None

                else:
                    summary["Frames"] += 1
                    aro_per_time_step[frame] = None

    return summary, aro_per_time_step

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Detecting potential aromatic interactions between two specified amino acids over the course of simulation")
	parser.add_argument("-traj", type=validate_file, help="Trajectory: xtc trr cpt", required=True)
	parser.add_argument("-tpr", type=validate_file, help="Structure+mass(db): pdb", required=True)	
	parser.add_argument("-residue1", help="Specifying the id of first residue following the amino acid ids specified in flag -tpr \
					(example: TYR140A -> amino acid Tyr at amino acid id of 140 and chain id of A). \
     					Note: amino acid triple code must be written in capital", required=True)
	parser.add_argument("-residue2", help="Specifying the id of second residue following the amino acid ids specified in flag -tpr \
					(example: TYR140A -> amino acid Tyr at amino acid id of 140 and chain id of A). \
     					Note: amino acid triple code must be written in capital", required=True)
	args = parser.parse_args()


	################################################################
	################## Trajectory and topology #####################
	################################################################		

	trajectory = TRAJ_loader(args.traj, args.tpr)
	
	trajectory_residues = np.array([i.name for i in trajectory.top.residues])

	trajectory_aro = []

	# constructing resid resname database
	ResIDResname={}
	
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure('struct', args.tpr)  
	
	count = 0
	for model in structure:
		for chain in model:
			for residue in chain:
				ResIDResname.update({f"{residue.resname}{residue.id[1]}{chain.id}": count})
				count += 1
				
	# Aromatic interactions to search
	idx1=ResIDResname[args.residue1]
	idx2=ResIDResname[args.residue2]

	# TYROSINE
	if args.residue1[:3].upper() == "TYR":
		trajectory_aro.append(":{0}@CG,CD1,CD2,CE1,CE2,CZ".format(idx1 + 1))

	if args.residue2[:3].upper() == "TYR":
		trajectory_aro.append(":{0}@CG,CD1,CD2,CE1,CE2,CZ".format(idx2 + 1))

	# PHENYLALANINE
	if args.residue1[:3].upper() == "PHE":
		trajectory_aro.append(":{0}@CG,CD1,CD2,CE1,CE2,CZ".format(idx1 + 1))

	if args.residue2[:3].upper() == "PHE":
		trajectory_aro.append(":{0}@CG,CD1,CD2,CE1,CE2,CZ".format(idx2 + 1))

	# HISTIDINE
	if args.residue1[:3].upper() == "HIS":
		trajectory_aro.append(":{0}@CG,ND1,CD2,CE1,NE2".format(idx1 + 1))

	if args.residue2[:3].upper() == "HIS":
		trajectory_aro.append(":{0}@CG,ND1,CD2,CE1,NE2".format(idx2 + 1))

	# TRYPTOPHAN
	if args.residue1[:3].upper() == "TRP":
		trajectory_aro.append(":{0}@CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2".format(idx1 + 1))

	if args.residue2[:3].upper() == "TRP":
		trajectory_aro.append(":{0}@CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2".format(idx2 + 1))

	# Running
	SINAPs_aromatic_summary, SINAPs_aromatic_per_time_step =SINAPs_aro_aro(trajectory_aro, trajectory)
	
	# Write JSON object to file
	with open(f"{args.residue1}_{args.residue2}_summary.json", "w") as outfile:
		json.dump(SINAPs_aromatic_summary, outfile)
		
	with open(f"{args.residue1}_{args.residue2}_per_step.json", "w") as outfile:
		json.dump(SINAPs_aromatic_per_time_step, outfile)       
        
        
        
        
        
        
        
        
        
        
        
          
