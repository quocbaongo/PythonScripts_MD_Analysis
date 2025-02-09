import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio.PDB.PDBParser import PDBParser
import argparse
import json
import os

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f

if __name__ == "__main__":


	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Detecting type of aromatic interaction per time step")
	parser.add_argument("--input", type=validate_file, help="The file with the extension '_per_step.json' generated from executing 'Detecting_aromatic_interactions.py' script", required=True)
	parser.add_argument('--dt', type=int, help='Time step (ps)')
	
	args = parser.parse_args()

	# Open and read the JSON file
	with open(f"{args.input}", "r") as file:
		SINAPs_aromatic_per_time_step = json.load(file)
		
	#print(SINAPs_aromatic_per_time_step)
	
	# Pi-Stacking
	SINAPs_pi_stacking=[]
	SINAPs_T_shape=[]
	SINAPs_L_shape=[]
	
	for key, value in SINAPs_aromatic_per_time_step.items():
		if value == "Pi-Stacking":
			SINAPs_pi_stacking.append(1)
			SINAPs_T_shape.append(0)
			SINAPs_L_shape.append(0)
			
		elif value == "T-Shape":
			SINAPs_pi_stacking.append(0)
			SINAPs_T_shape.append(1)
			SINAPs_L_shape.append(0)		
	
		elif value == "L-Shape":
			SINAPs_pi_stacking.append(0)
			SINAPs_T_shape.append(0)
			SINAPs_L_shape.append(1)
			
		elif value == None:
			SINAPs_pi_stacking.append(0)
			SINAPs_T_shape.append(0)
			SINAPs_L_shape.append(0)			


	print(sum(SINAPs_pi_stacking))
	print(sum(SINAPs_T_shape))
	print(sum(SINAPs_L_shape))


	# Plotting
	
	# Set up xtick label
	TimeStep=[(float(i) * float(args.dt)) / 1000 for i in range(len(SINAPs_pi_stacking))]
		
	fig,axes = plt.subplots(1, 3, figsize=(17, 5), layout='constrained',sharey=True)

	# Pi-stacking
	axes[0].plot(TimeStep, SINAPs_pi_stacking, color='blue')
	axes[0].text(0.40, 0.90, "Pi-stacking", transform = axes[0].transAxes, fontsize=20)
	axes[0].set_ylabel("Number of contact",size=20, labelpad=20)
	axes[0].set_xlabel("Time (ns)",size=20, labelpad=20)
	
	# T-shape
	axes[1].plot(TimeStep,SINAPs_T_shape, color='red')
	axes[1].text(0.40, 0.90, "T-shape", transform = axes[1].transAxes, fontsize=20)
	axes[1].set_xlabel("Time (ns)",size=20, labelpad=20)
		
	# L-shape
	axes[2].plot(TimeStep,SINAPs_L_shape, color='green')
	axes[2].text(0.40, 0.90, "L-shape", transform = axes[2].transAxes, fontsize=20)
	axes[2].set_xlabel("Time (ns)",size=20, labelpad=20)
	
	# Set y axis limit
	axes[0].set_ylim([-1,1.5])
	axes[1].set_ylim([-1,1.5])
	axes[2].set_ylim([-1,1.5])
	
	axes[0].set_yticks([0,1])
	
	# Save figure
	fig.tight_layout()
	plt.savefig(f"Aromatic_interactions_evolution.png",dpi=600,bbox_inches="tight")
	


