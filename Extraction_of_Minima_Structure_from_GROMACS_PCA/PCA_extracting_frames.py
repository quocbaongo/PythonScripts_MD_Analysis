import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import distance
import MDAnalysis
import subprocess
import argparse
import os

coords = []

def onclick(event):

	global coords
	
	print("button=%d, x=%d, y=%d, xdata=%f, ydata=%f"%(
        event.button, event.x, event.y, event.xdata, event.ydata))
	
	coords.append((float(format(event.xdata, ".5f")), float(format(event.ydata, ".5f"))))

def validate_file(f):
	if not os.path.exists(f):
		# Argparse uses the ArgumentTypeError to give a rejection message like:
		# error: argument input: x does not exist
		raise argparse.ArgumentTypeError("{0} does not exist".format(f))
	return f


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Interacting with the principle component map to extract favorable conformational states")
	parser.add_argument("--input", type=validate_file, help="Input file (2dproj.xvg)")
	parser.add_argument("--traj", help="Trajectory (traj.xtc)", required=True)
	parser.add_argument("--tpr", help="Topology (topol.tpr)", required=True)
	parser.add_argument("--dt", type=int, help="Time step (ps). Default: 10 ps", required=False)
	

	args = parser.parse_args()
	
	PCA = []
	
	with open(args.input, "r") as f:
		for line in f.readlines():
			if not (line.startswith("#") or line.startswith("@")):
				PCA.append(line)	
				
	PCA = [i.strip().split() for i in PCA]
	
	PCA1 = [round(float(i[0]), 5) for i in PCA]
	PCA2 = [round(float(i[1]), 5) for i in PCA]
	
	# Plotting PCA
	fig = plt.figure()
	plt.plot(PCA1,PCA2,c="k")
	plt.xlabel("PCA1")
	plt.ylabel("PCA2")
	cid = fig.canvas.mpl_connect("button_press_event", onclick)

	plt.show()
	
	# Start extracting frame
	commands = []
	
	for i in coords:
		distances = [distance.euclidean(i, [PCA1, PCA2]) for PCA1, PCA2 in zip(PCA1, PCA2)]	
		index = distances.index(min(distances))
		print(f"Selected coordinate: PCA1={PCA1[index]} and PCA2={PCA2[index]}")
		
		if args.dt == None:
			TimeStep=10
	
			commands.append(f"Suggested Gromacs commands to extract chosen conformations: gmx trjconv -f {args.traj} -s {args.tpr} -dump {int(index) * TimeStep} -o frame_{PCA1[index]}_{PCA2[index]}_{int(index) * 10}ps.pdb")
			
		else:		
			commands.append(f"Suggested Gromacs commands to extract chosen conformations: gmx trjconv -f {args.traj} -s {args.tpr} -dump {int(index) * 10} -o frame_{PCA1[index]}_{PCA2[index]}_{int(index) * 10}ps.pdb")


	for i in commands:
		print(i)
		
	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
