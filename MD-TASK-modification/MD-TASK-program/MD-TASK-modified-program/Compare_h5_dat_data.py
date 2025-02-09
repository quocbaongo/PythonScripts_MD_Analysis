#!/usr/bin/env python
import numpy as np
import argparse
import h5py
import os, sys

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description = "Python script to compare the result obtained from original MD-TASK calc_network.py \
							and the modified calc_network.py allowing parallel processing")
	parser.add_argument("--input_h5", help="Input .h5 file generated from modified calc_network.py script", required=True)					
	parser.add_argument("--input_dat", help="Input .dat file generated from original calc_network.py script", required=True)
	parser.add_argument("--calc_L", help="Comparing _L files", action='store_true', default=False)
	parser.add_argument("--calc_BC", help="Comparing _bc files", action='store_true', default=False)
	parser.add_argument("--calc_avg_L", help="Comparing _avg_L files", action='store_true', default=False)
    
	args = parser.parse_args()
	
	
	if not args.calc_BC and not args.calc_L and not args.calc_avg_L:
		print("Either the --calc_BC, --calc_L or --calc_avg_L flags must be set.")
		sys.exit(1)
	elif (args.calc_BC and args.calc_L) or (args.calc_BC and args.calc_avg_L) or (args.calc_L and args.calc_avg_L):
		print("Can not set two flags simultaneously.")
		sys.exit(1)
		
	# Open .dat file
	Dat_Array=np.loadtxt(args.input_dat)
	# Open .h5 file
	with h5py.File(f"{args.input_h5}", "r") as file:
		H5_ArrayKey=os.path.basename(f"{args.input_h5}").split("_")
		H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
			
		H5_Array=np.array(file[H5_ArrayKey])


	if args.calc_BC or args.calc_avg_L:
		# Compare two arrays
		if np.array_equal(Dat_Array, H5_Array[0]):
			print(f"Data from {args.input_h5} and {args.input_dat} are similar")
		else:
			print(f"Data from {args.input_h5} and {args.input_dat} are different")
	elif args.calc_L:
		# Compare two arrays
		if np.array_equal(Dat_Array, H5_Array):
			print(f"Data from {args.input_h5} and {args.input_dat} are similar")
		else:
			print(f"Data from {args.input_h5} and {args.input_dat} are different")




