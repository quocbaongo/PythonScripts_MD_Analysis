#!/usr/bin/env python
#
# Calculate the change in the valuses (BC or L) of each residue in a protein over the
# course of an MD simulation
#
# Script distributed under GNU GPL 3.0
#
# Author: David Brown
# Date: 17-11-2016

from natsort import natsorted

from lib.cli import CLI
from lib.utils import Logger
from lib.strategies import normalization

import numpy as np

import os, sys, argparse, matplotlib

import h5py

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import multiprocessing as mp
from multiprocessing import Pool


def calc_delta_L_seq(reference_file, alternative_files, normalizer, result_dir, generate_plots=False):
    
    # Open reference .h5 file
    with h5py.File(reference_file, "r") as file:
        H5_ArrayKey=os.path.basename(reference_file).split("_")
        H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
        reference=np.array(file[H5_ArrayKey])

    num_nodes = reference.shape[0]

    label = normalizer.get_label()

    alternatives = natsorted(alternative_files)
    
    log.info("Calculating %s for %d networks...\n" % (label, len(alternatives)))

    for i, alternative in enumerate(alternatives):
        log.info("Calculating %s (%d/%d)\r" % (label, i + 1, len(alternatives)))

        title = f"{result_dir}/{os.path.basename(alternative).split('.')[0]}"
                
        # Open alternative .h5 file
        with h5py.File(alternative, "r") as file:
            H5_ArrayKey=os.path.basename(alternative).split("_")
            H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
            alternative=np.array(file[H5_ArrayKey])
        
        difference = alternative - reference
        difference = normalizer.normalize(difference, reference)
        prefix = "%s_%s_delta_%s" % (title, normalizer.get_prefix(), normalizer.matrix_type)
        
        
        # Save output to .h5 file
        np.savetxt("%s.dat" % prefix, difference)

        if generate_plots:
            node_axis = range(1, num_nodes + 1)
            plt.plot(node_axis, difference)
            plt.axhline(0, color='black')
            plt.title("%s %s" % (title, label), fontsize=18)
            plt.xlabel('Residue Numbers', fontsize=16)
            plt.ylabel(label, fontsize=16)
            plt.savefig("%s.png" % prefix, dpi=300, bbox_inches="tight")
            plt.close()

    log.info("\n")

def calc_delta_BC_seq(reference_file, alternative_files, normalizer, result_dir, generate_plots=False):
    
    # Open reference .h5 file
    with h5py.File(reference_file, "r") as file:
        H5_ArrayKey=os.path.basename(reference_file).split("_")
        H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
        reference=np.array(file[H5_ArrayKey])

    num_nodes = reference[0].shape[0]

    label = normalizer.get_label()

    alternatives = natsorted(alternative_files)
    
    log.info("Calculating %s for %d networks...\n" % (label, len(alternatives)))

    for i, alternative in enumerate(alternatives):
        log.info("Calculating %s (%d/%d)\r" % (label, i + 1, len(alternatives)))

        title = f"{result_dir}/{os.path.basename(alternative).split('.')[0]}"
                
        # Open alternative .h5 file
        with h5py.File(alternative, "r") as file:
            H5_ArrayKey=os.path.basename(alternative).split("_")
            H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
            alternative=np.array(file[H5_ArrayKey])
        
        difference = alternative[0] - reference[0]
        difference = normalizer.normalize(difference, reference[0])
        prefix = "%s_%s_delta_%s" % (title, normalizer.get_prefix(), normalizer.matrix_type)
        
        
        # Save output to .h5 file
        np.savetxt("%s.dat" % prefix, difference)

        if generate_plots:
            node_axis = range(1, num_nodes + 1)
            plt.plot(node_axis, difference)
            plt.axhline(0, color='black')
            plt.title("%s %s" % (title, label), fontsize=18)
            plt.xlabel('Residue Numbers', fontsize=16)
            plt.ylabel(label, fontsize=16)
            plt.savefig("%s.png" % prefix, dpi=300, bbox_inches="tight")
            plt.close()

    log.info("\n")

def calc_delta_L_single_file(reference_file, alternative_file, normalizer, result_dir, generate_plots=False):
    
    # Open reference .h5 file
    with h5py.File(reference_file, "r") as file:
        H5_ArrayKey=os.path.basename(reference_file).split("_")
        H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
        reference=np.array(file[H5_ArrayKey])

    num_nodes = reference.shape[0]

    label = normalizer.get_label()

    title = f"{result_dir}/{os.path.basename(alternative_file).split('.')[0]}"
                
    # Open alternative .h5 file
    with h5py.File(alternative_file, "r") as file:
        H5_ArrayKey=os.path.basename(alternative_file).split("_")
        H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
        alternative=np.array(file[H5_ArrayKey])
        
    difference = alternative - reference
    difference = normalizer.normalize(difference, reference)
    prefix = "%s_%s_delta_%s" % (title, normalizer.get_prefix(), normalizer.matrix_type)
        
        
    # Save output to .h5 file
    np.savetxt("%s.dat" % prefix, difference)

    if generate_plots:
        node_axis = range(1, num_nodes + 1)
        plt.plot(node_axis, difference)
        plt.axhline(0, color='black')
        plt.title("%s %s" % (title, label), fontsize=18)
        plt.xlabel('Residue Numbers', fontsize=16)
        plt.ylabel(label, fontsize=16)
        plt.savefig("%s.png" % prefix, dpi=300, bbox_inches="tight")
        plt.close()

def calc_delta_BC_single_file(reference_file, alternative_file, normalizer, result_dir, generate_plots=False):
    
    # Open reference .h5 file
    with h5py.File(reference_file, "r") as file:
        H5_ArrayKey=os.path.basename(reference_file).split("_")
        H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
        reference=np.array(file[H5_ArrayKey])

    num_nodes = reference[0].shape[0]

    label = normalizer.get_label()

    title = f"{result_dir}/{os.path.basename(alternative_file).split('.')[0]}"
                
    # Open alternative .h5 file
    with h5py.File(alternative_file, "r") as file:
        H5_ArrayKey=os.path.basename(alternative_file).split("_")
        H5_ArrayKey=f"{H5_ArrayKey[0]}_{H5_ArrayKey[1]}"
        alternative=np.array(file[H5_ArrayKey])
        
    difference = alternative[0] - reference[0]
    difference = normalizer.normalize(difference, reference[0])
    prefix = "%s_%s_delta_%s" % (title, normalizer.get_prefix(), normalizer.matrix_type)
        
        
    # Save output to .h5 file
    np.savetxt("%s.dat" % prefix, difference)

    if generate_plots:
        node_axis = range(1, num_nodes + 1)
        plt.plot(node_axis, difference)
        plt.axhline(0, color='black')
        plt.title("%s %s" % (title, label), fontsize=18)
        plt.xlabel('Residue Numbers', fontsize=16)
        plt.ylabel(label, fontsize=16)
        plt.savefig("%s.png" % prefix, dpi=300, bbox_inches="tight")
        plt.close()


def get_normalizer(matrix_type, mode):
    type_map = {
        "L": "standard",
        "BC": "plusone"
    }

    if matrix_type not in type_map:
        log.info("ERROR: --matrix-type must be specified and must be either 'BC' or 'L'")
        sys.exit(1)

    normalization_mode = mode if mode else type_map[matrix_type]

    return getattr(normalization, normalization_mode, "standard")(matrix_type)


def main(args):

    if args.mode == "Seq":
        if args.matrix_type == "L":
            if args.normalize:
                normalizer = get_normalizer(args.matrix_type, args.normalization_mode)
            else:
                normalizer = normalization.none(args.matrix_type)

            calc_delta_L(args.reference, args.alternatives, normalizer, args.directory, args.generate_plots)
        
        elif args.matrix_type == "BC":
            if args.normalize:
                normalizer = get_normalizer(args.matrix_type, args.normalization_mode)
            else:
                normalizer = normalization.none(args.matrix_type)

            calc_delta_BC(args.reference, args.alternatives, normalizer, args.directory, args.generate_plots)

    elif args.mode == "Multi":
        # Multiprocessing
        pool = mp.Pool(mp.cpu_count())
        
        if args.matrix_type == "L":
            if args.normalize:
                normalizer = get_normalizer(args.matrix_type, args.normalization_mode)
            else:
                normalizer = normalization.none(args.matrix_type)

            pool.starmap(calc_delta_L_single_file, [(args.reference, h5_file, normalizer, args.directory, args.generate_plots) for h5_file in args.alternatives])

        elif args.matrix_type == "BC":
            if args.normalize:
                normalizer = get_normalizer(args.matrix_type, args.normalization_mode)
            else:
                normalizer = normalization.none(args.matrix_type)

            pool.starmap(calc_delta_BC_single_file, [(args.reference, h5_file, normalizer, args.directory, args.generate_plots) for h5_file in args.alternatives])
          
log = Logger()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--matrix-type", help="The type of values in the matrices i.e. BC or L", default=None)
    parser.add_argument("--directory", help="Path to directory, where the output is stored",default=".")
    parser.add_argument("--reference", help="The reference matrix (.dat)")
    parser.add_argument("--alternatives", help="The alternative matrices (.dat)", nargs="*")
    parser.add_argument("--normalize", help="Normalizes the values", action='store_true', default=False)
    parser.add_argument('--normalization-mode', help="Method used to normalize (default for L = standard, default for BC = plusone)", default=None)
    parser.add_argument("--generate-plots", help="Plot results - without setting this flag, no graph will be generated", action='store_true', default=False)
    parser.add_argument("--mode", help="Run the program in sequential or in parallel. Type 'Seq' for sequential run and 'Multi' for parallel run",default="Seq")

    CLI(parser, main, log)
