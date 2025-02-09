#!/usr/bin/env python
#
# Calculate residue interaction network for frames in an MD trajectory and
# determine betweenness centrality and average shortest path for residues in
# the frames
#
# Script distributed under GNU GPL 3.0
#
# Author: David Brown
# Date: 17-11-2016

from lib.cli import CLI
from lib.utils import Logger
from lib.trajectory import load_trajectory, calc_distance

import numpy as np
import networkx as nx

import os, sys, argparse, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import multiprocessing as mp
from multiprocessing import Pool

import h5py

def construct_graph(frame, ligands=None, prefix="frame", threshold=6.7, save_graph=True):
    atom_filter = "(name CB and protein) or (name CA and resname GLY)"
    if ligands:
        ligands = ligands.split(",")

        for ligand in ligands:
            arr = ligand.split(":")
            atom_filter += " or (name %s and resname %s)" % (arr[1], arr[0])

    atoms = frame.topology.select(atom_filter)

    nodes_range = len(atoms)

    nodes = range(0, len(atoms))
    edges = []

    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            dist = calc_distance(frame, atoms[i], atoms[j]) * 10
            if dist < threshold:
                edges.append((i, j))

    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)

    if save_graph:
        nx.write_gml(protein_graph, "%s_graph.gml" % prefix)
        nx.write_graphml(protein_graph, "%s_graph.graphml" % prefix)

    return protein_graph


def calc_shortest_paths(frame, traj_name, args):


    prefix = "%s_%d" % (".".join(traj_name.split(".")[:-1]), frame.time[0])

    pg = construct_graph(frame, args.ligands, prefix, args.threshold, args.discard_graphs)
            
    calc_shortest_path(pg, prefix, args.generate_plots, args.xmgrace)
            
def calc_shortest_path(protein_graph, prefix, generate_plots=True, xmgrace=False):

    num_nodes = len(protein_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)

    path_dict = dict(nx.all_pairs_shortest_path_length(protein_graph))
    dj_path_matrix = np.zeros((num_nodes, num_nodes))

    for i in range(num_nodes):
        for j in range(num_nodes):
            try:
                dj_path_matrix[i,j] = path_dict[i][j]
            except KeyError as ke:
                raise nx.exception.NetworkXNoPath("\nERROR::type=orphan_node:message=No link between %d and %d:exception=%s\n" % (i, j, str(ke)))

    # Write dictionary to hdf5
    with h5py.File(f"{prefix}_L.h5", "w") as hdf:
        hdf.create_dataset(f"{prefix}", data=dj_path_matrix)

    avg_L_per_node = np.sum(dj_path_matrix, axis=0)/(num_nodes - 1)

    if generate_plots:
        plt.plot(nodes_axis, avg_L_per_node)
        plt.title("%s L" % prefix, fontsize=18)
        plt.xlabel('Node Indices', fontsize=16)
        plt.ylabel('L', fontsize=16)
        plt.savefig("%s_L.png" % prefix, dpi=300, bbox_inches='tight',format="png")
        plt.close()

    avg_L_per_node = avg_L_per_node.reshape(1, num_nodes)
    
    # Write dictionary to hdf5
    with h5py.File(f"{prefix}_avg_L.h5", "w") as hdf:
        hdf.create_dataset(f"{prefix}", data=avg_L_per_node)    

    # if xmgrace:
    #     dat2xmgrace(avg_L_per_node, prefix, "L", traj=traj)

    return dj_path_matrix


def calc_centralities(frame, traj_name, args):


    prefix = "%s_%d" % (".".join(traj_name.split(".")[:-1]), frame.time[0])
    
    pg = construct_graph(frame, args.ligands, prefix, args.threshold, args.discard_graphs)

    calc_BC(pg, prefix, args.generate_plots)


def calc_BC(protein_graph, prefix, generate_plots=True):
    bc = nx.betweenness_centrality(protein_graph, normalized=False)
    bc = np.asarray(list(bc.values()))

    num_nodes = len(protein_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)

    if generate_plots:
        plt.plot(nodes_axis, bc)
        plt.title("%s BC" % prefix, fontsize=18)
        plt.xlabel('Node Indices', fontsize=16)
        plt.ylabel('BC', fontsize=16)
        plt.savefig("%s_BC.png" % prefix, dpi=300, bbox_inches='tight',format="png")
        plt.close()

    bc = bc.reshape(1, num_nodes)
    
    # Write dictionary to hdf5
    with h5py.File(f"{prefix}_bc.h5", "w") as hdf:
        hdf.create_dataset(f"{prefix}", data=bc)     
   
    return bc

def main(args):
    if not args.calc_BC and not args.calc_L:
        log.error("At least one of the --calc-BC or --calc-L flags must be set.")
        sys.exit(1)

    global traj
    traj_name = os.path.basename(args.trajectory)
    traj, total_frames = load_trajectory(args.trajectory, args.topology, args.step, args.lazy_load)


    # Multiprocessing
    pool = mp.Pool(mp.cpu_count())
    
    if args.calc_BC:
        # Calculate betweeness centrality
        pool.starmap(calc_centralities, [(frame, traj_name, args) for current, frame in enumerate(traj)]) 

    if args.calc_L:
        # Calculate shortest paths
        pool.starmap(calc_shortest_paths, [(frame, traj_name, args) for current, frame in enumerate(traj)]) 
        
log = Logger()
traj = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Topology PDB file (required if trajectory does not contain topology information)")
    parser.add_argument("--ligands", help="Specify any ligands that should be included in the network", default=None)
    parser.add_argument("--threshold", help="Maximum distance threshold in Angstroms when constructing graph (default: 6.7)", default=6.7, type=float)
    parser.add_argument("--step", help="Size of step when iterating through trajectory frames", default=1, type=int)
    parser.add_argument("--generate-plots", help="Generate figures/plots", action='store_true', default=False)
    parser.add_argument("--calc-L", help="Calculate delta L", action='store_true', default=False)
    parser.add_argument("--calc-BC", help="Calculate delta BC", action='store_true', default=False)
    parser.add_argument("--discard-graphs", help="Discard calculated networks when complete (default: save networks in graphml and gml formats)", action='store_false', default=True)
    parser.add_argument("--lazy-load", help="Read frames as they are needed (memory efficient - use for big trajectories)", action='store_true', default=False)
    parser.add_argument("--xmgrace", help="Generate xmgrace compatible format", action='store_true', default=False)

    CLI(parser, main, log)
