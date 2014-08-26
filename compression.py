import pymatgen as pm
import numpy as np
import tetpack
import copy
import sys
import os

def main():
    #Filename of starting structure:
    filename = 'mp-1368.mson'

    #Cutoff for tetrahedron distortion -- higher numbers will accept more distorted tetrahedra
    std_dev_cutoff = 0.50

    #Initial factor for increasing all axes
    initial_increase_factor = 2.5

    #Compression increment -- compress by fixed percentage:
    compression_factor = -0.01
    
    #Tetrahedra become distorted due compounding numerical error. Re-regularize every n steps:
    normalization_frequency = 5

    #Save structure every n steps:
    save_frequency = 5

    #Controls how much tetrahedra can jostle during packing
    temp = 150.

    #How many tries to fit a tetrahedra before skipping
    resolution_max = 10000

    #How far a tet can travel randomly
    distance_max = 1.0

    #Initialize sturcture and tetrahedra
    print '\nLoading initial structure...',
    sys.stdout.flush()
    initial_structure = pm.read_structure(filename)
    path = initial_structure.composition.alphabetical_formula.replace(' ', '')
    if not os.path.exists(path):
        os.mkdir(path)
    print initial_structure.composition.alphabetical_formula + ' loaded.'
    print '\nExtracting tetrahedra...',
    sys.stdout.flush()
    tet_str, tet_reg = tetpack.tetrahedra_from_structure(initial_structure, stdcutoff=std_dev_cutoff)
    print str(len(tet_reg)) + ' initial tetrahedra extracted.'

    #Expand structure initally
    print '\nExpanding cell axes by factor of ' + str(initial_increase_factor) + '...',
    sys.stdout.flush()
    current_tet_str = tetpack.adjust_axes(tet_str, initial_increase_factor)
    current_tet_reg = map(tetpack.tetrahedron, [current_tet_str[5*i:5*i+5] for i in range(len(current_tet_str)/5)])
    print 'done: \na = ' + str(current_tet_str.lattice.a) + '\nb = ' + str(current_tet_str.lattice.b) + '\nc = '  + str(current_tet_str.lattice.c) 

    print '\nRelaxing structure via Ewald summation...'
    sys.stdout.flush()
    current_tet_str = tetpack.ewald_relaxation(current_tet_str, max_steps = 5)

    print '\nBeginning compression loop:'

    #Loop until collision
    collision = False
    step = 0
    while(not collision):
        if np.mod(step, normalization_frequency) == 0:
            print 'Normalizing tetrahedra...',
            sys.stdout.flush()
            [tet.regularize() for tet in current_tet_reg]
            print 'done.'
        if np.mod(step, save_frequency) == 0:
            print 'Writing structure...',
            sys.stdout.flush()
            pm.write_structure(current_tet_str, os.path.join(path, str(step) + '.cif'))
            print 'done.'
        current_tet_str, current_tet_reg = compress(current_tet_str, current_tet_reg, compression_factor)
        print 'Step '+ str(step) + ' packing fraction: ' + str(tetpack.packing_density(current_tet_str)) + '...',
        sys.stdout.flush()
        failed = check_and_resolve_collisions(current_tet_str, current_tet_reg, temp, distance_max, resolution_max)
        if failed:
            print 'Relaxing structure...',
            sys.stdout.flush()
            current_tet_str, current_tet_reg = compress(current_tet_str, current_tet_reg, -compression_factor)
            print 'done. Packing fraction: ' + str(tetpack.packing_density(current_tet_str))
            print 'Single-step Ewald relaxation...'
            sys.stdout.flush()
            current_tet_str = tetpack.ewald_relaxation(current_tet_str, max_steps = 1)
            print 'done.'
            failed = False
        step += 1

#Compress structure by fixed percentage
def compress(current_str, current_tet, compression_factor):
    compressed_str = tetpack.adjust_axes(current_str, compression_factor)
    compressed_tet = map(tetpack.tetrahedron, [current_str[5*i:5*i+5] for i in range(len(current_str)/5)])
    return compressed_str, compressed_tet

def check_and_resolve_collisions(current_str, current_tet, temp, distance_max, resolution_max):
    supercell = tetpack.tet_supercell(current_str, current_tet)
    nearby_tet_indices = tetpack.get_nearby_tets(current_tet, supercell, radius=distance_max + 1.)
    coll = [tetpack.tetrahedron_collision(current_tet[i], [supercell[j] for j in nearby_tet_indices[i]]) for i in range(len(current_tet))]
    failed = False
    if any(coll):
        print 'collisions detected!'
        collision_indices = np.where(np.array(coll)==True)[0]
        np.random.shuffle(collision_indices)
        for tet_index in collision_indices:
            original_center = current_tet[tet_index].center
            original_vertices = current_tet[tet_index].fit_regular_tetrahedron
            print 'Resolving collision (tet# ' + str(tet_index) + ')...',
            sys.stdout.flush()
            neighbors = [supercell[j] for j in nearby_tet_indices[tet_index]]
            resolution_iter = 0
            traveled = 0.0 
            while tetpack.tetrahedron_collision(current_tet[tet_index], neighbors):
                if resolution_iter > resolution_max:
                    print 'Couldn\'t resolve after ' + str(resolution_max) + ' tries.  Resetting.'
                    current_tet[tet_index].center = original_center
                    current_tet[tet_index].fit_regular_tetrahedron = original_vertices
                    resolution_iter = 0
                    failed = True
                    break
                if np.linalg.norm(original_center - current_tet[tet_index].center) > distance_max:
                    print 'Tet has moved too far (' + str(np.linalg.norm(original_tet.center - current_tet[tet_index].center) ) + '). Resetting.'
                    current_tet[tet_index].center = original_center
                    current_tet[tet_index].fit_regular_tetrahedron = original_vertices
                    resolution_iter = 0 
                    failed = True
                    break
                #attempt to move tet away from neighbors
                #closest_neighbor = tetpack.nearest(current_tet[tet_index], neighbors)
                #direction = current_tet[tet_index].center - closest_neighbor.center
                #current_tet[tet_index].translate(direction/(10.*np.linalg.norm(direction)))
                #otherwise, random movements
                current_tet[tet_index].jostle(T = temp,rotation_only=True)
                resolution_iter += 1
            print 'resolved.'
    else:
        print 'no collisions detected.'
    return failed

if __name__ == "__main__":
    main()
