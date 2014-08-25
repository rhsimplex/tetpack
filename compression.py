import pymatgen as pm
import numpy as np
import tetpack

def main():
    #Filename of starting structure:
    filename = 'mp-1368.mson'

    #Cutoff for tetrahedron distortion -- higher numbers will accept more distorted tetrahedra
    std_dev_cutoff = 0.50

    #Initial factor for increasing all axes
    initial_increase_factor = 3.0

    #Compression increment -- compress by fixed percentage:
    compression_factor = -0.01

    #Initialize sturcture and tetrahedra
    print 'Loading initial structure...',
    initial_structure = pm.read_structure(filename)
    print initial_structure.composition.alphabetical_formula + ' loaded.'
    print 'Extracting tetrahedra...',
    tet_str, tet_reg = tetpack.tetrahedra_from_structure(initial_structure, stdcutoff=std_dev_cutoff)
    print str(len(tet_reg)) + ' initial tetrahedra extracted.'

    #Expand structure initally
    print 'Expanding cell axes by factor of ' + str(initial_increase_factor) + '...',
    current_tet_str = tetpack.adjust_axes(tet_str, initial_increase_factor)
    current_tet_reg = map(tetpack.tetrahedron, [current_tet_str[5*i:5*i+5] for i in range(len(current_tet_str)/5)])
    print 'done: \na = ' + str(current_tet_str.lattice.a) + '\nb = ' + str(current_tet_str.lattice.b) + '\nc = '  + str(current_tet_str.lattice.c) 

    print '\nBeginning compression loop:\n\n'

    #Loop until collision
    collision = False
    while(not collision):
        current_tet_str, current_tet_reg = compress(current_tet_str, current_tet_reg, compression_factor)
        print 'Packing fraction: ' + str(tetpack.packing_density(current_tet_str)) + '...',
        coll = check_collisions(current_tet_str, current_tet_reg)
        if any(coll):
            print 'collisions detected!'
            print 'Colliding tetrahedra: ' + ', '.join( [str(x) for x in np.where(np.array(coll)==True)])
            return current_tet_str, current_tet_reg
            break
        else:
            print 'no collisions detected.'

#Compress structure by fixed percentage
def compress(current_str, current_tet, compression_factor):
    compressed_str = tetpack.adjust_axes(current_str, compression_factor)
    compressed_tet = map(tetpack.tetrahedron, [current_str[5*i:5*i+5] for i in range(len(current_str)/5)])
    return compressed_str, compressed_tet

#Get vector of bools: True = tet at index i is colliding with some other tet(s), False = no collisions
def check_collisions(current_str, current_tet):
    supercell = tetpack.tet_supercell(current_str, current_tet)
    nearby_tet_indices = tetpack.get_nearby_tets(current_tet, supercell)
    coll = [tetpack.tetrahedron_collision(current_tet[i], [supercell[j] for j in nearby_tet_indices[i]]) for i in range(len(current_tet))]
    return coll

if __name__ == "__main__":
    main()
