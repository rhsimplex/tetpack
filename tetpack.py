import pymatgen as pm
import numpy as np
import csv
import copy
import warnings
import ewald
from itertools import combinations, combinations_with_replacement, permutations
from rotation_matrix import *
from scipy.optimize import minimize
from sys import maxint
"""
The tetpack module works like this:

1. Start with a solid-state structure, preferably something "tetrahedral close-packed (tcp)," e.g. gamma brass(Cu5Zn8)
2. Scan up to n nearest neighbors for each site (see function n_nearest_neighbors)
3. Pick out "most regular" tetrahedra formed by each site and any three neighbors (get_tets)
4. Remove duplicates (remove_duplicate_tets)
5. Build regular tetrahedra. They are aligned to best match the tetrahedral cell they are drawn from, which of course is not regular (see the tetrahedron class)
6. Add the perfect tetrahedra to a new periodic structure as methane -- C is the center, each H is a vertex (tetrahedron.add_to_structure)
7. TO BE IMPLEMENTED -- compression, annealing...

There is a convenience method tetrahedron_from_structure that combines steps 1-6. It takes a pymatgen structure object and returns a new structure with tetrahedra (methanes) and a list of tetrahedra objects:

Example:

import pymatgen as pm                                       #import pymatgen
import tetpack                                              #import this module
gamma = pm.read_structure('mp-1368.mson')                   #load the gamma brass structure
tet_str, tet_reg = tetpack.tetrahedra_from_structure(gamma) #generate a new pymatgen structure (periodic) and a list of tetrahedra objects
pm.write_structure(tet_str, 'tet_str.cif')                  #export our methane tetrahedra to a CIF so we can view it in the crystallohraphic structure program of our choice!
"""
def random_configuration(structure, n_tets):
    #randomly distributes n_tets tetrahedra in structure
    pass

def ewald_relaxation(structure, max_steps=3, motion_factor = 1.):
    #attempt to relax structure via computing ewald forces, moving coordinates, etc.
    ions = structure.copy()
    ions.remove_species('H')
    ions.add_oxidation_state_by_element({'C':1})

    es = ewald.EwaldSummation(ions)
    
    for i in range(max_steps):
        print "Total Energy: " + str(es.total_energy)
        disp_vec = motion_factor/10.*es.forces
        ions_pos_current = [x.coords for x in ions.sites] + disp_vec
        ions_updated = ions.copy()
        ions_updated.remove_oxidation_states()
        ions_updated.remove_species('C')
        for pos in ions_pos_current:
            ions_updated.append('C', pos, coords_are_cartesian = True)
        ions_updated.add_oxidation_state_by_element({'C':1})
        es = ewald.EwaldSummation(ions_updated)
        ions = ions_updated.copy()
    ions.remove_oxidation_states()
    hydrogens = structure.copy()
    hydrogens.remove_species('C')
    return_str = hydrogens.copy()
    return_str.remove_species('H')
    for i in range(len(ions.sites)):
        hydrogen_group = hydrogens[4*i:4*i+4]
        t = ions[i].coords - structure[5*i].coords
        return_str.append('C', ions[i].coords, coords_are_cartesian = True)
        for hydrogen in hydrogen_group:
            return_str.append('H', hydrogen.coords + t, coords_are_cartesian = True)
    return return_str

def n_nearest_neighbors(structure, n, max_dist = 5.0):
    t = structure.get_all_neighbors(max_dist)
    for neighbors in t:
        neighbors.sort(key = lambda d: d[-1])
    return [map(lambda d: d[0], neighbors[:n])  for neighbors in t]

def simplex_volume(sites):#triple product evaluation of tetrahedral volume
    if type(sites[0]) == pm.core.sites.PeriodicSite:
        return abs(np.dot((sites[0].coords - sites[3].coords), np.cross(sites[1].coords-sites[3].coords, sites[2].coords-sites[3].coords) ))/6.0
    else:
        return abs(np.dot((sites[0] - sites[3]), np.cross(sites[1]-sites[3], sites[2]-sites[3]) ))/6.0

def get_tets(structure, site, neighbors, stdcutoff = 0.10, volcutoff = 2.0):#relax stdcutoff to 0.50 for more initial tetrahedra
    tets=[]
    for face in combinations(neighbors, 3):
        edge_distances = [
                site.distance_from_point(face[0].coords),
                site.distance_from_point(face[1].coords),
                site.distance_from_point(face[2].coords),
                face[0].distance_from_point(face[1].coords),
                face[0].distance_from_point(face[2].coords),
                face[1].distance_from_point(face[2].coords)
                ]
        dstd = np.std(edge_distances)
        v = simplex_volume( (site,) + face )
        if dstd < stdcutoff and v > volcutoff:
            tets.insert(-1, (site,) + face)            
    return tets

def site_in_cell(site):
    return not any(map(lambda x: round(x,4) < 0. or round(x,4) >= 1., site.frac_coords))

def tets_in_cell(structure, tets):
    #takes a list of tetrahedra and removes tets which have centers outside the unitcell of structure
    dummy_structure = structure.copy()
    dummy_structure.remove_sites(range(len(dummy_structure.sites)))
    in_tets = []
    for tet in tets:
        dummy_structure.append('C', tet.center, coords_are_cartesian = True)
        if site_in_cell(dummy_structure[-1]):
            in_tets.append(tet)
    return in_tets

def tet_supercell(structure, tets, cutoff=2.0):
    #creates a 3x3 supercell, with the original tetrahedra in the middle cell. may overlap edge tetrahedra
    t = []
    supercell = []
    frac_cutoff_a = cutoff/structure.lattice.a
    frac_cutoff_b = cutoff/structure.lattice.b
    frac_cutoff_c = cutoff/structure.lattice.c
    empty_structure = structure.copy()
    empty_structure.remove_sites(range(len(empty_structure.sites)))
    for x in combinations_with_replacement([-1,0,1],3):
        for y in permutations(x,3):
            t.append(y)
    t = np.unique(t)
    t = np.dot(t, structure.lattice_vectors())
    for vec in t:
        for tet in tets:
            new_tet = copy.deepcopy(tet)
            new_tet.translate(vec)
            empty_structure.append('C', new_tet.center, coords_are_cartesian = True)
            if empty_structure.sites[-1].frac_coords[0] < 1. + frac_cutoff_a and\
                    empty_structure.sites[-1].frac_coords[0] >= -frac_cutoff_a and\
                    empty_structure.sites[-1].frac_coords[1] < 1. + frac_cutoff_b and\
                    empty_structure.sites[-1].frac_coords[1] >= -frac_cutoff_b and\
                    empty_structure.sites[-1].frac_coords[2] < 1. + frac_cutoff_c and\
                    empty_structure.sites[-1].frac_coords[2] >= -frac_cutoff_c:
                supercell.append(new_tet)
    return supercell

def get_nearby_tets(tets, supercell_tets, radius=2.0):
    #returns indices of tets in supercell whose centers lie within the radius of the centers of tets
    centers = np.zeros((len(supercell_tets), 3))
    i = 0
    neighbor_indices = []
    for tet in supercell_tets:
        centers[i] = tet.center
        i += 1
    for tet in tets:
        differences = centers - np.tile(tet.center, (len(centers),1))
        norms = np.apply_along_axis(np.linalg.norm, 1, differences)
        neighbor_indices.append(np.where((norms < radius) & (norms > 0.00001) )[0])
    return neighbor_indices
    
def tet_center(tet_sites):
    #returns CARTESIAN COORDINATES
    if type(tet_sites[0]) == pm.core.sites.PeriodicSite:
        return reduce(np.add, [site.coords for site in tet_sites])/4.
    else:
        return reduce(np.add, [site for site in tet_sites])/4.

def remove_duplicate_tets(tet_list):#remove duplicate tets from a list by center of gravity only
    eps = 0.1
    unique_tets = []
    while len(tet_list) > 0:
        temp = tet_list.pop()
        if len(unique_tets) <= 0:
            unique_tets += (temp,)
        elif not any([np.linalg.norm(tet_center(temp) - tet_center(tet)) < eps for tet in unique_tets]):
            unique_tets += (temp,)
    return unique_tets

def tetrahedra_from_structure(structure, n=12, maxdist=5.0, stdcutoff=0.10, volcutoff=2.0):
    neighbors = n_nearest_neighbors(structure, n, maxdist)
    raw_tets = []
    for i in range(len(structure.sites)):
        raw_tets += get_tets(structure, structure.sites[i], neighbors[i], stdcutoff, volcutoff)
    unique_tets = remove_duplicate_tets(raw_tets)
    regular_tetrahedra = [tetrahedron(tet) for tet in unique_tets]
    empty_structure = structure.copy()
    empty_structure.remove_sites(range(len(empty_structure.sites)))
    [tet.add_to_structure(empty_structure) for tet in regular_tetrahedra]
    return empty_structure, tets_in_cell(empty_structure, regular_tetrahedra)

def nearest(tet, neighbors):
#gets nearest neighbor to tet from list of tets
    nearest_tet = neighbors[0]
    distance = 100000.
    for neighbor in neighbors:
        temp_distance = np.linalg.norm(tet.center - neighbor.center)
        if temp_distance < distance:
            distance = temp_distance
            nearest_tet = neighbor
    return nearest_tet

def adjust_axes(structure, a_per, b_per=False, c_per=False, alpha_per=False, beta_per=False, gamma_per=False):
    if not b_per:
        b_per = a_per
    if not c_per:
        c_per = a_per
    
    """
    #Doesn't work well...use apply strain for now
    new_a = (1. + a_per) * structure.lattice.a
    new_b = (1. + b_per) * structure.lattice.b
    new_c = (1. + c_per) * structure.lattice.c
    if alpha_per:    
        new_alpha = (1. + alpha_per) * structure.lattice.alpha
    if beta_per:
        new_beta  = (1. + beta_per) * structure.lattice.beta
    if gamma_per:
        new_gamma = (1. + gamma_per) * structure.lattice.gamma
    """
    centers = structure.copy()
    old_centers = structure.copy()
    vertices = structure.copy()
    centers.remove_species('H')
    old_centers.remove_species('H')
    vertices.remove_species('C')

    tets = [vertices.sites[x:x+4] for x in range(0, len(vertices.sites), 4)]
    
    """
    centers.lattice.a = new_a
    centers.lattice.b = new_b
    centers.lattice.c = new_c
    
    if alpha_per:    
        centers.lattice.alpha = new_alpha 
    if beta_per:
        centers.lattice.beta = new_beta 
    if gamma_per:
        centers.lattice.gamma = new_gamma
    """
    centers.apply_strain([a_per, b_per, c_per])

    new_structure = centers.copy()
    new_structure.remove_species('C')

    for i in range(len(centers.sites)):
        displacement = centers.sites[i].coords - old_centers.sites[i].coords
        new_structure.append('C', centers.sites[i].coords, coords_are_cartesian = True)
        for vertex in tets[i]:
            new_structure.append('H', vertex.coords + displacement, coords_are_cartesian = True)

    return new_structure

def to_challenge_output(structure, output_path):
    with open(output_path, 'wb') as f:
        writer = csv.writer(f, delimiter = " ")
        for basis in structure.lattice_vectors():
            writer.writerow(basis)
        vertices = structure.copy()
        vertices.remove_species('C')
        tets = [vertices.sites[x:x+4] for x in range(0, len(vertices.sites), 4)]
        for tet in tets:
            row_block = np.array([vertex.frac_coords for vertex in tet])
            row = row_block.reshape(1,12)[0]
            writer.writerow(row)

def packing_density(structure):
    #assumes tetrahedra are properly formed, and unit volume etc. Simply divide the number of carbons by the unit cell volume
    return [site.specie.symbol for site in structure.sites].count('C')/structure.volume

def tetrahedron_collision(tet, neighbors):
    #checks if tet collides with any neighbors, using triangle collision detection
    #a modified implementation of triangle-triangle intersection available at http://hal.archives-ouvertes.fr/docs/00/07/21/00/PDF/RR-4488.pdf
    #assume triangles of tet are already canonical order (vertices (abc) are clockwise in sense: det(abc) > 0)
    def bracket_det(a,b,c,d):
        return np.linalg.det(np.concatenate((np.array([a,b,c,d]), np.ones((4,1))),axis = 1))
    if neighbors == []:
        return False
    for face in tet.triangles:
        for neighbor in neighbors:
            for neighbor_face in neighbor.triangles:
                if bracket_det(face[0], face[1], neighbor_face[0], neighbor_face[1]) <= 0.:
                    if bracket_det(face[0], face[2], neighbor_face[2], neighbor_face[0]) <= 0.:
                        return True
    return False

class tetrahedron:
    def __init__(self, tetsites):
        #if structure is all methanes, use alternate constructor
        if(all(np.unique([site.specie.symbol for site in tetsites]) == ['C', 'H'])):
            self.match_methane(tetsites)
            return

        #else fit tets to structure
        self.center = tet_center(tetsites)
        self.coordinates = np.array([site.coords for site in tetsites])
        #regular tetrahedron of unit volume to align with tetrahedral cell from structure
        platonic_tetrahedron = 3**(1/3.)/2.*np.array([   [ 1.0,  1.0,  1.0],
                                            [ 1.0, -1.0, -1.0],
                                            [-1.0,  1.0, -1.0],
                                            [-1.0, -1.0,  1.0]])

        #subtract out the center
        centered_tet = self.coordinates - self.center

        #initialize a rotation matrix
        R = np.zeros([3,3])

        #shuffle input coordinates, just in case there is incoming bias
        np.random.shuffle(centered_tet)
        
        #find the rotation matrix
        R_2vect(R, platonic_tetrahedron[0], centered_tet[0])

        #rotate the unit tetrahedron to align with first point
        aligned_platonic_tet = R.dot(platonic_tetrahedron.T).T

        #rotate around the primary axis (axis from origin to first point) to minimize MSE
        #--------------------------------------------------------------------------------

        #finds the minimum error over all choice of matching:
        #   A       1
        #  / \     / \
        # B---C   2---3
        #
        #the triangles to be matched are exact up to labeling of coordinates. to find the most similarity, we should try matching any order of points (123) (132) (213) ... 
        def min_config(plat_base, coords_base):
            min_err = maxint
            min_config = (0,0,0)
            for comb in permutations(np.arange(3), 3):
                err = sum(np.linalg.norm(plat_base[list(comb)] - coords_base, axis = 1))
                if err < min_err:
                    min_err = err
                    min_config = comb
            return err, comb
        
        #returns sum of errors from alignment per angle. returns the min error over all configurations
        def error(angle):
            Rr = np.zeros([3,3])
            R_axis_angle(Rr, centered_tet[0], angle)
            rotated = Rr.dot(aligned_platonic_tet.T).T
            err, comb = min_config(aligned_platonic_tet[1:], centered_tet[1:])
            return err
            
        #minimize error through rotation
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            min_angle = minimize(error, 0.0, bounds=(0.0, 2*pi)).x

        #rotate coords
        R_axis_angle(R, centered_tet[0], min_angle)
        best_fit_tet = R.dot(aligned_platonic_tet.T).T
        #un-center tet and store coordinates
        self.fit_regular_tetrahedron = best_fit_tet + self.center
        #compute canonical order for triangles (for collision detection)
        self.triangles = self.canonical_triangles()
        
    def add_to_structure(self, structure):
        #adds a tet as a methane molecule to a sturcture IF TETRAHEDRON CENTER IS WITHIN UNIT CELL
        structure.append('C',self.center, coords_are_cartesian = True)
        if(site_in_cell(structure.sites[-1])):
            for coord in self.fit_regular_tetrahedron:
                structure.append('H', coord, coords_are_cartesian = True)
        else:
            structure.remove_sites([len(structure.sites)-1])

    def canonical_triangles(self):
        #gets coordinates in canonical order for each face. important for collision detection
        triangles = []
        for x in combinations(range(4), 3):
            triangles.append(self.fit_regular_tetrahedron[list(x)])
            if np.linalg.det(triangles[-1]) < 0:
                triangles[-1] = np.flipud(triangles[-1])
        return triangles

    def translate(self, vec):
        #translates regular tet and center by adding vec
        self.center += vec
        for coord in self.fit_regular_tetrahedron:
            coord += vec
        for triangle in self.triangles:
            for coord in triangle:
                coord += vec
        self.triangles = self.canonical_triangles()

    def match_methane(self, sites):
        #turns a CH4 into a tetrahedron. does no check for regularity.
        assert(len(sites) == 5)
        regular_tetrahedron = []
        for site in sites:
            if site.specie.symbol == 'C':
                self.center = site.coords
            else:
                regular_tetrahedron.append(site.coords)
        self.fit_regular_tetrahedron = np.array(regular_tetrahedron)
        self.triangles = self.canonical_triangles()

    def regularize(self):
        #tetrahedra will become distorted due to compounding numerical error. this function restores a tet to regularity
        #regular tetrahedron of unit volume to align with tetrahedral cell from structure
        platonic_tetrahedron = 3**(1/3.)/2.*np.array([   [ 1.0,  1.0,  1.0],
                                            [ 1.0, -1.0, -1.0],
                                            [-1.0,  1.0, -1.0],
                                            [-1.0, -1.0,  1.0]])

        #subtract out the center
        centered_tet = self.fit_regular_tetrahedron - self.center

        #initialize a rotation matrix
        R = np.zeros([3,3])

        #shuffle input coordinates, just in case there is incoming bias
        np.random.shuffle(centered_tet)
        
        #find the rotation matrix
        R_2vect(R, platonic_tetrahedron[0], centered_tet[0])

        #rotate the unit tetrahedron to align with first point
        aligned_platonic_tet = R.dot(platonic_tetrahedron.T).T

        #rotate around the primary axis (axis from origin to first point) to minimize MSE
        #--------------------------------------------------------------------------------

        #finds the minimum error over all choice of matching:
        #   A       1
        #  / \     / \
        # B---C   2---3
        #
        #the triangles to be matched are exact up to labeling of coordinates. to find the most similarity, we should try matching any order of points (123) (132) (213) ... 
        def min_config(plat_base, coords_base):
            min_err = maxint
            min_config = (0,0,0)
            for comb in permutations(np.arange(3), 3):
                err = sum(np.linalg.norm(plat_base[list(comb)] - coords_base, axis = 1))
                if err < min_err:
                    min_err = err
                    min_config = comb
            return err, comb
        
        #returns sum of errors from alignment per angle. returns the min error over all configurations
        def error(angle):
            Rr = np.zeros([3,3])
            R_axis_angle(Rr, centered_tet[0], angle)
            rotated = Rr.dot(aligned_platonic_tet.T).T
            err, comb = min_config(aligned_platonic_tet[1:], centered_tet[1:])
            return err
            
        #minimize error through rotation
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            min_angle = minimize(error, 0.0, bounds=(0.0, 2*pi)).x

        #rotate coords
        R_axis_angle(R, centered_tet[0], min_angle)
        best_fit_tet = R.dot(aligned_platonic_tet.T).T
        #un-center tet and store coordinates
        self.fit_regular_tetrahedron = best_fit_tet + self.center
        #compute canonical order for triangles (for collision detection)
        self.triangles = self.canonical_triangles()

    def rotate(self, vertex_index, theta): 
        #subtract out the center
        centered_tet = self.fit_regular_tetrahedron - self.center

        #initialize a rotation matrix
        R = np.zeros([3,3])

        #rotate coords
        R_axis_angle(R, centered_tet[vertex_index], theta)
        rotated_tet = R.dot(centered_tet.T).T
        #un-center tet and store coordinates
        self.fit_regular_tetrahedron = rotated_tet + self.center
        #compute canonical order for triangles (for collision detection)
        self.triangles = self.canonical_triangles()

    def jostle(self, T=1., rotation_only = False):
        #translates and rotates a tetrahedron randomly, with magnitude according to temperature T
        pars = float(T)*np.random.randn(6)
        if not rotation_only:
            self.translate(pars[0:3])
        for i in range(3):
            vertices = np.random.choice(range(4), 3, replace=False)
            self.rotate(vertices[i], 2*np.pi*pars[3+i])
        self.regularize()

class dimer:
    def __init__(self, center, axial_direction, first_vertex_direction):
        pass

class barrelan:
    def __init__(self, center, axial_direction, first_vertex_direction):
        pass
