import pymatgen as pm
import numpy as np
from itertools import combinations, permutations
from rotation_matrix import *
from scipy.optimize import minimize
from sys import maxint

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

class tetrahedron:
    def __init__(self, tetsites):
        self.center = tet_center(tetsites)
        self.coordinates = np.array([site.coords for site in tetsites])
        #retular tetrahedron of unit volume to align with tetrahedral cell from structure
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
        min_angle = minimize(error, 0.0, bounds=(0.0, 2*pi)).x

        #rotate coords
        R_axis_angle(R, centered_tet[0], min_angle)
        best_fit_tet = R.dot(aligned_platonic_tet.T).T
        #un-center tet and store coordinates
        self.fit_regular_tetrahedron = best_fit_tet + self.center
        
    def add_to_structure(self, structure):
        #adds a tet as a methane molecule to a sturcture
        structure.append('C',self.center, coords_are_cartesian = True)
        for coord in self.fit_regular_tetrahedron:
            structure.append('H', coord, coords_are_cartesian = True)
    

