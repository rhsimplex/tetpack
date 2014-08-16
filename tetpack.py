import pymatgen as pm
import numpy as np
from itertools import combinations

def n_nearest_neighbors(structure, n, max_dist = 5.0):
    t = structure.get_all_neighbors(max_dist)
    for neighbors in t:
        neighbors.sort(key = lambda d: d[-1])
    return [map(lambda d: d[0], neighbors[:n])  for neighbors in t]

def simplex_volume(sites):
    return abs(np.dot((sites[0].coords - sites[3].coords), np.cross(sites[1].coords-sites[3].coords, sites[2].coords-sites[3].coords) ))/6.0

def get_tets(structure, site, neighbors):    
    for face in combinations(neighbors, 3):
        edge_distances = [
                site.distance_from_point(face[0].coords),
                site.distance_from_point(face[1].coords),
                site.distance_from_point(face[2].coords),
                face[0].distance_from_point(face[1].coords),
                face[0].distance_from_point(face[2].coords),
                face[1].distance_from_point(face[2].coords)
                ]
        print np.std(edge_distances), simplex_volume( (site,) + face )

