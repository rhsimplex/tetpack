tetpack
=======

Packing tetrahedra based on tetrahedrally-close-packed (TCP) chemical structures

This program attempts to pack tetrahdra starting from real solid-state structures. The algorithm is roughly as follows:

1. Starting from a real structure, add tetrahedra to the tetrahedral voids (subject to some regularity constraint).  Expand structure.
2. Use Ewald summation to relax structure.  Tetrahedra orientations are preserved, but pushed apart.
3. Compress cell axes by small fixed percentage.
4. Check for collisions.
5. If there are collisions, attempt to resolve by randomly walking the tetrahedron in configuration space for N steps (six degrees of freedom: three translational, three rotational).
6. If collisions cannot be resolved, relax structure -- increase cell axes and relax via Ewald summation.
7. Go to step 4.

Currently, this algorithm hits a wall around 10% packing density.  A number of possible improvements are described at the end of this document.

=======
Installation and running:

You will need the pymatgen library (pymatgen.org) installed.  Clone this repository. You should be able to run the tetrahedra extraction and collision function on the default structure (Cu5Zn8) with:

$ python compression.py

Or if you prefer to run in ipython:

In[1]: import compression
In[2]: compression.main()

The compression script runs indefinitely. Resolved packings will be output into a subdirectory named after the original structure (.cif and .csv).

Three structures are included (in pymatgen .mson format):
mp-1368.mson -- Cu5Zn8 "cell-centered packing"
mp-196.mson -- Al5Co2 "face-centered packing"
mp-30784.mson -- Mg2Zn11 "vertex-centered packing"
(Missing: edge-centered packing)

To run with an alternative structure use, e.g.:

$ python compression.py mp-196.mson

or in ipython:

In[1]: compression.main('mp-196.mson')

=======
Missing Features/Future improvements:

1. Currently the compression loop contains no cooling profile -- that is, as the structure becomes more compressed, the allowed random motion/rotation of tetrahedra does not decrease. Neither does the Ewald relaxation step get smaller.
2. Cell axes/angles cannot vary freely in this implementation, but only for simplicity.  Since symmetry is immediately broken by tetrahedra moving around, axes/angles ought to be able to move freely.
3. Compression could be more Metropolis-like: any configuration resulting in collisions is currently rejected, but could be accepted with some probability related to the packing density change.
4. ...

=======
Attributions:

1. Starting from real TCP crystal structures was inspired by my own thesis research (https://www.dropbox.com/s/vte7zpx59gaqdme/thesis.pdf?dl=0).
2. This program relies heavily on the pymatgen (pymatgen.org) library, and the initial structures were pulled from the Materials Project (materialsproject.org).
3. rotation_matrix.py is part of the relax package written by Edward d'Auvergne.
4. ewald.py is part of pymatgen but not included with the default build, so it's included separately here. Written by Shyue Ping Ong and William Davidson Richard.
