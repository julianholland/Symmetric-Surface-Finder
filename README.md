# Symmetric Surface Finder (SSF)

This code is meant to find symmetric surfaces for any given periodic crystal and Miller index. 

## The Theory

To find a symmetric surface we create all possible cuts for a given Miller index and periodic crystal (automatically reduced to a primitive cell) within a single layer (where a layer represents the minimum size surface slab that contains sufficient information to replicate the bulk should it be repeated infinitely perpendicular to the direction that it was cut). For each cut a new slab is generated which is at least two layers thick. We then cut on the side opposite to our chose cut and check for symmetry operations using spglib. Should any operation of our surface have a "-z" operator present then our surface is symmetric and we stop there and create our surface.

## Use
We have provided several examples of crystals in the ".cif" format in the input files folder. However, in principle, any periodic crystal in a file format that ASE can read will work. Edit the "bulk=" line in term.py to point to your file.

then simply call the following in your terminal `python term.py`

## Package and Version Requirements

Here are the packages required and their respective versions that we have used to test this code
- python (v3.11.5)
- ase (v3.22.1)
- numpy (v1.26.2
- spglib (v2.2.0)
- tqdm (v4.66.1)


## Acknowledgement

Tom Demeyere was instrumental in getting this code working although he may not approve of some of the code quirks.
