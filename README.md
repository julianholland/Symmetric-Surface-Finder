# Symmetric Surface Finder (SSF)

This code is meant to find symmetric surfaces for any given periodic crystal and Miller index. 

## The Theory

To find a symmetric surface we create all possible cuts for a given Miller index and periodic crystal (automatically reduced to a primitive cell) within a single layer (where a layer represents the minimum size surface slab that contains sufficient information to replicate the bulk should it be repeated infinitely perpendicular to the direction that it was cut). For each cut a new slab is generated which is at least two layers thick. We then cut on the side opposite to our chosen cut and check for symmetry operations using spglib. Should any operation of our surface have a "-z" operator present then our surface is symmetric and we stop there and create our surface.

## Use
We have provided several examples of crystals in the ".cif" format in the input files folder. However, in principle, any periodic crystal in a file format that ASE can read will work. Edit the `bulk=` line in term.py to point to your file.

The `minimum_thickness=` parameter can be tuned to any value thickness you desire. Setting it to 0 will give you the smallest possible layer that still contains all information to reproduce the bulk.

`miller_max=` will generate all allowed Miller planes up to the value provided. E.g. `miller_max=2` will produce a list of Miller planes from \[0,0,1\] to \[2,2,2\]. Alternatively, you can provide a list of desired Miller planes in the `allowed_miller_planes=` variable directly.

then simply call the following in your terminal `python term.py`

## Package and Version Requirements

Here are the packages required and their respective versions that we have used to test this code
- python (v3.11.5)
- ase (v3.22.1)
- numpy (v1.26.2)
- spglib (v2.2.0)
- tqdm (v4.66.1)


## Outstanding known bugs/issues
- Unable to reconcile by hand determination or Tom's number of unique cuts (typically get significantly more) this likely means we are missing or generating too many surfaces
- different number of slabs produced depending on minimum thickness

## Features to be implemented
+ Add an option to generate surfaces of multiple thicknesses in a single run for F+M
+ Add a means to check similarity between different miller planes (i.e. are these miller planes the same for this material) see felix's notes on that
+ Add feature to find what atoms are present on each surface


## Acknowledgement

Tom Demeyere was instrumental in getting this code working although he may not approve of some of the code quirks. Dr. Felix Hanke had the original idea of the "-z" symmetry check.
