from glob import glob
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import surface
from ase.io import read, write
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from spglib import find_primitive, get_symmetry
from tqdm import tqdm

"""
Description: A code to find all possible symmetric surfaces for a given bulk (periodic) crystal and a given miller plane
Author: Julian Holland, Email: Julian.Holland@soton.ac.uk
v0.2.2 (2023/12/18)

Change log
+ Added minimum thickness requirement
+ Changed outpu

Outstanding known bugs/issues
- Unable to reconcile by hand determination or Tom's number of unique cuts (typically get significantly more) this likely means we are missing or generating too many surfaces
- low minimum thicknesses (below 5) sometimes produce slabs of thicknesses slightly smaller than desired

Features to be implemented
+ Add an option to generate surfaces of multiple thicknesses
+ Add a means to check similarity between different miller planes (i.e. are these miller planes the same for this material)
+ Add feature to find what atoms are present on each surface
"""


def find_smallest_sym_surface(surf, surface_layers, vacuum, position_accuracy):
    """For any given surface this function will produce a symmetric surface should one exist"""
    cell = surf.cell  # find cell parameters
    frac_cell_value = (
        ((cell[2, 2] - 2 * vacuum) * surface_layers) / (surface_layers + 1)
    ) + vacuum  # find a value for minimum space we want to explore for this cut
    z_coords = surf.positions[:, 2]  # create list of all z-coords
    unique_surfaces = np.unique(
        z_coords.round(decimals=position_accuracy)
    )  # find all unique cuts
    mask = (
        unique_surfaces > frac_cell_value
    )  # find unique cuts for region we wish to explore
    frac_unique_surfaces = unique_surfaces[mask]
    for i in range(len(frac_unique_surfaces)):
        mask = z_coords < frac_unique_surfaces[i]
        new_surf = surf[mask]
        new_surf.center()
        rotational_symmetry = get_symmetry(new_surf, symprec=1e-5)["rotations"]
        z_rotational_operations = np.stack(rotational_symmetry)
        if np.any(z_rotational_operations[:, 2, 2] == -1):
            sym_surf = new_surf
            break
        else:
            sym_surf = False
    return sym_surf


def make_surface(thickness, miller_plane, bulk, vacuum):
    """Create a surface with 10 angstrom of vacuum provided a bulk, miller plane and layer thickness"""
    surf = surface(bulk, miller_plane, thickness)
    surf.center(vacuum=vacuum, axis=2)  # centre the surface and add vacuum
    single_layer_z_length = (surf.cell[2, 2] - (2 * vacuum)) / thickness
    return surf, single_layer_z_length


def permute_surface(bulk, miller_plane, dirname, minimum_thickness):
    """Permutes a given surface for all unique z-axis coordinates and creates an xyz file"""
    vacuum = 10
    position_accuracy = 4  # modify this depending on how accurate your atomic measurements are. A higher value generally gives more surfaces as less atoms lie on the 'same' z-plane
    surface_layers = 1
    nudge = 10 ** (-1 * position_accuracy)
    dummy_surf, single_cell_z_length = make_surface(1, miller_plane, bulk, vacuum)
    min_thick_w_buffer = minimum_thickness + (single_cell_z_length / 2)
    while ((single_cell_z_length) * (surface_layers)) < min_thick_w_buffer:
        surf, single_cell_z_length = make_surface(
            thickness=(surface_layers + 2),
            miller_plane=miller_plane,
            bulk=bulk,
            vacuum=vacuum,
        )
        surface_layers += 1
    z_coords = surf.positions[:, 2]  # create list of all z-coords
    unique_surfaces = np.unique(
        z_coords.round(decimals=position_accuracy)
    )  # create list of all unique z coordinates
    mask = unique_surfaces < single_cell_z_length + vacuum
    single_cell_unique_surfaces = unique_surfaces[mask]
    surf_lengths_list = []
    print(
        "There are "
        + str(len(single_cell_unique_surfaces))
        + " unique surface cuts for the "
        + str(miller_plane)
        + " miller plane. One layer of this surface is "
        + str(single_cell_z_length.round(position_accuracy))
        + " Angstrom thick"
    )
    for i in tqdm(range(len(single_cell_unique_surfaces))):
        surf_cut = (
            unique_surfaces[i] + nudge
        )  # define where we cut the surface on the left hand side
        ensure_thickness = int(
            i + (((surface_layers - 1) / surface_layers) * len(unique_surfaces))
        )  # integer value of number of unique values rwquired to maintain minimum thickness
        surf_term = (
            unique_surfaces[ensure_thickness] + nudge
        )  # difine where we cut the surface on the right hand side
        newmask = (z_coords > surf_cut) & (
            z_coords < surf_term
        )  # create T/F matrix for appropriate slice
        cut_surf = surf[newmask]  # apply mask to surface
        cut_surf.center(vacuum=vacuum, axis=2)  # centre and keep vacuum the same

        preceding_zero_i = "{:03d}".format(i)
        filename = (
            dirname
            + str(miller_plane[0])
            + str(miller_plane[1])
            + str(miller_plane[2])
            + "_cut_"
            + preceding_zero_i
            + "_tllzo.xyz"
        )
        # useful for debugging comment out everything else in the loop after this and uncomment here to see what the surfaces look like before being symmetrised
        # write(filename, cut_surf)
        sym_surf = find_smallest_sym_surface(
            cut_surf, surface_layers, vacuum, position_accuracy
        )  # find symmetric surfaces
        if sym_surf == False:
            continue
        surf_lengths_list.append(sym_surf.cell[2, 2] - (2 * vacuum))
        if check_surf_unique(dirname, sym_surf) > 0:
            continue
        write(filename, sym_surf)
    surf_list = glob(dirname + "/*.xyz")  # find all xyz files in a directory
    if len(surf_list) == 0:
        print("Found no symmetric structures for this Miller plane.\n")
    else:
        print(
            "Found "
            + str(len(surf_list))
            + " unique, symmetric surfaces for "
            + str(miller_plane)
            + " between "
            + str(min(surf_lengths_list).round(decimals=position_accuracy))
            + " and "
            + str(max(surf_lengths_list).round(decimals=position_accuracy))
            + " Angstrom thick ("
            + str(surface_layers)
            + " layers thick)"
            + "\n"
        )


def check_surf_unique(directory, surf):
    """Checks all surfaces in a given directory are unique"""
    surf_list = glob(directory + "/*.xyz")  # find all xyz files in a directory
    if len(surf_list) == 0:
        return 0
    atoms_dir = [read(i) for i in surf_list]
    for i in range(len(surf_list)):
        if atoms_dir[i].get_chemical_formula(
            mode="metal", empirical=False
        ) == surf.get_chemical_formula(mode="metal", empirical=False):
            new_surf_sym = get_symmetry(surf)
            old_surf_sym = get_symmetry(atoms_dir[i])
            if old_surf_sym["rotations"].all() == new_surf_sym["rotations"].all():
                return 1
            else:
                continue
    return 0


# Read in bulk cell and convert to primitve
bulk = read("./input_files/Li7La3Zr2O12.cif")

prim_spglib = find_primitive(bulk)
prim = Atoms(
    scaled_positions=prim_spglib[1], cell=prim_spglib[0], numbers=prim_spglib[2]
)

# Generate all structures up to a given miller maximum
miller_max = 2
all_miller_planes = np.mgrid[
    0 : miller_max + 1, 0 : miller_max + 1, 0 : miller_max + 1
].T.reshape(-1, 3)
allowed_miller_planes = np.delete(all_miller_planes, 0, 0)

# Alternatively provide a list of desired Miller planes to generate surfaces for
# allowed_miller_planes = [[1,0,0], [0,0,1], [2,1,0], [1,1,0]]
minimum_thickness = 8
for miller_combo in allowed_miller_planes:
    dirname = (
        "./"
        + "surfaces/"
        + str(miller_combo[0])
        + "_"
        + str(miller_combo[1])
        + "_"
        + str(miller_combo[2])
        + "/"
    )
    Path(dirname).mkdir(parents=True, exist_ok=True)
    permute_surface(prim, miller_combo, dirname, minimum_thickness)

write("./surfaces/primitive_cell.xyz", prim)
