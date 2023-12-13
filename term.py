from ase.io import read, write
from ase.build import surface
from ase.spacegroup import get_spacegroup, Spacegroup
from ase import Atoms
from spglib import get_spacegroup, get_symmetry, find_primitive
import numpy as np
from ase.visualize import view
from tqdm import tqdm
from glob import glob
from pathlib import Path

def find_smallest_sym_surface(surf):
    """ For any given surface this function will produce a symmetric surface should one exist"""
    cell = surf.cell  # find cell parameters
    halfcell_value = cell[2, 2] / 2
    z_coords = surf.positions[:, 2]  # create list of all z-coords
    unique_surfaces = np.unique(z_coords.round(decimals=4))
    mask = unique_surfaces > halfcell_value
    half_unique_surfaces = unique_surfaces[mask]
    for i in range(len(half_unique_surfaces)):
        mask = z_coords < half_unique_surfaces[i]
        new_surf = surf[mask]
        new_surf.center()
        rotational_symmetry = get_symmetry(new_surf, symprec=1e-5)["rotations"]
        z_rotational_operations = np.stack(rotational_symmetry)
        if np.any(z_rotational_operations[:, 2, 2] == -1):
            sym_surf=new_surf
            break
        else:
            sym_surf=False
    return sym_surf
        
def make_surface(thickness, miller_plane, bulk):
    """Create a surface with 10 angstrom of vacuum provided a bulk, miller plane and layer thickness"""
    surf = surface(bulk, miller_plane, thickness) 
    surf.center(vacuum=10, axis=2)  # centre the surface and add vacuum
    return surf


def permute_surface(bulk, miller_plane, dirname):
    """ Permutes a given surface for all unique z-axis coordinates and creates an xyz file """
    surf = make_surface(3, miller_plane, bulk)
    cell = surf.cell  # find cell parameters
    single_cell_distance = (cell[2, 2]-20) / 3
    z_coords = surf.positions[:, 2]  # create list of all z-coords
    position_accuracy=4 # modify this depending on how accurate your atomic measurements are. A higher value generally gives more surfaces as less atoms lie on the 'same' z-plane
    unique_surfaces = np.unique(z_coords.round(decimals=position_accuracy))
    mask = unique_surfaces < single_cell_distance +10
    single_cell_unique_surfaces = unique_surfaces[mask]
    print(print("There are " + str(len(single_cell_unique_surfaces)) +
          " unique surfaces for the " + str(miller_plane) + " miller plane"))
    for i in tqdm(range(len(single_cell_unique_surfaces))):
        surf_cut = unique_surfaces[i]+0.0001
        ensure_2_thick = int(i + ((2 / 3) * len(unique_surfaces)))
        surf_term = unique_surfaces[ensure_2_thick]+0.0001
        newmask = (z_coords > surf_cut) & (z_coords < surf_term)
        cut_surf = surf[newmask]
        cut_surf.center(vacuum=10, axis=2)
        sym_surf=find_smallest_sym_surface(cut_surf)
        if sym_surf==False:
            continue
        preceding_zero_i = "{:03d}".format(i)
        filename = dirname + str(miller_plane[0]) + str(miller_plane[1]) + str(miller_plane[2]) + "_cut_" + preceding_zero_i + "_tllzo.xyz"
        write(filename, sym_surf)

def check_surf_unique(directory):
    """ Checks all surfaces in a given directory are unique """
    surf_list=glob(directory + '/*.xyz') # find all
    atoms_dir = [read(i) for i in surf_list]
    atoms_list= [read(i) for i in surf_list] #change here
    same_count=0
    for i in range(len(surf_list)):
        for j in range(len(surf_list)):
            if i >= j:
                continue
            if len(atoms_list[i]) != len(atoms_dir[j]):
                continue
            if np.allclose(atoms_list[i].positions, atoms_dir[j].positions):
                print(str(surf_list[i]) + ' same as ' + str(surf_list[j]))
                same_count+=1
    if same_count==0:
        print("All surfaces are unique")

# read in bulk cell and convert to primitve
bulk = read('Li7La3Zr2O12.cif')
#bulk = read('Fe2O3_3.cif')
prim_spglib=find_primitive(bulk)
prim=Atoms(scaled_positions=prim_spglib[1], cell=prim_spglib[0], numbers=prim_spglib[2])

#generate all structures up to a given miller maximum
""" miller_max=2
all_miller_planes = np.mgrid[0:miller_max+1, 0:miller_max+1, 0:miller_max+1].T.reshape(-1, 3)
allowed_miller_planes=np.delete(all_miller_planes, 0, 0)
 """
allowed_miller_planes = [[1,0,0]]
for miller_combo in allowed_miller_planes:
    dirname='./' + 'surfaces/' + str(miller_combo[0]) + '_' + str(miller_combo[1]) + '_' + str(miller_combo[2]) +'/'
    Path(dirname).mkdir(parents=True, exist_ok=True)
    print("Finding surfaces for the Miller plane: " + str(miller_combo))
    permute_surface(prim, miller_combo, dirname)

write('./surfaces/primitive_cell.xyz', prim)
check_surf_unique('.')

