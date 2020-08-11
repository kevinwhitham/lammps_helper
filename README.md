# lammps_helper
## Install

`pip install lammps_helper`

## About

**lammps_helper** is python code to help create input files for and extract output data from the molecular dynamics package LAMMPS. It contains just a few features that I needed for my specific project and is not meant to be general in any way.

You're probably asking yourself 'Why do I need to download another package? Why doesn't this already exist?' So here is a table comparing lammps_helper to some other packages.


|Package          | Features | Issues |
|-----------------|----------|---------|
| **lammps_helper** | Adds topology to LAMMPS data | Limited to simple structures
|                   | Extracts info about dipoles, thermo data. |
|**lammps-interface** | Converts cif to LAMMPS data     | Doesn't apply symmetry operations when building supercells. Works for cif files with P1 symmetry.|
|                     | Sophisticated topology detection |  |
| **atomsk**          | Converts all kinds of structure files including LAMMPS | Doesn't create topology (bonds, angles, dihedrals) |
|                     | Manipulates structures (create, delete, and move atoms)|
| **topotools**       | All kinds of manipulations on atoms, molecules, and topology. | Can either set bonds for all atoms or single bonds. |
|                     | Great for organics, polymers, proteins. | A bit awkward to use with crystals.  |
| **Vesta**           | Converts and displays all kinds of files. | Doesn't output topology. |
| **atomman**         | Easy interaction with LAMMPS from python. Setup and run simulations, extract output. Imports cif. | Doesn't deal with topology.
| **lammps-cython** | Another python wrapper for LAMMPS. Setup and run, extract output. | No structural manipulation.
| **pymatgen** | General suite of tools for handling materials computations via python. Includes some LAMMPS input/output and control. Has a function that generates bonds. | No way to set specific bonds, uses a global distance parameter.
| **moltemplate** | Generates LAMMPS data using molecule building instructions. | Can't import .cif, .xyz, etc. 



Obviously there is more than one workflow to create structure data for LAMMPS. The simplest is to write the LAMMPS data file by hand. That works for very simple structures up to about 10 atoms and/or bonds per unit cell. If there are no explicit bonds in your structure you can simply convert a .cif, .xyz or similar file into LAMMPS format with (from easiest to hardest) atomsk, topotools, lammps-interface, atomman, or pymatgen. If your structure has bonds and you don't want to write out all the bonds, angles, and dihedrals by hand you have a couple options. 

1. Write out instructions to build the structure with moltemplate. Requires a lot of time up front, but makes it easy to change the structure in future.

2. Convert a .cif, .xyz, etc. file to a LAMMPS data file. Fast and easy, but tweaking the structure is tedious. 

If you choose to convert a .cif, .xyz, etc. it can be done several ways. One way is to use atomsk to convert the cif file to a LAMMPS file like this:

`atomsk crystal.cif lmp`

That will give you a LAMMPS data file with a list of all the atoms in the unit cell. You can then add bonds, angles, and dihedrals using lammps_helper like this:

```
import lammps_helper as lh
import pandas as pd

bond_pairs = pd.DataFrame(data = dict(element1 = ['N',  'C'],
                                       element2 = ['C',  'C'],
                                       cutoff   = [1.5, 1.52]))

bonds =  lh.add_bond_data(f'my_structure.lmp', bond_pairs)
angles = lh.add_angle_dihedral_data(f'my_structure.lmp', bonds)
```

With lammps_helper you must specify each type of bond by the names of the atoms and a cutoff distance (in whatever units are used in the LAMMPS file, probably Angstroms). Any atoms of those types within that distance will be bonded. This code is very simple and probably can't handle cyclic molecules. If you have a structure with rings I suggest you try moltemplate or lammps-interface.

lammps_helper also has code useful for extracting data about molecular dipoles and for plotting distributions of dipole directions over time. See [the documentation](https://lammps-helper.readthedocs.io/en/latest/) for more info. lammps_helper is available on [github.](https://github.com/kevinwhitham/lammps_helper)

