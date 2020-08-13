Topology
========

Importing structure
-------------------

To create bond topology you first need a LAMMPS data file with the `Masses` and `Atoms`
sections. This can be generated from any structure file (cif, xyz, pdb, vesta, etc.) by
the command-line tool `atomsk`::

    atomsk structure.cif lmp

It can also be generated other ways such as `topotools` or `pymatgen`::

    import pymatgen.io.lammps.data as pgld
    from pymatgen import Structure

    pg_structure = Structure.from_file('structure.cif')
    ld_object = pgld.LammpsData.from_structure(pg_structure)
    ld_object.write_file('structure.lmp')

The `Masses` section must have a comment after each mass with the name
of the element. Element names are added automatically by `atomsk` but not by `pymatgen`. For example::

    Masses

    1   207.200000000 # Pb
    2    14.006700000 # N
    3    12.010700000 # C
    4   126.904470000 # I

The `Atoms` section lists all atom coordinates. It can optionally have a column with the
atomic charges.

Creating bond topology
----------------------
First specify pairs of atoms to bond and the maximum distance for bonding. The following
will create bonds between nitrogen and carbon atoms up to 1.5 Angstroms and
between carbon atoms up to 1.52 Angstroms::

    import pandas as pd

    bond_pairs = pd.DataFrame(data = dict(element1 = ['N',  'C'],
                                          element2 = ['C',  'C'],
                                          cutoff   = [1.5, 1.52]))


To add bond, angle and dihedral information to the LAMMPS file::

    import lammps_helper as lh
    bonds =  lh.add_bond_data('lammps_data.lmp', bond_pairs)
    angles, dihedrals = lh.add_angle_dihedral_data('lammps_data.lmp', bonds)

Keep in mind that the topology must be very simple. This code probably cannot handle
cyclic molecules. If your structure has molecules with rings, I suggest you use
moltemplate to generate the structure or lammps-interface to convert the structure from
a .cif file.