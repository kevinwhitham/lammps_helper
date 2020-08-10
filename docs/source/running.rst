Running LAMMPS with lammps_helper
=================================

lammps_helper provides two things to help run LAMMPS simulations from python. One to help
create LAMMPS input files and one to run them. I suggest you use an IPython environment
such as Jupyter Lab.

Creating LAMMPS input files:
----------------------------
At the beginning of your Jupyter Lab notebook, use this code::

    import lammps_helper as lh

    ip = get_ipython()
    ip.register_magics(lh.lammps_magics)

Then you can write LAMMPS files in your notebook and export them rather than writing
them separately. Using the `%%writetemplate` cell magic, any python variables
in your LAMMPS code will be replaced by their values during export. For example::

    temperature = 300
    output_base_name = f'my_lammps_output_{temperature}K'
    simulation_steps = 20000
    list_element_names = 'Pb I N C H'

::

    %%writetemplate in.my_lammps_input_file

    # This is a LAMMPS input file

    units      real
    dimension  3
    boundary   p p p
    atom_style full

    pair_style      lj/cut 12.500
    bond_style harmonic
    angle_style harmonic
    dihedral_style charmm

    read_data data.my_lammps_data_file

    thermo_style custom step temp etotal spcpu cpuremain

    # equillibrate at temperature
    fix fix_npt_equilibrate all npt temp {temperature} {temperature} 100 iso 1.0 1.0 500
    dump dump_trajectory_equilibrate all dcd 100 {output_base_name}_equilibration.dcd
    dump_modify dump_trajectory_equilibrate unwrap yes
    run {simulation_steps}

    # dump the system to check geometry
    write_dump all xyz {output_base_name}_after_equilibrate.xyz modify element {list_element_names}

    unfix fix_npt_equilibrate
    undump dump_trajectory_equilibrate

    # write restart file after equilibration
    write_restart {output_base_name}_after_equilibrate.restart

Will create a file named `in.my_lammps_input_file` that will run a simulation at 300 K
for 20000 timesteps and the data files created by LAMMPS will share a common naming
scheme that includes the temperature.

Running LAMMPS
--------------

To run your input file in LAMMPS, do this ::

    lh.run_lammps('in.my_lammps_input_file', f'log_{output_base_name}_equilibrate.lammps')

This simply issues the lmp_serial command. You may add variables here as well.
See :func:`run_lammps` for more info. You will get output like this ::

    Running calculation...
    Writing to log_my_lammps_output_300K_equilibrate.lammps.
    Calculation started at 08-06-20 16:10:45
    Calculation complete at 08-06-20 16:10:45 in 0:00:00.076904.

