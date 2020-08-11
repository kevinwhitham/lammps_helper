try:
    shell = get_ipython().__class__.__name__
    if shell is not None:
        from IPython.core.magic import (Magics, magics_class, line_magic, cell_magic, line_cell_magic,
                                        needs_local_scope)


        @magics_class
        class lammps_magics(Magics):
            """
            use %%writetemplate instead of %%writefile to substitute Python variables
            e.g. 'dump 1 all xyz {output_base_name}.xyz'
            output_base_name will be replaced with the value of the variable
            """

            @line_cell_magic
            @needs_local_scope
            def writetemplate(self, line, cell=None, local_ns=None):
                with open(line, 'w') as f:
                    args = eval('globals()', self.shell.user_ns, local_ns)
                    f.write(cell.format(**args))
                    print(f'Overwriting {line}')

except NameError:
    print('IPython not detected. %%writetemplate not available.')
        

import subprocess 
import shlex
import time
import datetime
import numpy as np
import pandas as pd
import mmap
from io import StringIO
import re
import plotly.graph_objects as go

#------------------------------------ Running LAMMPS -------------------------------------------------
def run_lammps(input_file, log_file, variables = None):
    """
    Run lmp_serial -in input_file variables -log log_file
    Provides start, end, elapsed time.

    Parameters
    ----------
    input_file : str
        path to LAMMPS input file
    log_file : str
        path to LAMMPS log file
    variables : list of str
        variables to feed LAMMPS on the command line, e.g. ['x = 1', 'y = 2']

    Returns
    -------

    """
    print('Running calculation...')
    print(f'Writing to {log_file}.')
    
    if variables:
        vars_string = ''.join(('-var {} {} '.format(*i) for i in variables.items()))
    else:
        vars_string = ''
    
    start_time = time.perf_counter()
    
    print(f'Calculation started at {time.strftime("%m-%d-%y %H:%M:%S")}')
    
    command_line = f'lmp_serial -in {input_file} {vars_string} -log {log_file}'
    args = shlex.split(command_line)
    p = subprocess.Popen(args)
    result = p.wait()
    
    end_time = time.perf_counter()

    print(f'Calculation complete at {time.strftime("%m-%d-%y %H:%M:%S")} in {str(datetime.timedelta(seconds=(end_time - start_time)))}.')

#------------------------------------ Log File -------------------------------------------------
def get_thermo_data_from_log(filename):
    """
    Converts the data output by the LAMMPS thermo_style command to an array.

    Parameters
    ----------
    filename : str
        path to LAMMPS log file

    Returns
    -------
    numpy array
        thermo data with columns matching the LAMMPS thermo_style setting

    """
    with open(filename, 'r+') as f:
        file_txt = mmap.mmap(f.fileno(), 0)
        result = re.search(b'(Step Temp.*?\n)(.*?)(?=Loop?)', file_txt, re.DOTALL)
        names = re.findall('\S+', result.group(1).decode('ascii'))
        log_txt = StringIO(result.group(2).decode('ascii'))

        print(f'Names of variables in log data: {names}')

        # debug
        # print(log_txt.readline())
        # print(log_txt.readline())

        # don't raise an exception if there are other lines of output in the thermo data
        return np.genfromtxt(log_txt, names=names, invalid_raise=False)

#------------------------------------ Dipoles -------------------------------------------------
def get_dipole_data(dipole_data, dipole_file_name, location_file_name, temperature, append = False):
    """
    Convert output from LAMMPS compute_dipole to an array.

    Parameters
    ----------
    dipole_data : numpy array
        An array to hold dipole data.
        If the array is pre-allocated with enough rows and cols, will attempt to fill the array. (maybe faster)
        Otherwise will create a new array and grow it as data is imported. (maybe slower)
    dipole_file_name : str
        path to LAMMPS output file containing compute_dipole data
    location_file_name : str
        path to LAMMPS output file containing compute_com data
    temperature : str, int
        temperature of the simulation for book-keeping
    append : bool
        If True, append imported data to the dipole_data array.
        Useful for adding data from multiple runs, i.e. at different temperatures.

    Returns
    -------
    DataFrame
        Pandas DataFrame with columns
            temperature : temperature during simulation
            timestep : timestep when the dipole was sampled
            molecule : index of the molecule with the dipole
            Dx : x-component of the dipole moment vector
            Dy : y-component of the dipole moment vector
            Dz : z-component of the dipole moment vector
            moment : total moment
            com_x : x-coordinate of center of mass
            com_y : y-coordinate of center of mass
            com_z : z-coordinate of center of mass
            cos_theta : angle of the dipole relative to the +z-axis
            phi : angle of the dipole relative to the +x-axis

    int
        number of dipole measurements (rows in the DataFrame)

    """
    
    # Assume that if the shape is wrong, we can discard any data in the array
    if dipole_data.shape[1] != 12:
        dipole_data = np.empty((dipole_data.shape[0],12))
    
    with open(dipole_file_name, 'r+') as dipole_file:
        with open(location_file_name, 'r+') as location_file:
            
            # StringIO and mmap execute in about the same time
            dipole_txt   = mmap.mmap(dipole_file.fileno(), 0)    #StringIO(dipole_file.read())
            location_txt = mmap.mmap(location_file.fileno(), 0)  #StringIO(location_file.read())
            
            # skip header
            for i in range(3):
                dipole_txt.readline()
                location_txt.readline()

            done = False
            
            if append:
                row = dipole_data.shape[0]
            else:
                row = 0

            while not done:

                try:
                    # get timestep and number of dipoles
                    #dipole_chunk_data = dipole_chunk_data_iterator.get_chunk(1)
                    dipole_chunk_data = dipole_txt.readline()
                    location_chunk_data = location_txt.readline()
                    
                    dipole_items   = re.findall(b'(\S+)', dipole_chunk_data)
                    location_items = re.findall(b'(\S+)', location_chunk_data)
                    
                    if len(dipole_items) and len(location_items):
                        
                        # This assumes the data rows contain more than two columns
                        if (len(dipole_items) == 2) and (len(location_items) == 2):
                            
                            dipole_timestep, num_dipoles = dipole_items
                            location_timestep, num_locations = location_items

                            assert location_timestep == dipole_timestep
                            assert num_locations == num_dipoles

                        else:
                            
                            # make we're reading the same molecule from each file
                            assert dipole_items[0] == location_items[0]
                            
                            # if the dipole moment is 0, there is no dipole, discard it
                            if np.float(dipole_items[4]) != 0.0:
                                
                                Dx     = np.float(dipole_items[1])
                                Dy     = np.float(dipole_items[2])
                                Dz     = np.float(dipole_items[3])
                                moment = np.float(dipole_items[4])
                                
                                cos_theta = Dz / moment
                                phi   = np.arctan2(Dy, Dx) * 180.0 / np.pi
                                    
                                data_row = np.array([temperature,
                                                     int(dipole_timestep),
                                                     int(dipole_items[0]),
                                                     Dx,
                                                     Dy,
                                                     Dz,
                                                     moment,
                                                     float(location_items[1]),
                                                     float(location_items[2]),
                                                     float(location_items[3]),
                                                     cos_theta,
                                                     phi])
                                    
                                if row < dipole_data.shape[0]:
                                    
                                    dipole_data[row] = data_row
                                    
                                else:

                                    dipole_data = np.vstack((dipole_data, data_row))
                                    
                                row = row + 1
                                
                                # debug
                                #if (dipole_data.shape[0] % 1000) == 0:
                                #    print(f'Data rows: {dipole_data.shape[0]}')
                            
                    else:
                        dipole_txt.close()
                        location_txt.close()
                        
                        # Convert dipole_data to a pandas DataFrame
                        dipole_data = pd.DataFrame(dipole_data, columns=['temperature',
                                                                         'timestep',
                                                                         'molecule',
                                                                         'Dx',
                                                                         'Dy',
                                                                         'Dz',
                                                                         'moment',
                                                                         'com_x',
                                                                         'com_y',
                                                                         'com_z',
                                                                         'cos_theta',
                                                                         'phi'])
                        
                        print(f'Finished reading {dipole_file_name} and {location_file_name}.')
                        done = True


                except Exception as e:
                    print(f'Exception: {e}')
                    dipole_txt.close()
                    location_txt.close()
                    print(f'Error reading {dipole_file_name} and {location_file_name}.')
                    done = True
            
    return (dipole_data, row)

def make_dipole_contour_plot(dipole_data, title = None, subtitle = None):
    """
    Creates a 2D histogram of dipole orientations.
    Axes are cos_theta and phi, the orientation in polar coordinates.
    Use .show() on the returned figure to view.

    Parameters
    ----------
    dipole_data : array-like
        must contain columns 'cos_theta' and 'phi'

    title : str
        plot title (optional)

    subtitle : str
        plot subtitle (optional)

    Returns
    -------
    plotly Figure

    """
    
    fig = go.Figure()
    
    fig.add_trace(go.Histogram2dContour(
        x = dipole_data['phi'],
        y = dipole_data['cos_theta'],
        colorscale = 'Blues',
        showscale = False,
        reversescale = True,
        xaxis = 'x',
        yaxis = 'y',
        hovertemplate = 'cos(theta): %{y:.2f}<br>'+
                            'phi: %{x}<br>'+
                            'count: %{z}<extra></extra>'
    ))

    fig.add_trace(go.Histogram(
            y = dipole_data['cos_theta'],
            xaxis = 'x2',
            marker = dict(
                color = 'rgba(0,0,0,1)'
            )
        ))
    fig.add_trace(go.Histogram(
            x = dipole_data['phi'],
            yaxis = 'y2',
            marker = dict(
                color = 'rgba(0,0,0,1)'
            ),
        ))

    fig.update_layout(
        autosize = False,
        title_text = title,
        annotations=[
                     dict(x = 0,    y = 1.1, showarrow = False, xref = 'paper', yref = 'paper', text = subtitle),   # subtitle
                     dict(x = -180, y = 1.05,   showarrow = False, xref = 'x',     yref = 'paper',    text = '-x'),
                     dict(x =  -90, y = 1.05,   showarrow = False, xref = 'x',     yref = 'paper',    text = '-y'),
                     dict(x =    0, y = 1.05,   showarrow = False, xref = 'x',     yref = 'paper',    text = '+x'),
                     dict(x =   90, y = 1.05,   showarrow = False, xref = 'x',     yref = 'paper',    text = '+y'),
                     dict(x =  180, y = 1.05,   showarrow = False, xref = 'x',     yref = 'paper',    text = '-x'),
                     dict(x =  1.1, y = 1,      showarrow = False, xref = 'paper', yref = 'y',        text = '+z'),
                     dict(x =  1.1, y = 0,      showarrow = False, xref = 'paper', yref = 'y',        text = 'xy'),
                     dict(x =  1.1, y = -1,     showarrow = False, xref = 'paper', yref = 'y',        text = '-z')
        ],
        xaxis = dict(
            title = 'phi (deg.)',
            zeroline = False,
            domain = [0,0.85],
            showgrid = False
        ),
        yaxis = dict(
            title = 'cos(theta)',
            zeroline = False,
            domain = [0,0.85],
            showgrid = False
        ),
        xaxis2 = dict(
            zeroline = False,
            domain = [0.85,1],
            showgrid = False,
            showticklabels = False
        ),
        yaxis2 = dict(
            zeroline = False,
            domain = [0.85,1],
            showgrid = False,
            showticklabels = False
        ),
        height = 600,
        width = 600,
        bargap = 0,
        hovermode = 'closest',
        showlegend = False
    )
    
    return fig

def plot_mean_dipole_orientation(dipole_data):
    """
    Creates a 3D plot with each molecule represented by a cone.
    The cone points in the mean direction of the dipole.
    The size of the cone is proportional to the total dipole moment.
    Call .show() on the returned figure to view.

    Parameters
    ----------
    dipole_data : pandas DataFrame
        generated by get_dipole_data
        must contain the columns molecule, Dx, Dy, Dz, com_x, com_y, com_z

    Returns
    -------
    plotly figure

    """
    mean_uvw = dipole_data.pivot_table(index='molecule', values=['Dx','Dy','Dz'], aggfunc='mean')
    mean_xyz = dipole_data.pivot_table(index='molecule', values=['com_x', 'com_y', 'com_z'], aggfunc='mean')

    fig = go.Figure(go.Cone(x = mean_xyz['com_x'],
                            y = mean_xyz['com_y'],
                            z = mean_xyz['com_z'],
                            u = mean_uvw['Dx'],
                            v = mean_uvw['Dy'],
                            w = mean_uvw['Dz']))
    
    return fig

def plot_mean_dipole_angles(dipole_data, title = 'Mean Dipole Angles'):
    """
    Create 3D plots of dipole angles cos(theta) and phi. The angles are represented using
    a color scale in a volume plot.

    Parameters
    ----------
    dipole_data : pandas DataFrame
        generated by get_dipole_data
        must contain columns molecule, com_x, com_y, com_z, cos_theta, phi

    title : str

    Returns
    -------
    plotly figures
        cos(theta)
        phi
    """

    fig_cos = go.Figure()
    fig_phi = go.Figure()

    means = dipole_data.pivot_table(index='molecule', values=['com_x', 'com_y', 'com_z', 'cos_theta', 'phi'],
                                           aggfunc='mean')

    xmin = dipole_data['com_x'].min()
    xmax = dipole_data['com_x'].max()
    ymin = dipole_data['com_y'].min()
    ymax = dipole_data['com_y'].max()
    zmin = dipole_data['com_z'].min()
    zmax = dipole_data['com_z'].max()


    grid_x, grid_y, grid_z = np.mgrid[xmin:xmax:((xmax-xmin)/5j), ymin:ymax:((ymax-ymin)/5j), zmin:zmax:((zmax-zmin)/5j)]

    fig_cos.add_trace(go.Volume(x=grid_x.flatten(),
                            y=grid_y.flatten(),
                            z=grid_z.flatten(),
                            value=means['cos_theta'].to_numpy(),
                            opacity=0.2,
                            surface_count=20,
                            coloraxis='coloraxis'))

    fig_phi.add_trace(go.Volume(x=grid_x.flatten(),
                            y=grid_y.flatten(),
                            z=grid_z.flatten(),
                            value=means['phi'].to_numpy(),
                            opacity=0.2,
                            surface_count=20,
                            coloraxis='coloraxis'))

    fig_cos.update_layout(title=title,
                      height=600,
                      width=1200,
                      coloraxis=dict(colorscale='RdBu_r', colorbar_title_text='cos(theta)'),
                      font_size=16)

    fig_phi.update_layout(title=title,
                          height=600,
                          width=1200,
                          coloraxis=dict(colorscale='RdBu_r', colorbar_title_text='phi'),
                          font_size=16)

    return fig_cos, fig_phi


# -----------------------------------  Topology ------------------------------
def add_bond_data(lammps_xyz_file, bond_pairs):
    """
    Add bond information to a LAMMPS structure data file.

    If you have a file that gives atomic coordinates (a .cif or .xyz file),
    use atomsk to convert to a .lmp (LAMMPS) data file.

    This function uses atomic coordinates in a LAMMPS file to find and add bonds.
    Bonds are made according to the bond_pairs parameter.

    Parameters
    ----------
    lammps_xyz_file : str
        path to a LAMMPS data file with atom types and atomic coordinates
    bond_pairs : pandas DataFrame
        Each row contains
            element1 : str
                element name (e.g. 'N')
            element2 : str
                element name (e.g. 'C')
            cutoff : float
                maximum distance for a bond between atoms of element1 and element2

    Returns
    -------
    numpy array
        column:  0                1           2           3       4          5          6         7  8  9 10
        bond_index bond_type_number atom1_index atom2_index comment atom1_type atom2_type Boundary: nx ny nz

        bond_index : int
            index of the bond
        bond_type_number : int
            bond type
        atom1_index : int
            index of one of the atoms in the bond
        atom2_index : int
            index of the other atom in the bond
        comment : str
            a hash (#) character followed by the names of the atoms in the bond, e.g. # N - C
            all columns after this one are part of the comment
        atom1_type : int
            type of the atom with index atom1_index
        atom2_type : int
            type of the atom with index atom2_index
        Boundary: : str
            the string 'Boundary:' signifies that a bond crosses the cell/box boundary
        nx : int
            if a bond crosses the cell/box boundary, part of the molecule will be translated by nx unit
            vectors in the x-direction
        ny : int
            if a bond crosses the cell/box boundary, part of the molecule will be translated by ny unit
            vectors in the y-direction
        nz : int
            if a bond crosses the cell/box boundary, part of the molecule will be translated by nz unit
            vectors in the z-direction


    """
    
    with open(lammps_xyz_file, 'r') as data_file:
        file_text = data_file.read()

    element_names = list()
    atom_types = list()

    # bonds structure
    # column:  0                1           2           3       4          5          6         7  8  9 10
    # bond_index bond_type_number atom1_index atom2_index comment atom1_type atom2_type Boundary: nx ny nz
    bonds = np.empty((0, 11))

    # bonds could be a structured array for better code reuse
    #dtype=[('bond_index', np.int), 
    #                                 ('bond_type',  np.int), 
    #                                 ('atom1_index',np.int),
    #                                 ('atom2_index',np.int),
    #                                 ('comment',    np.str),
    #                                 ('atom1_type', np.int),
    #                                 ('atom2_type', np.int),
    #                                 ('boundary_label',np.str),
    #                                 ('nx',np.int),
    #                                 ('ny',np.int),
    #                                 ('nz',np.int)]
    
    # get unit cell vectors, assuming the box size is the unit cell size
    #  0.00000000      27.68200000  xlo xhi
    #  0.00000000       8.87670000  ylo yhi
    #  0.00000000       8.69620000  zlo zhi
    
    directions = ['x', 'y', 'z']
    unit_vectors = np.array([0.0,0.0,0.0], dtype=np.float32)
    
    for direction in directions:
        search_result = re.search(f'([\d|\.]+)\s+([\d|\.]+)\s+{direction}lo\s+{direction}hi', file_text)
        
        if search_result:
            unit_vectors[directions.index(direction)] = abs(float(search_result.group(2)) - float(search_result.group(1)))
            
            
    # debug
    print(f'Unit vectors: {unit_vectors}')
    
    # Get the coordinates of all atoms, add three columns for translation flags
    # atom_coords structure
    # column:  0         1 2 3 4      5      6      7
    # atom_index atom_type x y z x_flag y_flag z_flag
    atom_coords, charges = get_atom_coordinates(file_text)
    atom_coords = np.hstack((atom_coords, np.zeros((atom_coords.shape[0],3))))

    for bond_type_index, pair in bond_pairs.iterrows():

        print(f'Next pair: {pair["element1"]} - {pair["element2"]} cutoff: {pair["cutoff"]}')

        # Get atom type number that corresponds to the element name

        for element in [pair['element1'], pair['element2']]:

            if element not in element_names:

                search_result = re.search(f'(\d+).*?#\s*{element}', file_text)

                if search_result:

                    element_names.append(element)
                    atom_types.append(search_result.group(1))

                else:
                    print(f'No atom type found for atom named {element}. Add comment with name in Masses section.')
                    assert 0


        assert len(element_names) == len(atom_types)

        # Get the two types of atoms to look for
        type1 = int(atom_types[element_names.index(pair['element1'])])
        type2 = int(atom_types[element_names.index(pair['element2'])])
        
        type1_coords = atom_coords[atom_coords[:,1] == type1]

        # Get the coordinates of all type2 atoms
        if type1 == type2:
            type2_coords = type1_coords
        else:
            type2_coords = atom_coords[atom_coords[:,1] == type2]

        # debug
        #print(f'Found {type1_coords.shape[0]} {pair["element1"]} (type {type1}) atoms at (index, type, x,y,z):')
        #print(type1_coords)
        #print(f'Found {type2_coords.shape[0]} {pair["element2"]} (type {type2}) atoms at (index, type, x,y,z):')
        #print(type2_coords)

        # For each atom of type1, check distance to every type2 atom
        # if the distance is less than the cutoff, add a bond
        for atom1, atom1_row_num in zip(type1_coords, range(type1_coords.shape[0])):

            if type1 == type2:

                # when the atoms are the same type, don't double count
                type2_search_array = type2_coords[(atom1_row_num+1):]

            else:

                type2_search_array = type2_coords

            for atom2 in type2_search_array:

                # atom1, atom2 data structure:
                # column:  0          1 2  3  4      5      6      7
                # atom_index, atom_type x, y, z x_flag y_flag z_flag
                
                # construct a matrix to translate atom coordinates by unit cell vectors
                cell_transforms = np.array([[0,0,0],
                                            [1,0,0],
                                            [0,1,0],
                                            [0,0,1],
                                            [1,1,0],
                                            [0,1,1],
                                            [1,0,1],
                                            [1,1,1]], dtype=np.int)
                
                # transform_directions looks the same as cell_transforms, but with 0 -> 1 and 1 -> -1
                transform_directions = cell_transforms + 1
                np.subtract(transform_directions, 3 * cell_transforms, out = transform_directions, where = (cell_transforms == 1))

                for transform in cell_transforms:
                    
                    for transform_direction in transform_directions:
                        
                        atom2_x = (atom2[2] - transform[0] * unit_vectors[0] * transform_direction[0])
                        atom2_y = (atom2[3] - transform[1] * unit_vectors[1] * transform_direction[1])
                        atom2_z = (atom2[4] - transform[2] * unit_vectors[2] * transform_direction[2])

                        distance = np.sqrt( (atom1[2] - atom2_x)**2 +
                                            (atom1[3] - atom2_y)**2 +
                                            (atom1[4] - atom2_z)**2   )
                        
                        # debug
                        #print(f'Atom {atom1[0]} and {atom2[0]} distance: {distance}, transform: {transform * unit_vectors * transform_direction}')

                        if distance < pair['cutoff']:
                            
                            # Need to set the flags for all atoms in the molecule on one side of the boundary
                            # not just one of the two atoms that form a bond across the boundary
                            #atom_coords[atom_coords[:,0] == atom2[0],5:8] = transform_flags
                            #print(f'Setting flags for atom {atom2[0]}: {transform_flags}')
                            # store the transforms made to atoms across the box boundary
                            transform_flags = np.array(-transform * transform_direction, dtype='<i8')

                            # save the bond index, bond type number, indices of the bonded atoms, and a comment with the atom names and types
                            
                            bonds = np.vstack((bonds, [bonds.shape[0]+1,
                                                       bond_type_index+1,
                                                       int(atom1[0]),
                                                       int(atom2[0]),
                                                       f'# {pair["element1"]} - {pair["element2"]}',
                                                       type1,
                                                       type2,
                                                       'Boundary:',
                                                       int(transform_flags[0]),
                                                       int(transform_flags[1]),
                                                       int(transform_flags[2])]))

                            break

        print(f'Bonds found: {bonds.shape[0]}')
       
    # debug
    #print(bonds)
    #print('Translation flags')
    #print(atom_coords[:,-3:])
    
    # Set translation flags for bonds that cross the boundary
    for bond, bond_index in zip(bonds, range(bonds.shape[0])):
        
        # check if this bond goes through the boundary
        if np.any(bond[8:11].astype(np.int) != 0):
            
            up_atoms = [int(bond[2])]
            down_atoms = [int(bond[3])]
            
            # make two lists of atoms on either side of the boundary bond
            find_connected_atoms(int(bond[2]), np.delete(bonds, bond_index, axis=0), up_atoms)
            find_connected_atoms(int(bond[3]), np.delete(bonds, bond_index, axis=0), down_atoms)
            
            #print(f'up_atoms: {up_atoms}')
            #print(f'down_atoms: {down_atoms}')
            
            # flag the smaller part of the molecule to be translated
            if len(up_atoms) < len(down_atoms):
                
                for atom_index in up_atoms:
                    
                    # translation flags are set for the second atom, so negate the flags here
                    atom_coords[atom_coords[:,0] == atom_index, 5:8] = -1 * bond[8:11].astype(np.int)
                    
            else:
                
                for atom_index in down_atoms:
                    
                    atom_coords[atom_coords[:,0] == atom_index, 5:8] = bond[8:11].astype(np.int)

    # Replace atomic coordinates section of the data file to include the translation flags nx, ny, nz
    # This is preferable to translating the coordinates outside the box because 
    # LAMMPS will move all atoms inside the box on import
    
    string_stream = StringIO()
    
    # Example LAMMPS data file line of atomic coordinates
    # index type  charge       x                  y               z            nx ny nz
    # 1     1     0.00         13.84100000       4.43835000       4.34810000   0  0  0
    
    if charges.shape[0]:
        
        np.savetxt(string_stream, np.insert(atom_coords, 3, charges, axis=1), fmt=['%i', '%i', '%f', '%f', '%f', '%f', '%i', '%i', '%i'])
        
    else:
        
        np.savetxt(string_stream, atom_coords,  fmt=['%i', '%i', '%f', '%f', '%f', '%i', '%i', '%i'])
    
    file_text_with_flags = re.sub('(?s)(\nAtoms.*?\n)(.*?)((\n\w.*)|$)', f'\g<1>\n{string_stream.getvalue()}\n\n\g<3>', file_text)
    string_stream.close()
    
    if len(file_text_with_flags) != len(file_text):
        
        # rewrite the lammps file
        with open(lammps_xyz_file, 'w') as data_file:
            data_file.write(file_text_with_flags)
            data_file.write('\nBonds\n\n')
            np.savetxt(data_file, bonds, fmt='%s')
    else:
        # something went wrong with the substitution of atomic coordinates
        print('Could not replace atomic coordinates section, just adding bonds section.')
        
        # write bonds section to the lammps file
        with open(lammps_xyz_file, 'a') as data_file:
            data_file.write('\nBonds\n\n')
            np.savetxt(data_file, bonds, fmt='%s')


    # Add number of bonds to file header
    with open(lammps_xyz_file, 'r') as data_file:
        file_text = data_file.read()

    file_text_modified = re.sub('\n(\s*)(\d+)(\s+)(atoms.*?)',
                                f'\n\g<1>\g<2>\g<3>\g<4>\g<1>{bonds.shape[0]}\g<3>bonds',
                                file_text)

    if len(file_text_modified) == len(file_text):
        print('Could not add number of bonds to file header')
    else:
        with open(lammps_xyz_file, 'w') as data_file:
            data_file.write(file_text_modified)
    
    print(f'Added {bonds.shape[0]} bonds to {lammps_xyz_file}')
    
    return bonds
        
def get_atom_coordinates(text):
    """

    Extracts atom coordinates from a LAMMPS data file.
    Returns an array with atom coordinates.

    If the Atoms section contains charge data, returns a 1-D array with charges
    in the same order. Otherwise the charge data array is empty.

    Parameters
    ----------
    text : str
        contents of a LAMMPS data file with atomic coordinates

    Returns
    -------
    numpy array
        column:  0         1 2 3 4
        atom_index atom_type x y z

    numpy array
        column:   0
        atom_charge


    """
    
    # coords structure
    # column:  0         1 2 3 4
    # atom_index atom_type x y z
    coords = np.empty((0,5), dtype=np.float32)
                                   
    charge_data = np.empty((0,1), dtype=np.float32)
        
    # Example LAMMPS data file line of atomic coordinates
    # index type  charge       x                  y               z
    # 1     1     0.00         13.84100000       4.43835000       4.34810000

    # find atoms section, everything from the word Atoms to the next line starting with another word
    # or the end of the string
    search_result = re.search('\nAtoms.*?\n(.*)?((\n\w)|$)', text, re.DOTALL)
    
    if search_result:
        
        atom_section_text = search_result.group(1)
        
        search_result = re.findall(f'^\s*(\d+)\s+(\d+)\s+(.*)', atom_section_text, re.MULTILINE)

        if search_result:

            # re.findall returns a list of entries, one for each match
            for entry in search_result:

                # entry is a list, each member is a group in the regex pattern
                # 0 = atom index
                # 1 = atom type
                # 2 = x y z coordinates

                atom_index = int(entry[0])

                atom_type = int(entry[1])

                data_search_result = re.findall('[\d|\.]+', entry[2])

                if data_search_result:

                    # Column format: x, y, z, ...
                    charge_index = None
                    x_index = 0

                    if len(data_search_result) == 4:

                        # Column format: charge, x, y, z
                        charge_index = 0
                        x_index = 1
                                   
                    if charge_index != None:
                        charge = float(data_search_result[charge_index])
                        charge_data = np.vstack((charge_data, [charge]))

                    x = float(data_search_result[x_index])
                    y = float(data_search_result[x_index+1])
                    z = float(data_search_result[x_index+2])

                    coords = np.vstack((coords,[atom_index, atom_type, x, y, z]))
                    
        else:
            print('No lines that look like atomic coordinates found.')
    else:
        print('No text that looks like an atomic coordinates section found.')
                
    return (coords, charge_data)

def find_connected_atoms(atom_index, bonds, atom_list):
    """
    Find all other atoms in the molecule containing one specific atom.
    Atoms are added to atom_list.

    Parameters
    ----------
    atom_index : int
        index of an atom
    bonds : numpy array
        array of bonds generated by add_bond_data()
    atom_list : list
        list of atoms in the same molecule as the atom with atom_index

    Returns
    -------

    """
    
    for bond, bond_index in zip(bonds, range(bonds.shape[0])):
        
        # bonds structure
        # column:  0                1           2           3       4          5          6         7  8  9 10
        # bond_index bond_type_number atom1_index atom2_index comment atom1_type atom2_type Boundary: nx ny nz
        
        # if the atom is in this bond, add the other atom to the list
        if int(bond[2]) == atom_index:
            
            atom_list.append(int(bond[3]))
            find_connected_atoms(int(bond[3]), np.delete(bonds, bond_index, axis=0), atom_list)

        elif int(bond[3]) == atom_index:
        
            atom_list.append(int(bond[2]))
            find_connected_atoms(int(bond[2]), np.delete(bonds, bond_index, axis=0), atom_list)

def add_angle_dihedral_data(file, bonds):
    """
    Adds angles and dihedrals to a LAMMPS data file.

    Parameters
    ----------
    file : str
        path to LAMMPS data file
    bonds : numpy array
        array of bonds generated by add_bond_data()

    Returns
    -------
    numpy array
        array of angles

    numpy array
        array of dihedrals

    """
    
    # data structure of bonds
    # col:      0        1           2           3       4          5          6
    # bond_index bond_type atom_index1 atom_index2 comment atom_type1 atom_type2
    #          1         1           1           2 #Na-Cl           1          2
    
    # constants
    atom_index_offset = 2
    atom_type_offset = 5
    
    # columns:
    # angle_index, angle_type, atom1, atom2, atom3, comment
    #           1           1      1      2      3  # 122
    angles = np.empty((0,6))
    
    # columns:
    # angle_type, angle_type_str
    #          1           '122'
    angle_types = np.empty((0,2))
    
    # columns:
    # dihedral_index, dihedral_type, atom1, atom2, atom3, atom4, comment
    #              1              1      1      2      3      4  # 1222
    dihedrals = np.empty((0,7))
    
    # columns:
    # dihedral_type, dihedral_type_str
    #             1           '1222'
    dihedral_types = np.empty((0,2))
    
    for bond1, bond1_index in zip(bonds, range(bonds.shape[0])):
        
        for bond2, bond2_index in zip(bonds[bond1_index+1:], range(bonds[bond1_index+1:].shape[0])):
            
            angle_atoms = None
            dihedral_atoms = None

            for a in range(2):
                for b in range(2):
                    
                    if bond1[atom_index_offset + a] == bond2[atom_index_offset + b]:
                        
                        angle_atoms = [bond1[atom_index_offset + abs(a-1)],
                                       bond1[atom_index_offset + a],
                                       bond2[atom_index_offset + abs(b-1)]]
                        
                        angle_type_str = str(bond1[atom_type_offset + abs(a-1)]) + \
                                         str(bond1[atom_type_offset + a])        + \
                                         str(bond2[atom_type_offset + abs(b-1)])
                        
                        # debug
                        #print('angle atoms 1, 2, 3')
                        #print(angle_atoms)
                        
                        # look for a third bond forming a dihedral
                        for bond3, bond3_index in zip(bonds[bond1_index+1:], range(bonds[bond1_index+1:].shape[0])):
                            
                            if bond3_index != bond2_index:
                                
                                for a in range(2):
                                    for b in range(2):
                                        
                                        if bond2[atom_index_offset + a] == bond3[atom_index_offset + b]:
                                            
                                            #print('Found third bond for dihedral')
                                            
                                            # Check the 3 bonds are linear, not branched
                                            if bond2[atom_index_offset + a] != angle_atoms[1]:
                                            
                                                dihedral_atoms = angle_atoms
                                                
                                                # append is in-place, no assignment necessary
                                                dihedral_atoms.append(bond3[atom_index_offset + abs(b-1)])
                                                
                                                dihedral_type_str = angle_type_str + str(bond3[atom_type_offset + abs(b-1)])
                                                
                                                # debug
                                                #print(f'Dihedral atoms: {dihedral_atoms}')
                                                
                                                break
                            
                                
                        break
                
            if angle_atoms:
                
                # get the angle type number from the string of atom types
                angle_type_index = np.argwhere(angle_types == angle_type_str)
                
                # angle types have two equivalent configurations i.e. 'N-C-C' == 'C-C-N'
                if not angle_type_index.any():
                    angle_type_index = np.argwhere(angle_types == angle_type_str[::-1])

                if angle_type_index.any():
                    # this angle type has already been found
                    angle_type = angle_types[angle_type_index[0][0]][0]
                else:
                    # add a new angle type
                    angle_type = angle_types.shape[0]+1
                    angle_types = np.vstack((angle_types, [angle_type, angle_type_str]))
                    
                # debug
                #print(f'Adding angle: {angle_atoms}, type = {angle_type} ({angle_type_str})')

                angles = np.vstack((angles, [angles.shape[0]+1,
                                             angle_type,
                                             angle_atoms[0],
                                             angle_atoms[1],
                                             angle_atoms[2],
                                             f'# {angle_type_str}']))
                
            if dihedral_atoms:
                
                # get the dihedral type number from the string of atom types
                # e.g. the 4 atoms types are '1222'
                dihedral_type_index = np.argwhere(dihedral_types == dihedral_type_str)
                
                # dihedral types have two equivalent configurations i.e. 'NCCC' == 'CCCN'
                if not dihedral_type_index.any():
                    dihedral_type_index = np.argwhere(dihedral_types == dihedral_type_str[::-1])
                    
                if dihedral_type_index.any():
                    # this dihedral type has already been found
                    dihedral_type = dihedral_types[dihedral_type_index[0][0]][0]
                else:
                    # add new dihedral type
                    dihedral_type = dihedral_types.shape[0] + 1
                    dihedral_types = np.vstack((dihedral_types, [dihedral_type, dihedral_type_str]))
                    
                # debug
                #print(f'Adding dihedral: {dihedral_atoms}, type = {dihedral_type} ({dihedral_type_str})')
                
                dihedrals = np.vstack((dihedrals, [dihedrals.shape[0] + 1,
                                                   dihedral_type,
                                                   dihedral_atoms[0],
                                                   dihedral_atoms[1],
                                                   dihedral_atoms[2],
                                                   dihedral_atoms[3],
                                                   f'# {dihedral_type_str}']))
            
    with open(file, 'a') as data_file:
        data_file.write('\nAngles\n\n')
        np.savetxt(data_file, angles, fmt='%s')
        
        data_file.write('\nDihedrals\n\n')
        np.savetxt(data_file, dihedrals, fmt='%s')

    # Add number of angles, dihedrals to file header
    with open(file, 'r') as data_file:
        file_text = data_file.read()

    file_text_modified = re.sub('\n(\s*)(\d+)(\s+)(bonds.*?)',
                                f'\n\g<1>\g<2>\g<3>\g<4>\n\g<1>{angles.shape[0]}\g<3>angles\n\g<1>{dihedrals.shape[0]}\g<3>dihedrals',
                                file_text)

    if len(file_text_modified) == len(file_text):
        print('Could not add number of bonds to file header')
    else:
        with open(file, 'w') as data_file:
            data_file.write(file_text_modified)
        
    print(f'Added {angles.shape[0]} angles to {file}')
    print(f'Added {dihedrals.shape[0]} dihedrals to {file}')
        
    return angles, dihedrals