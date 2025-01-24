# RECONEXUS Object Oriented Script

## RECONX
This object oriented method is applied to 3D data files obtained from modified RECONX (Glocer et al, 2015). RECONX was modified to define near-null points, trace them into the solar wind, and output the trace results as readable 3D data files. RECONX is a private git repository and requires approval from original developer for access. 

## Main RECONEXUS Script (reconx.py)
The current iteration of the RECONXUS method is the python script `reconx.py` which reads a null file corresponding to a simulation time as a NullGroup object, identifies all positive and negative null pairs, which are represented as NullPair objects within the group. For each NullPair, the algorithm finds the associated trace files, retrieves the relevant IMF conditions from the SWMF run directory, and calculates the geoeffective potential. Additionally, it retrieves CPCP from the ionosphere file, all stored as a key within the NullPair object for further comparison.

## Test Files and Functionality

### Test Files
Test files can be found in the directory `Test`. RECONX null and trace output includes a negative null file `NegNulls_e20150321-020000.dat` and a positive null file `PlusNulls_e20150321-020000.dat` and there associated field line traces `null_line_neg_o01_001_e20150321-020000.dat` and `null_line_pls_o01_001_e20150321-020000.dat`. 

### reconx.py Usage

#### NullGroup
NullGroup is the main function for manipulating null and trace objects. Given a null file the algorithm finds the associated null pair, trace files, the relevant IMF conditions from the SWMF run directory (or given imffile), calculates the geoeffective potential, and stores all relevent information as keys.

example usage of class: a = NullGroup('NegNulls_e20150321-020000.dat', 'rundir', 'imffile')
Some functionality:
check keys: a.keys()
To call calculated geopotential: a[0]['geopot']


class NullGroup(list)

    Class for handling magnetic null pairs and all associated information.

    To instantiate,

    Paramters
    =========
    nullfile : str
        Path to null file, e.g., "NegNulls_e20150321-040000.dat". These files
        contain a list of nulls taken at a certain time.
    rundir : str, defaults to None
        If given, associates a SWMF run directory with this null set.
    imffile : str, defaults to None
        If given, associates an SWMF IMF input file with this null set.
        Will be used to set solar wind conditions when calculating potential
        drops across lines.

#### read_nulls
example usage: `read_nulls(NegNulls_e20150321-020000.dat, reorder=False)` - can read trace and null files

    Read a null file and return a dictionary-like data object.

    Parameters
    ----------
    filename : str
        Name of null file to open.
    reorder : bool, defaults to False
        Reorder the  points in the file such that they represent, in order,
        a trace along a magnetic field line.

    Returns
    -------
    data : dictionary
        Dictionary of nulls and null information.
    
example usage of class: NullPair(dict):

    Class for handling a single Null Pair.

    Parameters
    ----------
    x_pos, x_neg : dmarray
        The XYZ position of the positive and negative nulls, respectively
    inull : integer
        The number of the null pair.
    usw, bsw : 3-element numpy vectors
        The solar wind velocity and magnetic field corresponding to the time
        that the null was found. Units should be km/s and nT, respectively.
    path : str, defaults to empty string
        Path to file location, used for opening line files.
    time : datetime.datetime, defaults to None
        Date time corresponding to point in time when null found.
    iline : int, defaults to 1
        If there are more than 1 lines associated with this null, this
        kwarg sets which to read.
    '''

## First iteration of RECONXEXUS
The python scripts listed below are a set of original stand-alone code used to derive reconnection potential from geoeffective length (GEL). These scripts have since been iterated on to leverage an object oriented approach (see above).
### Scripts and Functionality
`get_footprints.py` 
`get_geolength.py`
`get_iono_drop.py`
