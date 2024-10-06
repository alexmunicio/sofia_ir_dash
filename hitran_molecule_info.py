from hapi import *
import pandas as pd

DEFAULT_STORAGE_DIR = "./HITRAN_Data"
db_begin(DEFAULT_STORAGE_DIR)

# HAPI reference temperature
REF_TEMP = 296

# add to convert from celsius to Kelvin
CELSIUS_TO_KELVIN = 273.15

# converts altitude from feet to kilometers
FEET_TO_KILOMETERS = 1/3280.84

# thermodynamic constant
C2 = 1.4387769


# function pulls data from HITRAN files and returns a dataframe (make sure the data has been fetched first)
def get_hitran_molecule_info(molecule_name, wavenumber_range = None):

    # intialize parameters
    param = ('nu','sw','gamma_air','gamma_self','local_iso_id','molec_id','elower') # the select() function uses a tuple (I decided to not use it)
    param_list = ["nu", "sw", "gamma_air", "gamma_self", "local_iso_id", "elower", "molec_id"] # the getColumns() function uses a list (I used this instead)

    # intialize lists containing HAPI information
    nu       = []
    wl       = []
    S        = []
    gair     = []
    gsel     = []
    iso_id   = []
    elower   = []
    molec_id = []

    nu, S, gair, gsel, iso_id, elower, molec_id = getColumns(molecule_name, param_list)
    wl = 10000.0 / nu

    # turn lists into a pandas dataframe
    df = pd.DataFrame({
        "wavenumber": nu,
        "wavelength": wl,
        "ref_trans_strength": S,
        "gamma_air": gair,
        "gamma_self": gsel,
        "iso_id": iso_id,
        "elower": elower,
        "molec_id": molec_id,
    })

    # check to see whether or not a specific wavenumber range was given to filter the dataframe
    if wavenumber_range is not None:

        # create filter to only include rows that fall within the specified wavenumber range
        conditions = (min(wavenumber_range) <= df["wavenumber"]) & (df["wavenumber"] <= max(wavenumber_range))
        
        # return the df with the filter applied
        return df[conditions]

    # if no wavenumber is specified, then the entire dataframe is returned (it might run very slowly)
    else:
        return df

# Calls "get_hitran_molecule_info" and adds a column to the returned dataframe.
# New column contains the transition strengths for a given temperature (experimental_temp)
def get_transition_strength_for_temp(molecule_name, experimental_temp, wavenumber_range = None): 

    # get molecule info
    molecule_df = get_hitran_molecule_info(molecule_name, wavenumber_range)

    # check to make sure that molecule_df is not empty
    if molecule_df.empty:
        raise MoleculeDataNotFound(f"Molecule information {molecule_name} information is not found in HITRAN_Data .data files", 404)
    
    # initialize partition sum columns for given temperatures, molecules, and isotopolouge
    molecule_df["Q_ref"] = molecule_df.apply(lambda row: partitionSum(row["molec_id"], row["iso_id"], REF_TEMP), axis = 1)

    # initialize experimental temperature in Kelvin
    experimental_temp_kelvin = experimental_temp + CELSIUS_TO_KELVIN
    molecule_df["Q"] = molecule_df.apply(lambda row: partitionSum(row["molec_id"], row["iso_id"], experimental_temp_kelvin), axis = 1)

    # calculate experimental transition strengths for the given experimental temperature
    term_1 = molecule_df["ref_trans_strength"]
    term_2 = molecule_df["Q_ref"] / molecule_df["Q"]
    term_3 = exp(-C2 * molecule_df["elower"] / experimental_temp_kelvin) / exp(-C2 * molecule_df["elower"] / REF_TEMP)
    term_4 = (1 - exp(-C2 * molecule_df["wavenumber"] / experimental_temp_kelvin)) / (1 - exp(-C2 * molecule_df["wavenumber"] / REF_TEMP))

    molecule_df["exp_trans_strength"] = term_1 * term_2 * term_3 * term_4

    return molecule_df

class MoleculeDataNotFound(Exception):
    def __init__(self, message, error_code):
        super().__init__(message)
        self.error_code = error_code