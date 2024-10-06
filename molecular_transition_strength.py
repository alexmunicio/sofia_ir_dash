import pandas as pd
from hitran_molecule_info import get_transition_strength_for_temp
from hitran_molecule_info import MoleculeDataNotFound
from atmospheric_info import get_atmosphere_info_for_altitude



O3_LATITUDE_CONFIG = ["O3_LAT_9", "O3_LAT_36", "O3_LAT_43", "O3_LAT_56"]

MOLECULE_CONFIG = {
    "H2O": {
        "is_lat_dependent": False,
        "is_mixed_gas": False,
        "column_density_header": "H2O",
        "concentration": 1,
        "iso_ids": [1,2,3,4,5,6,7],
    },
    "O3": {
        "is_lat_dependent": True,
        "column_density_header": O3_LATITUDE_CONFIG,
        "concentration": 1,
        "iso_ids": [1,2,3,4,5]
    },
    "HNO3": {
        "is_lat_dependent": False,
        "column_density_header": "MIX",
        "concentration": 10/1e9,
        "iso_ids": [1,2]
    },
    "N2O": {
        "is_lat_dependent": False,
        "column_density_header": "MIX",
        "concentration": 0.28/1e6,
        "iso_ids": [1,2,3,4,5]
    },
    "O2": {
        "is_lat_dependent": False,
        "column_density_header": "MIX",
        "concentration": 2.1/1e6,
        "iso_ids": [1,2,3]
    },
    "CO": {
        "is_lat_dependent": False,
        "column_density_header": "MIX",
        "concentration": 0.75/1e6,
        "iso_ids": [1,2,3,4,5,6]
    },
    "CO2": {
        "is_lat_dependent": False,
        "column_density_header": "MIX",
        "concentration": 400/1e6,
        "iso_ids": [1,2,3,4,5,6,7,9,0,11,12]
    },
    "CH4": {
        "is_lat_dependent": False,
        "column_density_header": "MIX",
        "concentration": 1.6/1e6,
        "iso_ids": [1,2,3,4]
    }
}

ISOTOPOLOGUE_CONFIG = {
    "H2O": {
        "H2O": 1,
        "H2(18O)": 2,
        "H2(17O)": 3,
        "HD(16O)": 4,
        "HD(18O)": 5,
        "HD(17O)": 6,
        "D2(16O)": 7
    },
    "O3": {
        "O3": 1,
        "O2(18O)": 2,
        "O(18O)O": 3,
        "O2(17O)": 4,
        "O(17O)O": 5
    },
    "HNO3": {
        "HNO3": 1,
        "H(15N)O3": 2
    },
    "N2O": {
        "N2O": 1,
        "N(15N)O": 2,
        "(15N)NO": 3,
        "N2(18O)": 4,
        "N2(17O)": 5
    },
    "O2": {
        "O2": 1,
        "O(18O)": 2,
        "O(17O)": 3
    },
    "CO": {
        "CO": 1,
        "(13C)O": 2,
        "C(18O)": 3,
        "C(17O)": 4,
        "(13C)(18O)": 5,
        "(13C)(17O)": 6
    },
    "CO2": {
        "CO2": 1,
        "(13C)O2": 2,
        "OC(18O)": 3,
        "OC(17O)": 4,
        "O(13C)(18O)": 5,
        "O(13C)(17O)": 6,
        "C(18O)2": 7,
        "(17O)C(18O)": 8,
        "C(17O)2": 9,
        "(13C)(18O)2": 0,
        "(18O)(13C)(17O)": 11,
        "(13C)(17O)2": 12
    },
    "CH4": {
        "CH4": 1,
        "(13C)H4": 2,
        "CH3D": 3,
        "(13C)H3D": 4
    }
}

ISOTOPOLOGUE_COLOR_CONFIG = {
        "H2O": "blue",
        "H2(18O)": "teal",
        "H2(17O)": "navy",
        "HD(16O)": "rgb(0, 83, 146)",
        "HD(18O)": "aqua",
        "HD(17O)": "midnight",
        "D2(16O)": "sky",

        "O3": "rgb(0, 99, 26)",
        "O2(18O)": "rgb(0, 176, 136)",
        "O(18O)O": "rgb(5, 170, 107)",
        "O2(17O)": "rgb(56, 124, 98)",
        "O(17O)O": "rgb(96, 172, 2)",

        "HNO3": "grey",
        "H(15N)O3": "grey",

        "N2O": "orange",
        "N(15N)O": "orange",
        "(15N)NO": "orange",
        "N2(18O)": "orange",
        "N2(17O)": "orange",

        "O2": "purple",
        "O(18O)": "indigo",
        "O(17O)": "purple",

        "CO": "green",
        "(13C)O": "green",
        "C(18O)": "green",
        "C(17O)": "green",
        "(13C)(18O)": "green",
        "(13C)(17O)": "green",

        "CO2": "rgb(234, 95, 2)",
        "(13C)O2": "rgb(222, 68, 125)",
        "OC(18O)": "plum",
        "OC(17O)": "rgb(225, 20, 64)",
        "O(13C)(18O)": "rgb(255, 0, 0)",
        "O(13C)(17O)": "rgb(224, 0, 0)",
        "C(18O)2": "rgb(174, 0, 0)",
        "(17O)C(18O)": "rgb(134, 0, 0)",
        "C(17O)2": "rgb(94, 0, 0)",
        "(13C)(18O)2": "rgb(85, 0, 0)",
        "(18O)(13C)(17O)": "rgb(76, 0, 0)",
        "(13C)(17O)2": "rgb(65, 0, 0)",

        "CH4": "brown",
        "(13C)H4": "rgb(150, 121, 105)",
        "CH3D": "rgb(180, 110, 13)",
        "(13C)H3D": "rgb(142, 99, 13)"
}

# total_isotopologue_count = 0
# for molecule, isotopologue_dict in ISOTOPOLOGUE_CONFIG.items():
#     for isotopologue in isotopologue_dict:
#         total_isotopologue_count += 1

# multiplies column density for experiment altitude and creates a new column in the dataframe returned by "get_transition_strength_for_temp"
def get_transition_strength_for_location(molecule_name, experimental_temp, altitude_km, latitude = None, wavenumber_range = None, cutoff = None) : 
    
    # initialize dataframe with column for temperature dependent transition strength
    try:
        molecule_df = get_transition_strength_for_temp(molecule_name, experimental_temp, wavenumber_range)

    # if the molecule data is not found (this is a normal occurrence when the wavelength is filtered)    
    except MoleculeDataNotFound as e:
        print(f"Caught an exception: {e}")
        print("Don't worry: this probably means that there just aren't any transitions for that specific molecule in the wavenumber range that was specified")
        return None

    
    # initialize dataframe row for specific altitude
    atm_info = get_atmosphere_info_for_altitude(altitude_km)

    # initialize molecule concentration
    molecular_concentration = MOLECULE_CONFIG[molecule_name]["concentration"]
    
    # check to see if the specific molecule is latitude dependent (essentially, whether or not it's ozone)
    if MOLECULE_CONFIG[molecule_name]["is_lat_dependent"]:

        # initialize the column densities for each latitude
        latitude_dependent_df = atm_info[MOLECULE_CONFIG[molecule_name]["column_density_header"]]

        # initialize column density based on what the latitudinal coordinate is
        if latitude <= 9:
            column_density = latitude_dependent_df["O3_LAT_9"]
        elif 9 < latitude <= 36:
            column_density = latitude_dependent_df["O3_LAT_36"]
        elif 36 < latitude <= 43:
            column_density = latitude_dependent_df["O3_LAT_43"]
        elif 43 < latitude <= 56:
            column_density = latitude_dependent_df["O3_LAT_56"]
        else:
            print("ozone column densities are not listed above 56 degrees")
            return None

        # calculate expected transition strength based on column density and molecular concetration
        molecule_df["col_den_trans"] = molecule_df["exp_trans_strength"] * column_density * molecular_concentration

    # if the molecule is not latitude dependent (i.e. it's not ozone), then we just calculate the new transition strength the easy way
    else: 
        # initialize column density value
        column_density = atm_info[MOLECULE_CONFIG[molecule_name]["column_density_header"]]
        molecule_df["col_den_trans"] = molecule_df["exp_trans_strength"] * column_density * molecular_concentration
    
    if cutoff:
        cutoff_condition = molecule_df >= cutoff
        filtered_molecule_df = molecule_df[cutoff_condition]
    
    # return the final dataframe
    return molecule_df
