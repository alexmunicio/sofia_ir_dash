from hapi import *
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

DEFAULT_STORAGE_DIR = "./HITRAN_Data"

# See: https://www.hitran.org/docs/molec-meta/
CONF = [
    {
        "table_name": "H2O",
        "hitran_molecule_number": 1,
        "hitran_isotope_list": [1, 2, 3, 4, 5, 6, 129],
    },
    {
        "table_name": "CH4",
        "hitran_molecule_number": 6,
        "hitran_isotope_list": [32, 33, 34, 35],
    },
    {
        "table_name": "CO",
        "hitran_molecule_number": 5,
        "hitran_isotope_list": [26, 27, 28, 29, 30, 31],
    },
    {
        "table_name": "CO2",
        "hitran_molecule_number": 2,
        "hitran_isotope_list": [7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122],
    },
    {
        "table_name": "O3",
        "hitran_molecule_number": 3,
        "hitran_isotope_list": [16, 17, 18, 19, 20],
    },
    {
        "table_name": "HNO3",
        "hitran_molecule_number": 12,
        "hitran_isotope_list": [47, 117],
    },
    {
        "table_name": "N2O",
        "hitran_molecule_number": 4,
        "hitran_isotope_list": [21, 22, 23, 24, 25],
    },
    {
        "table_name": "O2",
        "hitran_molecule_number": 7,
        "hitran_isotope_list": [36, 37, 38],
    },
]

# initialize minimum and maximum wavenumber
MIN_WAVENUMBER = 0
MAX_WAVENUMBER = 2000

# functino fetches molecule information from Hitran by using "fetch_by_ids" and CONFIG as arguments
def fetch_all_molecules_from_hitran(dir = DEFAULT_STORAGE_DIR, minWN=MIN_WAVENUMBER, maxWN=MAX_WAVENUMBER):

    # function signifies the directory where the molecule data will be stored
    db_begin(dir)

    # loop iterates through the molecule dictionaries in CONF
    # applies the
    for molecule in CONF:
        fetch_by_ids(
            molecule["table_name"], molecule["hitran_isotope_list"], minWN, maxWN
        )


if __name__ == "__main__":
    fetch_all_molecules_from_hitran()

# for i in CONF:
#     # initialize query arguments for the select function
#     param = ("nu", "sw", "gamma_air", "gamma_self", "local_iso_id", "molec_id", "elower")
#     cond = ("between", "nu", 705, 750)
#     dir = "hapi_temp_" + i["table_name"] + ".csv"
#     print(i["table_name"])

#     # query hapi using the select function with the query arguments
#     select(i["table_name"], ParameterNames = param, Conditions = cond, File = dir)
