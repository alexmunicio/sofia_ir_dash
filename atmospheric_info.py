import pandas as pd

ATMOSPHERIC_DATA_EXCEL_SHEET = "Model Atmosphere Table.xls"

MODEL_ATMOSPHERE_HEADER_RECONFIG = {
    "Col. Den. H20 (MOL/CM^2)": "H2O",
    "Col. Den. MIX GAS (MOL/CM^2)": "MIX",
    "Col. Den. 03 (MOL/CM^2)": "O3_LAT_9",
    "Unnamed: 10": "O3_LAT_36",
    "Unnamed: 11": "O3_LAT_43",
    "Unnamed: 12": "O3_LAT_56"}

def get_atmosphere_dataframe():

    # initialize pandas dataframe from excel sheet
    excel_df = pd.read_excel(ATMOSPHERIC_DATA_EXCEL_SHEET)

    # remove the first row, because it doesn't actually contain data (it just exists because of the excel format)
    atmosphere_df = excel_df.iloc[1:,:]

    # return atmosphere_df dataframe but with simpler headers
    return atmosphere_df.rename(columns = MODEL_ATMOSPHERE_HEADER_RECONFIG)

# calls "get_atmosphere_dataframe()" to pull column density for a given molecule, altitude, and latitudinal coordinate
def get_atmosphere_info_for_altitude(altitude):

    # initialize dataframe containing atmospheric data
    atm_df = get_atmosphere_dataframe()

    # filter for rows with 
    filtered_df = atm_df[atm_df["Alt (KM)"] <= altitude]
    # #return atm_df[filtered_df["Alt (KM)"]].
    # print()
    return filtered_df.loc[filtered_df["Alt (KM)"].idxmax()]

    #return atm_df.loc[idx]
# pulls a column density dataframe from "model atmosphere.xls"
