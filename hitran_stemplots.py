import pandas as pd
import plotly.graph_objects as go

from molecular_transition_strength import get_transition_strength_for_location
from molecular_transition_strength import MOLECULE_CONFIG, ISOTOPOLOGUE_CONFIG, ISOTOPOLOGUE_COLOR_CONFIG

# READ ME:
# As far as I am currently aware, there is no decent way to make a stemplot in the Plotly library, without making each stem it's own individual plot object.
# This however, makes the program run really slowly and blows up the graph's legend with a bunch of repetitive info.
#
# Plan: The two functions below are going to subvert this issue, by modifying dataframes to be compatible with a regular scatterplot function
# 
# They will modify the dataframe, such that the data can be put into a single plot object, while also making the plot look like separate stems.
# 
# Plotly.go will attempt to connect all of the stems together, because it is designed to connect every point in the same plot object.
# This should be fine, as long as we control how Plotly connects those points.
# 
# We will do this, by adding points in between each of actual data points in the dataframe.
# By adding a point at zero before and after the actual data points, we can make the graph appear like a group os distinct stemplots
# Each of the vertical lines will be connected horizontally at the x-axis, and thus be hidden from main part of the graph that the user is focused on

# TLDR: this function puts a row before and after each row in a dataframe, with a zero in the specified column for each row
def modify_row_in_dataframe_for_graphing_stemplot(row, y_coordinate_column):
    
    # make a copy of the row
    actual_row = row.copy()
    behind_row = row.copy()
    front_row = row.copy()

    # create a duplicate row, behind the given row

    # assign a value of zero to the column that hais meant for y coordinate values
    # since this is a dataframe object, we have to specify the row number that we want to change (even though there's only one row)
    behind_row[y_coordinate_column] = 0
    front_row[y_coordinate_column] = 0

    # repeat the process for the row in front of the original row

    return pd.DataFrame([behind_row, actual_row, front_row])

# TLDR: this function takes "modify_row_in_dataframe_for_graphing_stemplot()" and applies it to every row in an inputted dataframe
def modify_dataframe_for_graphing_stemplot(df, y_coordinate_column):
    return pd.concat(list(df.apply(lambda row: modify_row_in_dataframe_for_graphing_stemplot(row, y_coordinate_column), axis = 1)))


def get_isotopologue_dfs_from_molecule_transition_strengths_df(molecule_name, experimental_temp, altitude_km, latitude, wavenumber_range = None, cutoff = None, plotly_stemplot = True):
    isotopologue_dataframes = {}

    for iso_id in ISOTOPOLOGUE_CONFIG[molecule_name]:

        molecule_df = get_transition_strength_for_location(molecule_name, experimental_temp = experimental_temp, altitude_km = altitude_km, latitude = latitude, wavenumber_range = wavenumber_range, cutoff = cutoff)

        if molecule_df is None or molecule_df.empty:

            isotopologue_dataframes[iso_id] = None
            continue

        condition = (molecule_df["iso_id"] == ISOTOPOLOGUE_CONFIG[molecule_name][iso_id]) & (molecule_df["col_den_trans"] >= cutoff)
        filtered_molecule_df = molecule_df[condition]

        if not plotly_stemplot:
             
             isotopologue_dataframes[iso_id] = filtered_molecule_df

        else:
            
            if not filtered_molecule_df.empty: 
                
                isotopologue_dataframes[iso_id] = modify_dataframe_for_graphing_stemplot(filtered_molecule_df, "col_den_trans")
    
    return isotopologue_dataframes

def get_isotopologues_as_trace_object_stemplots(temperature, altitude_km, latitude, wavenumber_range = None, cutoff = 1e-4):

    hitran_list = []

    for molecule in MOLECULE_CONFIG:

        print("molecule:", molecule)

        iso_dict = get_isotopologue_dfs_from_molecule_transition_strengths_df(molecule, experimental_temp = temperature, altitude_km = altitude_km, latitude = latitude, wavenumber_range = wavenumber_range, cutoff = cutoff)
        
        for isotopologue, iso_df in iso_dict.items():

            print("isotopologue:", isotopologue)
            print(iso_df)

            if iso_df is None:
                print(f"there are no signficant transition strengths for {isotopologue} in the specified range of wavenumbers")
                continue

            hitran_list.append(
                    go.Scatter(
                    x=iso_df["wavenumber"],
                    y=iso_df["col_den_trans"],
                    mode = "lines+markers",
                    marker={"size": 3},
                    name = isotopologue,
                    yaxis = "y2",
                    line = {"color": ISOTOPOLOGUE_COLOR_CONFIG[isotopologue]},
                )
            )
            
            print(isotopologue, "trace appended successfully")

    return hitran_list