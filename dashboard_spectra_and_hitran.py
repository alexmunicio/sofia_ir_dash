import os
from dash import Dash, html, dcc, callback, Output, Input, State, dash_table
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import pandas as pd

# pd.set_option("display.max_rows", None)
import numpy as np
from scipy.signal import find_peaks
from exes_info import get_exes_file_data
from exes_info import get_endian_modified_exes_file_data
from hitran_stemplots import get_isotopologues_as_trace_object_stemplots

# MOLECULE_LIST = list(MOLECULE_CONFIG.keys())
MAP_CONFIG = {
    "scope": "north america",
    # "lonaxis": dict(range=[-155, -84], showticklabels = True),
    # "lataxis": dict(range=[10, 55], showticklabels = True),
    "lonaxis": {
        "range": [-155, -84],
        "showgrid": True,
        "gridcolor": "Black",
        "dtick": 5,
        "gridwidth": 0.25,
    },
    "lataxis": {
        "range": [10, 55],
        "showgrid": True,
        "gridcolor": "Black",
        "dtick": 5,
        "gridwidth": 0.25,
    }
}
MARKERSIZE = 2
EXTENSION = ".fits"
DIRECTORY = "EXES_Files"
# DIRECTORY = "AFGL_2136_search"

currently_loaded_data = {
    "file_name": None,
    "exes": None,
    "hitran": [],
    "smooth_width": None,
    "cutoff": None,
}

def get_all_fits_geographic_data(dir):

    geographic_info = []

    filelist = [file for file in os.listdir(dir) if file.endswith(EXTENSION)]

    for file in filelist:
        file_info = get_exes_file_data(file, dir=DIRECTORY)
        file_info["file_name"] = file
        geographic_info.append(file_info)

    return pd.DataFrame(geographic_info)

FITS_GEOGRAPHIC_INFO = get_all_fits_geographic_data(DIRECTORY)

def empty_spectra(graph_title):
    trace1 = go.Scatter(x=[], y=[], mode="markers", name="atran data")

    trace2 = go.Scatter(
        x=[], y=[], mode="markers", name="experimental data", yaxis="y2"
    )

    layout = go.Layout(
        title=graph_title,
        yaxis={"title": "y1 axis"},
        yaxis2={"title": "y2 axis", "overlaying": "y", "side": "right"},
    )

    # return new figure object
    return {"data": [trace1, trace2], "layout": layout}

def update_plot_axes_ranges(trace_list, selectedData):
    
    # if the dictionary selectedData contains the key, "lassoPoints", then it is data from the lasso tool
    if "lassoPoints" in selectedData:

        # if I had to guess, these are probably the coordinates outlining the lasso path (for some reason there's two sets of y coordinates)
        selected_lasso_points = pd.DataFrame(selectedData["lassoPoints"])

        # when the lasso tool is used the selected data points will be used to set the x and y bounds
        # the maximum/minimum x and y values will be used to set the x and y bounds
        x_lower = selected_lasso_points["x"].min()
        x_upper = selected_lasso_points["x"].max()
        y1_lower = selected_lasso_points["y"].min()
        y1_upper = selected_lasso_points["y"].max()
        y1_lower = selected_lasso_points["y2"].min()
        y1_upper = selected_lasso_points["y2"].max()

    # if the data was not from the lasso tool, then it is handled as data from the box select tool
    else:

        # intialize range dictionary containing the bounds for the box select tool
        selected_range = selectedData["range"]

        # take x, y1, and y2 bounds from the box select tool's range data
        x_lower = selected_range["x"][0]
        x_upper = selected_range["x"][-1]
        
        y1_lower = selected_range["y"][0]
        y1_upper = selected_range["y"][-1]

        # y2_lower = selected_range["y2"][0]
        if "y2" in selected_range:
            y2_upper = selected_range["y2"][-1]
        
        else:
            y2_upper = False

    for trace in trace_list:

        df = pd.DataFrame({
            "x": trace["x"],
            "y": trace["y"]
        })
        
        x_bounds = (x_lower <= df["x"]) & (df["x"] <= x_upper)
        y_bounds = True

        if trace["yaxis"] == "y":

            y_bounds = (y1_lower <= df["y"]) & (df["y"] <= y1_upper)
        
        if (trace["yaxis"] == "y2") and (y2_upper):
        
            # y_bounds = (y2_lower <= df["y"]) & (df["y"] <= y2_upper)
            y_bounds = df["y"] <= y2_upper
        

        filtered_df = df[x_bounds & y_bounds]

        trace["x"] = filtered_df["x"]
        trace["y"] = filtered_df["y"]

    return trace_list

def estimate_fwhm_for_exes_peak(exes_df, peak_heights_df):
    return None

def get_all_spectra_data(exes_file_name, smooth_width, cutoff):

    if (currently_loaded_data["file_name"] is None) or (currently_loaded_data["file_name"] != exes_file_name):

        currently_loaded_data["file_name"] = exes_file_name
        currently_loaded_data["exes"] = get_endian_modified_exes_file_data(exes_file_name, dir = DIRECTORY, smooth_width = smooth_width)
        currently_loaded_data["smooth_width"] = smooth_width

        temperature = currently_loaded_data["exes"]["temperature"]
        altitude_km = currently_loaded_data["exes"]["avg_altitude_km"]
        latitude = currently_loaded_data["exes"]["latitude"]
        wavenumbers = currently_loaded_data["exes"]["wavenumber"]

        currently_loaded_data["hitran"] = get_isotopologues_as_trace_object_stemplots(temperature, altitude_km, latitude, wavenumbers, cutoff)
        
    if currently_loaded_data["smooth_width"] != smooth_width:

        currently_loaded_data["exes"] = get_endian_modified_exes_file_data(exes_file_name, dir = DIRECTORY, smooth_width = smooth_width)

    # if currently_loaded_data["cutoff"] != cutoff:
        
    #     currently_loaded_data["hitran"] = get_isotopologues_as_trace_object_stemplots(temperature, altitude_km, latitude, wavenumbers)

    return currently_loaded_data["exes"], currently_loaded_data["hitran"]


app = Dash()

app.layout = [
    html.H1(children="SOFIA Experiment Dashboard", style={"textAlign": "center"}),
    dcc.Dropdown(FITS_GEOGRAPHIC_INFO["object"].sort_values().unique(), id = "observation_selection", multi = True),
    dcc.Graph(id="exps_map", config={"scrollZoom": False}),
    dcc.Input(id = "smooth_width_parameter", type = "number", placeholder = "EXES Spectra Smooth Width", value = 9),
    dcc.Input(id = "hitran_cutoff_parameter", type = "number", placeholder = "Hitran Transition Strength Cutoff", value = 1e-4),
    dcc.Graph(id="exp_spectra", config={"displayModeBar": True, "modeBarButtonsToAdd": ["select2d", "lasso2d"]},),
    html.H3(children="Line of Best Fit Tuning Parameters"),
    dcc.Input(id="baseline_parameter", type="number", placeholder="spectra baseline", value=1),
    dcc.Input(id="height_parameter", type="number", placeholder="height input", value=0.9),
    dcc.Input(id="prominence_parameter", type="number", placeholder="prominence input"),
    dcc.Input(id="distance_parameter", type="number", placeholder="distance input"),
    dcc.Graph(id="spectra_peaks", config={}),
    dash_table.DataTable(
        id="dynamic_table",
        data=[],
        columns=[
            {"id": "wavenumber", "name": "wavenumber"},
            {"id": "flux", "name": "flux"},
            {"id": "peak_model", "name": "peak_model", "presentation": "dropdown"},
        ],
        editable=True,
        dropdown={
            "peak_model": {
                "options": [
                    {"label": "gaussian", "value": "gaussian"},
                    {"label": "lorentzian", "value": "lorentzian"},
                ]
            }
        },
    ),
    dcc.Graph(id="spectra_fit"),
]


# callback for hoverdata, map formatting, and experiment info
@app.callback(
    Output("exps_map", "figure"),
    # Input("exps_map", "clickData"),
    Input("observation_selection", "value")
)
def on_geo_map_click(dropdown_data):

    if dropdown_data is not None:
        dropdown_data_condition = FITS_GEOGRAPHIC_INFO["object"].isin(dropdown_data)
    else:
        dropdown_data_condition = FITS_GEOGRAPHIC_INFO == FITS_GEOGRAPHIC_INFO

    fig = px.scatter_geo(
        FITS_GEOGRAPHIC_INFO[dropdown_data_condition],
        lat="latitude",
        lon="longitude",
        hover_name="file_name",
        hover_data={
            "object": True,
            "wavelength_range": True,
            "avg_altitude": True,
            "temperature": True,
        },
        color = "avg_altitude",
        color_continuous_scale = "Viridis",
        projection="mercator",  # scope = "usa" does not have a mercator projection, so nothing is displayed when this parameter is used
        title="SOFIA Experiment Locations",
    )
    fig.update_geos(**MAP_CONFIG)

    fig.update_layout(height=600, width=1600)

    return fig


# callback to plot a specific experiment's spectra, based on the experiment on the map that the user clicks on
@app.callback(
        Output("exp_spectra", "figure"),
        Input("exps_map", "clickData"),
        Input("smooth_width_parameter", "value"),
        Input("hitran_cutoff_parameter", "value")
)
def update_graph(clickData, smooth_width, hitran_cutoff):

    # this will check to see if a specific experiment on the map has been clicked and has relevant info
    if (not clickData) or (not smooth_width):

        return empty_spectra("Select an experiment from the map")

    # intialize dictionary from click data containing info for selected experiment
    experiment_info = clickData["points"][0]

    # pull the experiment file name from the dictionary
    experiment_file_name = experiment_info["hovertext"]

    

    # filter the datafrane of FITS file dictionaries to get the user selected experiment
    exes_dict, hitran_traces = get_all_spectra_data(experiment_file_name, smooth_width, hitran_cutoff)
    
    
    spectra_df = exes_dict["dataframe"]

    traces_list = []

    # graph spectra data for user selected experiment
    traces_list.append(go.Scatter(
        x=spectra_df["wavenumber"],
        y=spectra_df["atran"],
        mode="lines+markers",
        marker={"size": MARKERSIZE},
        name="atran data",
        yaxis="y1",
        line={"color": "red"},
    ))

    traces_list.append(go.Scatter(
        x=spectra_df["wavenumber"],
        y=spectra_df["smooth_flux"],
        mode="lines+markers",
        marker={"size": MARKERSIZE},
        name="experimental data",
        yaxis="y1",
        line={"color": "black"},
    ))

    traces_list += hitran_traces

    print("file name:", experiment_file_name)
    print("longitude:", exes_dict["longitude"])
    print("latitude:", exes_dict["latitude"])
    print("altitude:", exes_dict["avg_altitude"])
    print("temperature:", exes_dict["temperature"])
    print("wavenumber range:", exes_dict["wavenumber_range"])


    layout = go.Layout(
        title = experiment_file_name + ":  " + exes_dict["dashboard_spectrum_title"],
        yaxis={"title": "y1 axis"},
        yaxis2={"title": "y2 axis", "overlaying": "y", "side": "right", "type": "log"},
    )

    return {"data": traces_list, "layout": layout}


# callback to plot the user selected portion of the spectra onto a separate graph that will be used for plotting peaks, identified by scipy.signal.find_peaks()
@app.callback(
    # inputs are apparently assigned to the function in the same order as they are written in here
    Output("spectra_peaks", "figure"),
    Input("exp_spectra", "selectedData"),
    Input("exp_spectra", "figure"),
    Input("height_parameter", component_property="value"),
    Input("prominence_parameter", component_property="value"),
    Input("distance_parameter", component_property="value"),
    Input("baseline_parameter", "value"),
)
def update_spectra_peaks(selectedData, parent_spectra_figure, height, prominence, distance, baseline):

    if selectedData and height:

        # initialize figure values, so we can manipulate them to display the user selected data in the new graph
        spectra_layout = parent_spectra_figure["layout"]
        spectra_data = parent_spectra_figure["data"]

        updated_spectra_data = update_plot_axes_ranges(spectra_data, selectedData)

        exes_peaks_df = pd.DataFrame({

            "baseline": [baseline for i in range(len(updated_spectra_data[1]["x"]))],
            "exes_wavenumbers": updated_spectra_data[1]["x"],
            "exes_data": updated_spectra_data[1]["y"]
        
        })
        
        # we'll add in the baseline onto the plot
        spectra_data.append({
                "type": "scatter",
                "mode": "lines",
                "x": exes_peaks_df["exes_wavenumbers"],
                "y": exes_peaks_df["baseline"],
                "name": "Baseline Flux",
                "line": {"dash": "dash", "color": "green"}
        })

        # identify the peaks within filtered spectra dataframe
        # we need to invert the data, so that the dips become local maxima or "peaks"
        peak_indices, peak_heights = find_peaks(
            -1 * exes_peaks_df["exes_data"],
            height=-height,
            prominence=prominence,
            distance=distance,
        )

        # peaks_dataframe = exes_peaks_df.reset_index().loc[peak_indices]
        peaks_dataframe = exes_peaks_df.reset_index().loc[peak_indices]

        # now we must append two new elements into the spectra data list: the points marking where the peaks where identified AND a dashed line for an expected baseline flux
        spectra_data.append(
            {
                "type": "scatter",
                "mode": "markers",
                "x": peaks_dataframe["exes_wavenumbers"],
                "y": peaks_dataframe["exes_data"],
                "name": "Identified Peaks",
                "line": {"color": "red"}
            }
        )


        # manipulate x bounds for the layout object
        spectra_layout["xaxis"]["range"] = [min(exes_peaks_df["exes_wavenumbers"]), max(exes_peaks_df["exes_wavenumbers"])]

        # return new figure dictionary
        return {"data": spectra_data, "layout": spectra_layout}

    else:
        return empty_spectra("Select a section of the spectra with the box select tool or the lasso tool")

@app.callback(
        Output("dynamic_table", "data"),
        Input("spectra_peaks", "figure")
)
def update_table(parent_spectra_figure):

    if parent_spectra_figure["data"][-1]["name"] == "Identified Peaks":

        # initialize figure objects into something easier to read and access
        spectra_data = parent_spectra_figure["data"]

        peaks_info = spectra_data[-1]

        peaks_df = pd.DataFrame({
            "wavenumber": peaks_info["x"],
            "flux": peaks_info["y"],
            "peak_model": ["gaussian" for peak in peaks_info["x"]]
        })
        
        return peaks_df.to_dict("records")

    else:
        return []

@app.callback(
    Output("spectra_fit", "figure"),
    Input("spectra_peaks", "figure"),
    Input("dynamic_table", "data"),
    Input("baseline_parameter", "value")
)
def update_spectra_fits(parent_spectra_figure, table_data, baseline):

    # initialize figure values, so we can manipulate them to display the user selected data in the new graph
    spectra_layout = parent_spectra_figure["layout"]
    spectra_data = parent_spectra_figure["data"]

    if not table_data:
        
        return empty_spectra("Ensure the correct spectra peaks are identified and that the desired peak models are selected")

    exes_df = pd.DataFrame({
        "x": spectra_data[1]["x"],
        "y": spectra_data[1]["y"]
    })

    # intialize max height to determine the Full Width at Half Maximum (FWHM)
    maximum_peak_height_from_baseline = (baseline - exes_df["y"]).max()
    half_maximum = 0.5 * maximum_peak_height_from_baseline # actually no, we need to be doing this on the detected peaks, not the absolute maximum

    # add the first model (using the baseline as the amplitude of our Horizontal Line model) into the summation of models
    summation_model = models.Const1D(baseline)

    for peak_data in table_data:
        
        mean = peak_data["wavenumber"]
        amplitude = -1 * (baseline - peak_data["flux"])
        stddev = 0.025

        if peak_data["peak_model"] == "gaussian":
        
            model = models.Gaussian1D(amplitude = amplitude, mean = mean, stddev = stddev)
        
        else: #peak_data["peak_model"] == "lorentzian":

            model = models.Lorentz1D(amplitude = amplitude, x_0 = mean, fwhm = stddev)

        summation_model += model

    # # plot the summation model for testing purposes
    # spectra_data.append({
    #     "type": "scatter",
    #     "mode": "lines",
    #     "x": exes_df["x"],
    #     "y": summation_model(exes_df["x"]),
    #     "name": "inferred model",
    #     "line": {"color": "purple"}
    # })

    # take line of best fit, using the fitting object and the summation model
    fitter = fitting.TRFLSQFitter()
    fitted_model = fitter(summation_model, exes_df["x"], exes_df["y"])


    # add data for the line of best fit into the spectra data object, so that it can be graphed onto the plot
    spectra_data.append({
        "type": "scatter",
        "mode": "lines",
        "x": exes_df["x"],
        "y": fitted_model(exes_df["x"]),
        "name": "Line of Best Fit",
        "line": {"color": "rgb(127, 46, 231)"}
    })

    return {
        "data": spectra_data,
        "layout": spectra_layout
    }


if __name__ == "__main__":
    app.run_server(debug=True)
