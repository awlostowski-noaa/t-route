#!/usr/bin/env python
# coding: utf-8

"""NHD Network test cases

Test v02 routing on specific test cases

"""
## Parallel execution
import sys
import time
import numpy as np
import argparse
import pathlib
import pandas as pd
from functools import partial
from joblib import delayed, Parallel
from itertools import chain, islice
from operator import itemgetter

## network and reach utilities
import compute_nhd_routing_SingleSeg_v02 as nhd_compute
import troute.nhd_network_utilities_v02 as nnu
import mc_reach
import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io

ENV_IS_CL = False
if ENV_IS_CL:
    root = pathlib.Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = pathlib.Path("../..").resolve()


def build_test_parameters(
    test_name,
    supernetwork_parameters,
    run_parameters,
    output_parameters,
    restart_parameters,
    forcing_parameters,
    parity_parameters,
):

    if test_name == "pocono1":

        print("running test case for Pocono_TEST1 domain - NO RESERVOIRS")

        # File path to WRF Hydro data
        NWM_test_path = pathlib.Path(
            root, "test/input/geo/NWM_2.1_Sample_Datasets/Pocono_TEST1/"
        ).resolve()

        # Simulation domain RouteLink file
        routelink_file = pathlib.Path(
            NWM_test_path, "primary_domain", "DOMAIN", "Route_Link.nc",
        ).resolve()

        # Speficify WRF hydro restart file, name and destination
        time_string = "2017-12-31_06-00_DOMAIN1"
        wrf_hydro_restart_file = pathlib.Path(
            NWM_test_path, "example_RESTART", "HYDRO_RST." + time_string
        ).resolve()

        # specify supernetwork parameters
        supernetwork_parameters = {
            "title_string": "Pocono1_TEST",
            "geo_file_path": routelink_file,
            "columns": {
                "key": "link",
                "downstream": "to",
                "dx": "Length",
                "n": "n",  # TODO: rename to `manningn`
                "ncc": "nCC",  # TODO: rename to `mannningncc`
                "s0": "So",  # TODO: rename to `bedslope`
                "bw": "BtmWdth",  # TODO: rename to `bottomwidth`
                "waterbody": "NHDWaterbodyComID",
                "tw": "TopWdth",  # TODO: rename to `topwidth`
                "twcc": "TopWdthCC",  # TODO: rename to `topwidthcc`
                "musk": "MusK",
                "musx": "MusX",
                "cs": "ChSlp",  # TODO: rename to `sideslope`
            },
            "waterbody_null_code": -9999,
            "terminal_code": 0,
            "driver_string": "NetCDF",
            "layer_string": 0,
        }

        # specity output parameters
        output_parameters["nc_output_folder"] = None

        # specify restart parameters
        restart_parameters["wrf_hydro_channel_restart_file"] = wrf_hydro_restart_file
        restart_parameters["wrf_hydro_channel_ID_crosswalk_file"] = routelink_file
        restart_parameters["wrf_hydro_channel_ID_crosswalk_file_field_name"] = "link"
        restart_parameters[
            "wrf_hydro_channel_restart_upstream_flow_field_name"
        ] = "qlink1"
        restart_parameters[
            "wrf_hydro_channel_restart_downstream_flow_field_name"
        ] = "qlink2"
        restart_parameters["wrf_hydro_channel_restart_depth_flow_field_name"] = "hlink"

        # specify restart parameters
        forcing_parameters["qlat_input_folder"] = pathlib.Path(
            root,
            "test/input/geo/NWM_2.1_Sample_Datasets/Pocono_TEST1/example_CHRTOUT/",
        ).resolve()
        forcing_parameters["qlat_file_pattern_filter"] = "*.CHRTOUT_DOMAIN1"
        forcing_parameters["qlat_file_index_col"] = "feature_id"
        forcing_parameters["qlat_file_value_col"] = "q_lateral"

        # construct time domain (run_parameters)
        # qlat_files = glob.glob(
        # forcing_parameters["qlat_input_folder"] + forcing_parameters["qlat_file_pattern_filter"],
        # recursive=True
        # )
        # ql = nhd_io.get_ql_from_wrf_hydro_mf(
        # qlat_files,
        # forcing_parameters["qlat_file_index_col"]
        # )

        # wrf_time = ql.columns.astype("datetime64[ns]")
        # dt_wrf = (wrf_time[1] - wrf_time[0])
        # sim_duration = (wrf_time[-1] + dt_wrf) - wrf_time[0]

        # dt = 300
        # dt_routing = pd.Timedelta(str(dt) + 'seconds')

        # run_parameters["dt"] = dt
        # run_parameters["nts"] = round(sim_duration / dt_routing)
        # run_parameters["qts_subdivisions"] = dt_wrf/dt_routing

        run_parameters["dt"] = 300
        run_parameters["nts"] = 288
        run_parameters["qts_subdivisions"] = 12

        run_parameters["assume_short_ts"] = True

        parity_parameters["parity_check_input_folder"] = pathlib.Path(
            root,
            "test/input/geo/NWM_2.1_Sample_Datasets/Pocono_TEST1/example_CHRTOUT/",
        ).resolve()
        parity_parameters["parity_check_file_pattern_filter"] = "*.CHRTOUT_DOMAIN1"
        parity_parameters["parity_check_file_index_col"] = "feature_id"
        parity_parameters["parity_check_file_value_col"] = "streamflow"
        parity_parameters["parity_check_compare_node"] = 4186169

    return (
        supernetwork_parameters,
        run_parameters,
        output_parameters,
        restart_parameters,
        forcing_parameters,
        parity_parameters,
    )


def parity_check(parity_parameters, run_parameters, nts, dt, results):
    
    if "parity_check_input_folder" in parity_parameters:
        
        validation_files = pathlib.Path(parity_parameters["parity_check_input_folder"]).rglob(
            parity_parameters["parity_check_file_pattern_filter"]
        )

        # read validation data from CHRTOUT files
        validation_data = nhd_io.get_ql_from_wrf_hydro_mf(
            validation_files,
            parity_parameters["parity_check_file_index_col"],
            parity_parameters["parity_check_file_value_col"],
            # [compare_node],
        )
    
    if "parity_check_file" in parity_parameters:
        validation_data = pd.read_csv(parity_parameters["parity_check_file"], index_col=0)
        validation_data.index = validation_data.index.astype(int)
        validation_data.columns = validation_data.columns.astype("datetime64[ns]")
        validation_data = validation_data.sort_index(axis="index")
        
    compare_node = parity_parameters["parity_check_compare_node"]

    wrf_time = validation_data.columns.astype("datetime64[ns]")
    dt_wrf = wrf_time[1] - wrf_time[0]
    sim_duration = (wrf_time[-1] + dt_wrf) - wrf_time[0]

    # construct a dataframe of simulated flows
    fdv_columns = pd.MultiIndex.from_product([range(nts), ["q", "v", "d"]])
    if run_parameters.get("return_courant", False):
        flowveldepth = pd.concat(
            [pd.DataFrame(d, index=i, columns=fdv_columns) for i, d, c in results], copy=False
        )
    else:
        flowveldepth = pd.concat(
            [pd.DataFrame(d, index=i, columns=fdv_columns) for i, d in results], copy=False
        )
    flowveldepth = flowveldepth.sort_index()

    flows = flowveldepth.loc[:, (slice(None), "q")]
    flows = flows.T.reset_index(level=[0, 1])
    flows.rename(columns={"level_0": "Timestep", "level_1": "Parameter"}, inplace=True)
    flows["Time (d)"] = ((flows.Timestep + 1) * dt) / (24 * 60 * 60)
    flows = flows.set_index("Time (d)")

    # specify simulation calendar time
    dt_routing = pd.Timedelta(str(dt) + "seconds")
    time_routing = pd.period_range(
        (wrf_time[0] - dt_wrf + dt_routing), periods=nts, freq=dt_routing
    ).astype("datetime64[ns]")
    time_wrf = validation_data.columns.values

    # construct comparable dataframes
    trt = pd.DataFrame(
        flows.loc[:, compare_node].values,
        index=time_routing,
        columns=["flow, t-route (cms)"],
    )
    wrf = pd.DataFrame(
        validation_data.loc[compare_node, :].values,
        index=time_wrf,
        columns=["flow, wrf (cms)"],
    )

    # compare dataframes
    compare = pd.concat([wrf, trt], axis=1, sort=False, join="inner")
    print(compare)
