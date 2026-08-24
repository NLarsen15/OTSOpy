import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

from OTSO import cone


# ============================================================
# OTSO INPUT PARAMETERS
# ============================================================

date_inputs = {
    "year": 2020,
    "month": 1,
    "day": 1,
    "hour": 1,
    "minute": 0,
    "second": 0
}

magfield_inputs = {
    "internalmag": "IGRF",
    "externalmag": "TSY89c",
    "boberg": False,
    "bobergtype": "EXTENSION",
    "magnetopause": "Sphere",
    "spheresize": 10,
    "AdaptiveExternalModel": False
}

rigidity_inputs = {
    "startrigidity": 20,
    "endrigidity": 0,
    "rigiditystep": 0.01
}

transmission_inputs = {
    "transmission": False,
    "transmissionRstep": 0.001,
    "transmissionsamples": 20
}

solar_wind_inputs = {
    "vx": -500,
    "vy": 0,
    "vz": 0,
    "bx": 0,
    "by": 1,
    "bz": 1,
    "by_avg": 0,
    "bz_avg": 0,
    "density": 1,
    "pdyn": 0
}

geomagnetic_inputs = {
    "Dst": 0,
    "kp": 2,
    "n_index": 0,
    "b_index": 0,
    "sym_h_corrected": 0
}

tsyganenko_inputs = {
    "G1": 0,
    "G2": 0,
    "G3": 0,
    "W1": 0,
    "W2": 0,
    "W3": 0,
    "W4": 0,
    "W5": 0,
    "W6": 0
}

integration_inputs = {
    "intmodel": "Boris-Buneman",
    "gyropercent": 1,
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "mintrapdist": 6.6,
    "startaltitude": 20,
    "betaerror": 0.001,
    "totalbetacheck": True,
    "adaptivestep": True,
    "maxsteps": 100000
}

particle_inputs = {
    "Anum": 1,
    "anti": "YES",
    "zenith": 0,
    "azimuth": 0
}

computation_inputs = {
    "corenum": 1,
    "threadnum": 1,
    "Verbose": True,
    "delim": ";"
}

custom_field_inputs = {
    "g": None,
    "h": None,
    "max_degree": 13,
    "MHDfile": None,
    "MHDcoordsys": None
}

coord_inputs = {
    "coordsystem": "GEO",
    "inputcoord": "GDZ"
}

data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "OFF"
}


# ============================================================
# MAIN PROGRAM
# ============================================================

if __name__ == "__main__":

    # --------------------------------------------------------
    # Stations
    # --------------------------------------------------------

    stations_list = [
        "OULU",
        "CALG"
    ]


    # ========================================================
    # RUN OTSO CONE CALCULATION
    # ========================================================

    cone_results = cone(
        Stations=stations_list,
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params=rigidity_inputs,
        transmission_params=transmission_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params=integration_inputs,
        particle_params=particle_inputs,
        computation_params=computation_inputs,
        coordinate_params=coord_inputs,
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs
    )


    # ========================================================
    # EXTRACT RESULTS
    # ========================================================

    cone_df = cone_results[0]
    cutoff_df = cone_results[1]


    # ========================================================
    # PRINT RESULTS
    # ========================================================

    print()
    print("========== CONE RESULTS ==========")
    print(cone_df)

    print()
    print("========== CUTOFF RESULTS ==========")
    print(cutoff_df)


    # ========================================================
    # CREATE ROBINSON MAP
    # ========================================================

    fig = plt.figure(
        figsize=(14, 7)
    )

    ax = plt.axes(
        projection=ccrs.Robinson()
    )

    ax.set_global()


    # --------------------------------------------------------
    # Coastlines
    # --------------------------------------------------------

    ax.coastlines(
        linewidth=0.8
    )


    # --------------------------------------------------------
    # Country borders
    # --------------------------------------------------------

    ax.add_feature(
        cfeature.BORDERS,
        linewidth=0.3
    )


    # ========================================================
    # GRIDLINES
    # ========================================================

    # IMPORTANT:
    #
    # draw_labels=True caused the Shapely error:
    #
    # GEOSException:
    # Points of LinearRing do not form a closed linestring
    #
    # Therefore labels are disabled here.

    gl = ax.gridlines(
        draw_labels=False,
        color="gray",
        linestyle="--",
        linewidth=0.5
    )

    # Gridline spacing
    gl.xlocator = mticker.FixedLocator(
        range(-180, 181, 30)
    )

    gl.ylocator = mticker.FixedLocator(
        range(-90, 91, 30)
    )


    # ========================================================
    # COLOR SCALE
    # ========================================================

    norm = plt.Normalize(
        vmin=0,
        vmax=20
    )

    cmap = plt.cm.viridis


    # ========================================================
    # PLOT ASYMPTOTIC CONES
    # ========================================================

    sc = None


    # --------------------------------------------------------
    # Loop through all stations
    # --------------------------------------------------------

    for station in cone_df.columns[1:]:

        lats = []
        lons = []
        rigidities = []


        # ----------------------------------------------------
        # Loop through rigidity values
        # ----------------------------------------------------

        for _, row in cone_df.iterrows():

            # Rigidity
            rigidity = float(
                row["R [GV]"]
            )

            # OTSO result format:
            #
            # filter;latitude;longitude
            #
            # Example:
            #
            # 1;5.7192;281.4538

            value = row[station]


            try:

                filter_str, lat_str, lon_str = str(
                    value
                ).split(";")


                # --------------------------------------------
                # Filter value
                # --------------------------------------------

                filter_value = int(
                    filter_str
                )


                # --------------------------------------------
                # Only plot valid asymptotic directions
                # --------------------------------------------

                if filter_value != 1:
                    continue


                # --------------------------------------------
                # Latitude
                # --------------------------------------------

                lat = float(
                    lat_str
                )


                # --------------------------------------------
                # Longitude
                #
                # OTSO gives longitude in 0-360 degrees.
                # Convert to -180 to +180.
                # --------------------------------------------

                lon = float(
                    lon_str
                )

                if lon > 180:
                    lon -= 360


                # --------------------------------------------
                # Store data
                # --------------------------------------------

                lats.append(lat)
                lons.append(lon)
                rigidities.append(rigidity)


            except (
                ValueError,
                AttributeError
            ):
                continue


        # ====================================================
        # PLOT THIS STATION
        # ====================================================

        if len(lats) > 0:

            sc = ax.scatter(
                lons,
                lats,
                c=rigidities,
                cmap=cmap,
                norm=norm,
                s=10,
                transform=ccrs.PlateCarree(),
                alpha=0.75,
                label=station
            )


            print(
                f"{station}: "
                f"{len(lats)} valid asymptotic points"
            )

        else:

            print(
                f"{station}: "
                f"no valid asymptotic points"
            )


    # ========================================================
    # COLORBAR
    # ========================================================

    if sc is not None:

        cbar = plt.colorbar(
            sc,
            ax=ax,
            orientation="horizontal",
            fraction=0.08,
            pad=0.08
        )

        cbar.set_label(
            "Rigidity [GV]",
            fontsize=14
        )

        cbar.ax.tick_params(
            labelsize=12
        )


    # ========================================================
    # TITLE
    # ========================================================

    ax.set_title(
        "Asymptotic Cones\n",
        fontsize=15,
        pad=15
    )


    # ========================================================
    # SAVE FIGURE
    # ========================================================

    plt.tight_layout()

    plt.savefig(
        "asymptotic_cones.png",
        dpi=300,
        bbox_inches="tight"
    )


    # ========================================================
    # DISPLAY
    # ========================================================

    plt.show()