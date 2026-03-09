![Logo](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/OTSO_logo.png)

# OTSOpy
Python package version of the OTSO tool used for trajectory computations of charged particles in the Earth's magnetosphere.

OTSO is designed to be open-source; all suggestions for improvement are welcome, and please report any bugs you find. I welcome any help provided by the community in the development of OTSO.

__Supported Python Versions:__ 3.10, 3.11, 3.12, 3.13, and 3.14 
(I will endeavour to keep OTSO support as up to date as possible)

# Installation

Installation of OTSOpy is designed to be as simple as possible and can be done utilising pip. Users have two options when downloading OTSOpy.

## Option 1: PyPi
Users may install OTSO directly from PyPi using:

`pip install OTSO` 

This will install OTSO into your current Python environment.

## Option 2: Repository
Users may clone the repository and run the setup.py file within the main OTSOpy directory using:

`pip install .`

This will install OTSO into your current Python environment.

# Troubleshooting

## LINUX
Sometimes there are errors regarding libgfortran. Make a note of the libgfortran error message and then install the appropriate libgfortran version that is being requested. This should resolve the issue.

## MAC
The compiled fortran libraries can be flagged as potential malware. To resolve this you can attempt to compile the libraries yourself or in your settings grant permission for your computer to access the required .so file. 

# Functions

## Cutoff
Computes the geomagnetic cut-off rigidities for given locations around the Earth under user-inputted geomagnetic conditions.

![Cutoff](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/cutoffplot.png)
*Figure 1: Computation of the Oulu neutron monitor effective cut-off rigidity using the IGRF 2000 epoch and TSY89 model with kp index = 0. Penumbra is shown by the forbidden and allowed trajectories being black and white, respectively. The upper and lower cut-off values (Ru and Rl) are denoted in the legend, from which the effective cut-off (Rc) is computed.*

## Cone
Computes the asymptotic viewing directions for given locations around the Earth. Asymptotic latitudes and longitudes over a range of rigidity values are computed.
Asymptotic latitude and longitude can be given in any available coordinate system.

![Cones](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/coneplot.png)

*Figure 2: Asymptotic cones for the Oulu, Nain, South Pole, Thule, and Inuvik neutron monitors for the IGRF 2010 epoch and TSY89 model, with kp = 0. Latitudes and longitudes are in the geocentric coordinate system.*

## Trajectory
Computes and outputs the trajectory of a charged particle with a specified rigidity from a given start location on Earth. Positional information can be in any of the available coordinate systems.

![Trajectory](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/Trajectory_Plot.png)

*Figure 3: Computed trajectories of three cosmic rays of various rigidity values being backtraced from the Oulu neutron monitor for the IGRF 2000 and TSY89 model, with kp = 0. The 1GV particle is allowed (able to escape the magnetosphere); the 0.4GV particle is forbidden (it is trapped in the magnetosphere); and the 0.1GV is also forbidden (it returns to Earth).*

## Planet
Performs the cutoff function over a user-defined location grid, allowing for cutoffs for the entire globe to be computed instead of individual locations. There is the option to return the asymptotic viewing directions at each computed location by utilising a user-inputted list of rigidity levels.


![Planet](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/planetplot.png)
*Figure 4: Computed vertical effective cut-off rigidities across a 5°x5° grid of the Earth. These computations were done using the IGRF 2000 epoch and TSY89 model, with kp = 0.*

## Flight
Computes the cut-off rigidities along a user-defined path. The function is named Flight as it is primarily been developed for use in aviation tools, but any path can be entered. For example, the function can be applied to geomagnetic latitude surveys using positional data from a ship voyage, or it can be used to compute anisotropy and cut-off values for low-Earth orbit spacecraft. This function allows for changing altitude, location, and date values. 

![ISS](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/ISS_cutoffs.png)

*Figure 3: Computed effective vertical cut-off rigidities for the ISS between the 15th and 16th of March 2021. Geomagnetic parameters were extracted directly from OMNI for this period.*

## Trace
Traces the magnetic field lines around the globe or for a given location based on the geomagnetic configuration detailed by the user. It is useful for modelling the magnetosphere structure under disturbed conditions and for finding open magnetic field lines.

![Trace](https://raw.githubusercontent.com/NLarsen15/OTSOpy/main/src/images/traceplot.png)
*Figure 5: Computation of magnetic field line configuration in the X-Z plane on January 1st 2000 12:00:00. IGRF and TSY01 models used, and input variables were obtained using the server data option within OTSO.*

## Coordtrans
Converts input positional information from one coordinate system to another, utilising the [IRBEM](https://github.com/PRBEM/IRBEM) library of coordinate transforms.

## Magfield
Computes the total magnetic field strength at a given location depending on the user's input geomagnetic conditions. Outputs will be in the geocentric solar magnetospheric (GSM) coordinate system.

# Examples

## Cutoff

```python
import OTSO

if __name__ == '__main__':
    stations_list = ["OULU", "ROME", "ATHN", "CALG"]  # list of neutron monitor stations (using their abbreviations)

    cutoff = OTSO.cutoff(
        Stations=stations_list,
        computation_params={"corenum": 1},
        datetime_params={"year": 2000, "month": 1, "day": 1, "hour": 0}
    )

    print(cutoff[0])  # dataframe output containing Ru, Rc, Rl for all input locations
    print(cutoff[1])  # text output of input variable information
```

### Output
Ru = upper cut-off rigidity [GV]

Rc = effective cut-off rigidity [GV]

Rl = lower cut-off rigidity [GV]

```
    ATHN  CALG  OULU  ROME
Ru  8.98  1.14  0.72  6.35
Rc  8.33  1.07  0.71  6.07
Rl  6.31  1.02  0.59  5.37
```

## Cone

```python
import OTSO

if __name__ == '__main__':

    stations_list = ["OULU", "ROME", "ATHN", "CALG"]  # list of neutron monitor stations (using their abbreviations)

    cone = OTSO.cone(
        Stations=stations_list,
        computation_params={"corenum": 1},
        datetime_params={"year": 2000, "month": 1, "day": 1, "hour": 0}
    )

    print(cone[0])  # dataframe output containing asymptotic cones for all input locations
    print(cone[1])  # dataframe output containing Ru, Rc, Rl for all inputted locations
    print(cone[2])  # text output of input variable information
```

### Output
Showing only the cone[0] output containing the asymptotic viewing directions of the input stations. Result layout is: filter;latitude;longitude.
If the filter value is 1, then the particle of that rigidity has an allowed trajectory. If the filter value is NOT 1, then the particle of that rigidity has a forbidden trajectory.
```
      R [GV]                ATHN               CALG               OULU               ROME
0     20.000     1;-1.635;89.172   1;21.147;271.975    1;40.902;62.437     1;4.052;71.083
1     19.990     1;-1.661;89.200   1;21.131;271.973    1;40.890;62.435     1;4.027;71.101
2     19.980     1;-1.687;89.228   1;21.117;271.972    1;40.877;62.434     1;4.001;71.120
3     19.970     1;-1.713;89.256   1;21.101;271.970    1;40.865;62.432     1;3.975;71.138
4     19.960     1;-1.739;89.283   1;21.086;271.969    1;40.853;62.431     1;3.950;71.156
...      ...                 ...                ...                ...                ...
1995   0.050    -1;5.784;175.702   -1;38.217;16.098  -1;43.597;238.366   -1;3.370;219.840
1996   0.040   -1;20.165;207.726   -1;20.122;39.597  -1;24.492;229.951  -1;24.971;178.788
1997   0.030  -1;-32.777;224.618   -1;29.968;12.264  -1;17.566;214.996  -1;17.415;216.472
1998   0.020    -1;7.967;228.903   -1;33.543;96.634  -1;36.906;180.237  -1;30.685;186.313
1999   0.010   -1;19.485;224.643  -1;-4.295;338.997  -1;57.813;160.815  -1;26.726;219.327
```

## Trajectory

```python
import OTSO

if __name__ == '__main__':

    stations_list = ["OULU", "ROME", "ATHN", "CALG"]  # list of neutron monitor stations (using their abbreviations)

    trajectory = OTSO.trajectory(
        Stations=stations_list,
        rigidity_params={"startrigidity": 5},
        computation_params={"corenum": 1}
    )

    print(trajectory[0])  # dictionary output containing positional information for all trajectories generated starting
                         # from input stations
    print(trajectory[1])  # text output of input variable information

```

### Output
Showing the dataframe produced for the particle originating from Oulu. Other trajectories are within the trajectory[0] dictionary. Additionally the Filter value, letting you know if the trajectory is allowed or not, and the asymptotic latitude and longitude at the end point is included. 

```
{'NMname': 'OULU', 'trajectory':

       X_Re [GEO]  Y_Re [GEO]  Z_Re [GEO]
0      0.383531    0.182689    0.907014
1      0.383572    0.182709    0.907112
2      0.383618    0.182731    0.907221
3      0.383667    0.182756    0.907340
4      0.383722    0.182784    0.907471
..          ...         ...         ...
445  -28.578000   21.122000   -1.970810
446  -31.623400   21.590600   -1.878430
447  -34.974700   22.101100   -1.794190
448  -38.660200   22.672200   -1.729170
449  -42.710200   23.328800   -1.698580

[450 rows x 3 columns], 
'Filter': 1, 'AsymLat': 0.103, 'AsymLong': 170.53}
````

## Planet

```python
import OTSO

if __name__ == '__main__':

    planet = OTSO.planet(
        cutoff_comp="Vertical",
        computation_params={"corenum": 1},
        datetime_params={"year": 2000},
        rigidity_params={"rigiditystep": 0.1}
    )

    print(planet[0])  # dataframe containing cutoff results for planet grid
    print(planet[1])  # text output of input variable information
```

### Output
The default output is a 5°x5° grid of the Earth with no asymptotic viewing directions computed.

```
      Latitude  Longitude   Ru   Rc   Rl
0        -90.0        0.0  0.0  0.0  0.0
1        -90.0        5.0  0.0  0.0  0.0
2        -90.0       10.0  0.0  0.0  0.0
3        -90.0       15.0  0.0  0.0  0.0
4        -90.0       20.0  0.0  0.0  0.0
...        ...        ...  ...  ...  ...
1382      90.0      340.0  0.0  0.0  0.0
1383      90.0      345.0  0.0  0.0  0.0
1384      90.0      350.0  0.0  0.0  0.0
1385      90.0      355.0  0.0  0.0  0.0
1386      90.0      360.0  0.0  0.0  0.0
```

## Flight

```python
import OTSO
import datetime

if __name__ == '__main__':

    latitude_list = [10, 15, 20, 25, 30]  # [Latitudes]
    longitude_list = [10, 15, 20, 25, 30]  # [Longitudes]
    altitude_list = [30, 40, 50, 60, 80]  # [Altitudes] in km
    date_list = [datetime.datetime(2000, 10, 12, 8), datetime.datetime(2000, 10, 12, 9), datetime.datetime(2000, 10, 12, 10),
                 datetime.datetime(2000, 10, 12, 11), datetime.datetime(2000, 10, 12, 12)]  # [dates]

    flight = OTSO.flight(
        latitudes=latitude_list,
        longitudes=longitude_list,
        dates=date_list,
        altitudes=altitude_list,
        cutoff_comp="Vertical",
        computation_params={"corenum": 1}
    )
    
    print(flight[0])  # dataframe output containing Ru, Rc, Rl along flightpath
    print(flight[1])  # text output of input variable information
    print(flight[2])  # dataframe output of input variables
```

### Output
flight[0] dataframe output.

```
                  Date  Latitude  Longitude  Altitude     Ru     Rc     Rl
0  2000-10-12 08:00:00        10         10        30  14.86  14.86  14.86
1  2000-10-12 09:00:00        15         15        40  14.90  14.90  14.90
2  2000-10-12 10:00:00        20         20        50  14.48  14.48  14.48
3  2000-10-12 11:00:00        25         25        60  13.59  13.59  13.59
4  2000-10-12 12:00:00        30         30        80  12.18  11.49  10.39
```

## Trace

```python
import OTSO

if __name__ == '__main__':

    trace = OTSO.trace(
        computation_params={"corenum": 1},
        grid_params={"latstep": -10, "lonstep": 30}
    )

    print(trace[0])  # dictionary output containing positional information of magnetic field lines generated over
                    # the globe
    print(trace[1])  # text output of input variable information
```

### Output
Example output of one of the field line traces for the location latitude = 60° and longitude = 215°.
The L shell and Invariant Latitude are also computed from the magnetic field line tracing. 

```
'60_215': {'Trace':        
         X_GEO [Re]  Y_GEO [Re]  Z_GEO [Re]  Bx_GSM [nT]  By_GSM [nT]  Bz_GSM [nT]
0       -0.568632    0.010497   -0.820734    -0.000024    -0.000017    -0.000053
1       -0.569193    0.010380   -0.821270    -0.000024    -0.000017    -0.000053
2       -0.569754    0.010263   -0.821806    -0.000024    -0.000017    -0.000053
3       -0.570316    0.010146   -0.822342    -0.000024    -0.000017    -0.000053
4       -0.570877    0.010029   -0.822877    -0.000024    -0.000017    -0.000052
...           ...         ...         ...          ...          ...          ...
11862   -0.412834   -0.288995    0.866600     0.000047     0.000014    -0.000024
11863   -0.411355   -0.288183    0.864958     0.000047     0.000015    -0.000024
11864   -0.410863   -0.287913    0.864410     0.000047     0.000015    -0.000024
11865   -0.410370   -0.287643    0.863862     0.000047     0.000015    -0.000025
11866   -0.409878   -0.287372    0.863314     0.000047     0.000015    -0.000025

[11867 rows x 6 columns], 
'L_shell': 4.1267, 'Invariant_Latitude': 60.5105}
```

## Coordtrans

```python
import OTSO
import datetime

if __name__ == '__main__':

    lat_lon_alt_list = [[10, 10, 10]]  # [[Latitude,Longitude,Altitude]]
    date_list = [datetime.datetime(2000, 10, 12, 8)]  # [dates]
    
    Coords = OTSO.coordtrans(
        Locations=lat_lon_alt_list,
        dates=date_list,
        CoordIN="GEO",
        CoordOUT="GSM",
        corenum=1  # coordtrans uses individual parameters, not grouped ones
    )

    print(Coords[0])  # dataframe output of converted coordinates
    print(Coords[1])  # text output detailing the initial and final conversion coordinate system
```

### Output
Coords[0] output converting the [10,10,10] position from GEO coordinate system to GSM coordinate system. 

```
                  Date X_GEO [Re] Y_GEO [Re] Z_GEO [Re] X_GSM [Re] Y_GSM [Re] Z_GSM [Re]
0  2000-10-12 08:00:00    1.00157       10.0       10.0   7.508443   6.239794  10.280632
```

## Magfield

```python
import OTSO

if __name__ == '__main__':

    location_list = [[10, 10, 10]]  # [[X,Y,Z]] Earth radii Geocentric coordinates in this instance

    magfield = OTSO.magfield(
        Locations=location_list,
        coordinate_params={"coordsystem": "GEO"},
        computation_params={"corenum": 1}
    )

    print(magfield[0])  # dataframe of returned magnetic field vectors at input locations
    print(magfield[1])  # text output of input variable information

```

### Output
magfield[0] output showing the magnetic field vector at the input location in the GSM coordinate system. 

```
   X_GEO [Re]  Y_GEO [Re]  Z_GEO [Re]  GSM_Bx [nT]  GSM_By [nT]  GSM_Bz [nT]
0        10.0        10.0        10.0    10.735517    -2.413889    11.586277
```

# Acknowledgements
The fantastic IRBEM library has been used in the development of OTSO, which proved an invaluable asset and greatly sped up development. The latest release of the IRBEM library can be found at [https://doi.org/10.5281/zenodo.6867552](https://doi.org/10.5281/zenodo.6867552). Thank you to N. Tsyganenko for the development of the external magnetic field models and their code, which are used within OTSO.

Thank you to Don and Peggy Smart for their insightful discussion on the nature of cutoff computations and for providing me with a copy of their cutoff computation tool, from which I learnt a lot and adopted many of their inspired optimisation techniques.

A wider thanks goes to the space physics community who, through the use of the original [OTSO](https://github.com/NLarsen15/OTSO), provided invaluable feedback, advice on improvements, and bug reporting. All discussions and advice have aided in the continual development and improvement of OTSO, allowing it to fulfil its aim of being a community-driven open-source tool. The lessons learned from the initial OTSO versions have been incorporated into OTSOpy. Dr. Chris Davis was also instrumental in the development of OTSOpy with his suggestion of incorporating OTSO into the [AniMARIE](https://github.com/ssc-maire/AniMAIRE-public) tool, initiating the package development and providing help by expanding functionality and bug fixing. My personal thanks to Dr. Sergey Koldobsky for lending me his MacBook for MacOS Fortran compilations, expanding the number of available operating systems for OTSOpy.
OTSO was developed at the University of Oulu as part of the Academy of Finland QUASARE project. I would like to thank my colleagues at the University and the Academy of Finland for supporting the work.



# References
- **Larsen, N., Mishev, A., & Usoskin, I. (2023). A new open-source geomagnetosphere propagation tool (OTSO) and its applications. Journal of Geophysical Research: Space Physics, 128, e2022JA031061. https://doi.org/10.1029/2022JA031061**
