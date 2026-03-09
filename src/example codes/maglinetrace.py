from OTSO import trace

if __name__ == '__main__':

    # Example using grouped parameters
    trace_results = trace(
        computation_params={"corenum": 4},
        magfield_params={"externalmag":"NONE"},
        grid_params={"latstep": -10, "lonstep": 30, "maxlat": 60, "minlat": 0},
    )

    print(trace_results[0]) # dictionary output containing positional information magnetic field lines generated over
                    # the globe
    print(trace_results[1]) # text output of input variable information

