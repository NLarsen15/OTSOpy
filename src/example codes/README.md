# OTSO Example Codes

This directory contains example usage scripts for the OTSO package demonstrating the new parameter grouping features.

## New Parameter Grouping System

OTSO now uses parameter groups to organize related parameters, making function calls cleaner and more maintainable. Instead of passing 50+ individual parameters, you can now group related parameters together:

### Basic Usage Examples

**Traditional approach (still supported):**
```python
cutoff("OULU", corenum=1, year=2000, month=1, vx=-400, by=3.0)
```

**New grouped approach:**
```python
cutoff("OULU", 
       computation_params={"corenum": 1},
       datetime_params={"year": 2000, "month": 1},
       solar_wind={"vx": -400, "by": 3.0})
```

### Available Parameter Groups

- **solar_wind**: `vx, vy, vz, bx, by, bz, by_avg, bz_avg, density, pdyn`
- **geomagnetic**: `Dst, kp, n_index, b_index, sym_h_corrected`
- **tsyganenko**: `G1, G2, G3, W1, W2, W3, W4, W5, W6`
- **datetime_params**: `year, month, day, hour, minute, second`
- **magfield_params**: `internalmag, externalmag, boberg, magnetopause, spheresize`
- **integration_params**: `intmodel, gyropercent, minaltitude, maxdistance, maxtime`
- **particle_params**: `Anum, anti`
- **rigidity_params**: `startrigidity, endrigidity, rigiditystep, rigidityscan`
- **coordinate_params**: `coordsystem, inputcoord`
- **computation_params**: `corenum, Verbose`
- **data_retrieval_params**: `serverdata, livedata`
- **custom_field_params**: `g, h, MHDfile, MHDcoordsys`

## Running Tests

### Local Testing

Install pytest if you haven't already:
```bash
pip install pytest pytest-timeout
```

**No OTSO installation required!** The tests automatically add the project root to Python path.

Run all example tests:
```bash
# From project root
pytest "src/example codes" -v

# Or from this directory
cd "src/example codes"
pytest test_examples.py -v
```

Run only quick tests (trajectory + cone):
```bash
pytest "src/example codes" -v -k "trajectory or cone"
```

Run a specific example test:
```bash
pytest "src/example codes/test_examples.py::TestOTSOExamples::test_example_runs_successfully[trajectory]" -v
```

### Test Options

- `-v` : Verbose output showing each test
- `-s` : Show print statements from tests  
- `-x` : Stop on first failure
- `--tb=short` : Shorter traceback format
- `-k trajectory` : Run tests matching "trajectory"

## Available Tests

1. **Individual Example Tests** - Each `.py` file is tested automatically
2. **Import Test** - Verifies OTSO package imports work

## Adding New Examples

1. Create a new `.py` file in this directory
2. Follow existing example patterns
3. The test will automatically discover and test it
4. No additional test code needed!

## Troubleshooting

**Import errors:**
The tests automatically handle Python path setup, so no package installation is required. If you still see import errors, ensure you're running from the project root or the example codes directory.

**Timeout errors:**
```bash
# Run with longer timeout
pytest --timeout=600 "src/example codes"
```

**See detailed output:**
```bash
pytest "src/example codes" -v -s --tb=long
```