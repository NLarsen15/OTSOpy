# Data retrieval

Dictionary: `data_retrieval_params` · Class: `DataRetrievalParams`

**Used by:** all functions

Automatic space-weather data retrieval, as an alternative to manually filling in `solar_wind_params`/`geomagnetic_params`/`tsyganenko_params`.

| Key | Type | Default | Description |
|---|---|---|---|
| `serverdata` | str | `"OFF"` | `"ON"` downloads and uses the finalised OMNI database values for the given `datetime_params`, overriding manually-supplied solar wind/geomagnetic/Tsyganenko values. |
| `livedata` | str | `"OFF"` | `"ON"` pulls preliminary real-time data from NOAA instead. NOAA values are provisional and may later differ from the finalised OMNI database (`serverdata` takes precedence if both would apply). Not supported for `externalmag="TSY04"` or `"TA16_RBF"`. |

## Server data files and the `Source` column

With `serverdata="ON"`, OTSO builds one file per year, `{year}_TSY_Inputs.csv`, at 5-minute resolution (hourly before 1981). Each row has a `Source` column that records where that row's solar wind and IMF values came from. The derived Tsyganenko inputs (G1–G3, W1–W6, N-index and B-index) are computed from those values, so they share the row's source. Kp and Dst always come from hourly OMNI, and SYM-H from 5-minute OMNI.

When more than one applies, a row gets the highest code.

| `Source` | Meaning |
|---|---|
| 1 | OMNI 5-minute data. |
| 2 | OMNI hourly data: no valid 5-minute record exists, so the hourly value is repeated across the hour. |
| 3 | OMNI hourly data with the solar wind speed and/or density taken from SOHO CELIAS (from 1996), and Pdyn recomputed. |
| 4 | At least one solar wind or IMF value was linearly interpolated across a data gap of 3 hours or less. |
| 5 | At least one solar wind or IMF value is missing (data gap longer than 3 hours). |
| 6 | Before 1981: hourly OMNI Kp and Dst only. No solar wind or IMF data exist, and those values are set to 0. |

Codes 2–6 mark inputs that are not direct 5-minute OMNI measurements, which may be worth stating in publications. Files downloaded with OTSO versions before 1.3.9 have no `Source` column; delete them (`OTSO.clean`) to rebuild them with it.

G1–G3 are means over the preceding hour of the final 5-minute series: G1 = ⟨V h(B⊥) sin³(θ/2)⟩, G2 = ⟨V Bs⟩/200 (Tsyganenko 2002b) and G3 = ⟨N V Bs⟩/2000 (Tsyganenko et al. 2003).
