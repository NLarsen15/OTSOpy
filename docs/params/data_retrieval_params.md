# Data retrieval

Dictionary: `data_retrieval_params` · Class: `DataRetrievalParams`

**Used by:** all functions

Automatic space-weather data retrieval, as an alternative to manually filling in `solar_wind_params`/`geomagnetic_params`/`tsyganenko_params`.

| Key | Type | Default | Description |
|---|---|---|---|
| `serverdata` | str | `"OFF"` | `"ON"` downloads and uses the finalised OMNI database values for the given `datetime_params`, overriding manually-supplied solar wind/geomagnetic/Tsyganenko values. |
| `livedata` | str | `"OFF"` | `"ON"` pulls preliminary real-time data from NOAA instead. NOAA values are provisional and may later differ from the finalised OMNI database (`serverdata` takes precedence if both would apply). Not supported for `externalmag="TSY04"` or `"TA16_RBF"`. |
