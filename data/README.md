# Data folder

This folder holds the time-series inputs and small helper files required by the MATLAB pipeline.

---

## Expected `.mat` files

| Filename                      | Variable name(s)          | Shape (rows × cols)         | Units | Description / Notes |
|---                            |---                        |---                           |---    |---|
| `texas_2020_demand.mat`      | `area_load`               | `8760 × nBus`                | MW    | Hourly total demand **per bus** for the study year (UTC or local—see “Time conventions”). Each column maps to a bus index in the MATPOWER case. |
| `texas_2020_wind.mat`        | `wind_MW`                 | `8760 × nWind`               | MW    | Hourly available wind generation **per wind unit** (same order as `gen` rows tagged as wind in the case/XGD). Values typically represent PMAX profiles. |
| `texas_2020_solar.mat`       | `solar_MW`                | `8760 × nSolar`              | MW    | Hourly available solar generation **per solar unit** (PMAX profiles). |
| `texas_2020_hydro.mat`       | `hydro_MW`                | `8760 × nHydro`              | MW    | Hourly available hydro generation **per hydro unit** (PMAX profiles). |
| `texas_2020_generation.mat`  | `min_on`, `min_off`       | `nGen × 1`                   | hours | **Per-generator minimum on/off times**. Order **must** match the generator order in your MATPOWER case/XGD. |
| `reliability.mat`            | `reliability`             | `nBus × 1`                   | –     | Reliability-gate output, a vector of **qualified bus IDs**. |

### Time conventions & indexing
- **Indexing:** rows are hours `1…8760` for a non-leap year.
- **Time zone:** specify whether series are **local** or **UTC**. The default code treats `peakHours = 16:19` as **local**; if your data are UTC, adjust accordingly.
- **Missing values:** avoid `NaN`. If present, fill or document the imputation rule (forward-fill, zeros for availability PMAX, etc.).

---

## Quick shape check (MATLAB)

```matlab
% Hourly profiles
D = load('texas_2020_demand.mat');  assert(ismatrix(D.area_load) && size(D.area_load,1)==8760);
W = load('texas_2020_wind.mat');    assert(size(W.wind_MW,1)==8760);
S = load('texas_2020_solar.mat');   assert(size(S.solar_MW,1)==8760);
H = load('texas_2020_hydro.mat');   assert(size(H.hydro_MW,1)==8760);

% Per-generator constraints (static)
G = load('texas_2020_generation.mat'); 

% Reliability gate
R = load('reliability.mat');  % A vector of bus IDs
