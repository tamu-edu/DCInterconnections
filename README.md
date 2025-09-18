# Flexibility-Aware, Planner-Initiated Siting to Fast-Track Data Center Interconnections

This repository provides a **robust, resume-able MATLAB pipeline** to evaluate **1-GW data-center deployments** at candidate buses under **standardized flexibility envelopes** (Firm, Pause, Shift). The workflow follows a **planner-initiated, reliability-gated siting** approach: candidate nodes are screened for **N–1 feasibility** (can be extended to N-k reliability according to user needs), then **system-wide market impacts** are quantified via day-ahead SCUC and intra-day SCED, producing per-bus metrics for transparent shortlisting.

> **Highlights (research-aligned)**
> - **Planner-initiated siting.** Pre-defined operating envelopes and a vetted shortlist **before queue entry** to enable fast-track interconnection.
> - **Reliability-gated screening.** Ex-ante N–1 feasibility with configurable peak/off-peak windows.
> - **Market-impact quantification.** DA-SCUC + ID-SCED to summarize cross-node **LMP level/dispersion**, **binding-hour**, and **congestion rent**, with peak/off-peak decomposition.
> - **Scenario support.** Firm, Pause, and Shift envelopes for **flexibility-aware** comparisons at **1 GW** (easily adapted to other sizes).
> - **HPC-friendly reproducibility.** Atomic **day-level checkpoints**, per-bus folders, and **resume / aggregation-only** modes for crash-tolerant large runs.

---

## Repository layout

<pre>
DCInterconnections/
├─ README.md
├─ LICENSE
├─ matlab/
│  ├─ Flex_DC_Publish.m
│  └─ Reliability_screen.m
├─ scripts/
│  ├─ check_dependencies.m
│  └─ hprc_submit_example.slurm
└─ data/
   └─ README.md
</pre>

## Requirements

- **MATLAB** R2021a or later  
- **MATPOWER** & **MOST**  
- **YALMIP**  
- **GUROBI** (MATLAB interface on path, e.g., `gurobi.mexa64`)  
- Synthetic **Texas 2000-bus** case files on MATLAB path:
  - `YOUR MATPOWER Case File` (e.g., `DC_case2000`)
  - `YOUR Generator cost / commitment data` (e.g., `DC_GD_case2000`)
- Input time series (all `8760 × …` unless noted):
  - `YOUR_demand_profiles.mat` → variable `area_load` (MW, `8760 × nBus`)
  - `YOUR_wind_profiles.mat`   → variable `wind_MW`  (MW, `8760 × nWind`)
  - `YOUR_solar_profiles.mat`  → variable `solar_MW` (MW, `8760 × nSolar`)
  - `YOUR_hydro_profiles.mat`  → variable `hydro_MW` (MW, `8760 × nHydro`)
  - `YOUR_reliability.mat`     → vector of bus IDs that **pass the reliability gate**

> ⚠️ **Paths:** Replace `YOUR_…` placeholders in the MATLAB script and point to locations you control.

## Quick start

1. Clone the repo and open MATLAB at the repo root.
   ```matlab
   addpath(genpath(pwd));
2. Place required `.mat` inputs into `data/` or add their folder to the MATLAB path
3. (Optional) Edit in `matlab/Flex_DC_Publish.m`:
   - `start_day`, `end_day` (DOY integers)
   - `flexMode = 'firm'`, `sizeDC = 1000` (MW)
   - `peakHours = 16:19`
   - `headroom = 0.8`, `pausePCT = 0.85`, `shiftPCT = 0.80`
4. Run:
   ```matlab
   cd matlab
   run('Flex_DC_Publish.m');

## Outputs
- Per candidate bus **b**:
  ```bash
  <result_root>/bus####_firm_1000MW/
  scuc_sced/24h-SCUC-day-<DOY>.mat     # λ (LMP), μ (constraint duals), flows, Fmax
  metrics_summary.mat                  # per-bus price & congestion metrics
  DONE.txt
  ```
  
- Consolidated file:
  ```bash
  <result_root>/all_metrics_base_1000MW.mat
  - priceMetrics_struct : all-hour / peak / off-peak LMP stats
  - congMetrics_struct  : binding-hour & congestion rent
  ```

## Running on SLURM
- Edit and submit:
  ```bash
  sbatch scripts/hprc_submit_example.slurm
  ```
