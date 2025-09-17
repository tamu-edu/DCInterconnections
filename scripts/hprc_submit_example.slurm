#!/bin/bash
#SBATCH -J flex_dc_firm
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p processes
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x-%j.out
#SBATCH -e logs/%x-%j.err

# Load environment (EDIT for your cluster)
module load matlab/R2021a || true
module load gurobi || true

mkdir -p logs

# Run MATLAB non-interactively
MATLAB_SCRIPT="matlab/Flex_DC_Publish.m"

matlab -nodisplay -nosplash -r "try, addpath(genpath(pwd)); run('$MATLAB_SCRIPT'); catch ME, disp(getReport(ME)); exit(1); end; exit(0);"
