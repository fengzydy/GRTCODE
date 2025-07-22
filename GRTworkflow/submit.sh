#!/bin/bash

# Path to run script
runscript="$(pwd)/GRTworkflow/run-era5.sh"

# -----------------------------
# Available experiment options:
# These correspond to different GHG scenarios, each with a matching NetCDF file: $ghg_path/${exp}.nc
#
#   • Single-gas PI/4x/PD/2x/3x scenarios:
#     co2_PI, co2_2xPI, co2_3xPI, co2_4xPI, co2_PD, co2_4x
#     ch4_PI, ch4_2xPI, ch4_3xPI, ch4_4xPI, ch4_PD, ch4_4x
#     n2o_PI, n2o_2xPI, n2o_3xPI, n2o_4xPI, n2o_PD, n2o_4x
#     hfc134aeq_PI, hfc134aeq_2xPI, hfc134aeq_3xPI, hfc134aeq_4xPI, hfc134aeq_PD, hfc134aeq_4x
#     cfc12eq_PI, cfc12eq_2xPI, cfc12eq_3xPI, cfc12eq_4xPI, cfc12eq_PD, cfc12eq_4x
#
#   • Control experiments, time-varying well-mixed ghg following CMIP7:
#     control
#
#   • Reference experiment, holding well-mixed ghg at 1850:
#     PI
i
# To test or run one or more experiments, set the 'exps' array below:
# -----------------------------
exps=(PI)

# Directory containing greenhouse gas input files
ghg_path="/gpfs/f5/gfdl_m/world-shared/Jing.Feng/GHG"

# ERA5 input directory (options: era5_coarse, era5_coarse_noh2o, era5_fo3, era5_fo3strat)
era5_data="/gpfs/f5/gfdl_m/scratch/Jing.Feng/line-by-line/run/era5_noh2o"

# Loop over years
for year in $(seq 2010 2010); do
  # Loop over selected experiments
  for exp in "${exps[@]}"; do
    name_tag="${exp}_noh2o"  # Unique name tag for output

    # Construct and run command
    cmd="$runscript -p /ncrc/home1/Jing.Feng/scripts/grtcode $year $era5_data $ghg_path/${exp}.nc $name_tag"
    echo "Running: $cmd"
    eval "$cmd"
  done
done

