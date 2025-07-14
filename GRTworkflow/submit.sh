#!/bin/bash


runscript="/ncrc/home1/Jing.Feng/scripts/GRTworkflow/run-era5.sh"

# experiment list
exps=(control) 

# path to greenhouse gas input files.
# the filename of each input file is $data_dir/${exp}.nc
data_dir="/gpfs/f5/gfdl_m/world-shared/Jing.Feng/GHG"

# loop for years
for year in $(seq 2010 2010); do
    # loop for experiment names
    for exp in "${exps[@]}"; do
	# generate command line for job submission.
	# runscript is 
        cmd="$runscript -p /ncrc/home1/Jing.Feng/scripts/grtcode $year $data_dir/${exp}.nc ${exp}"
        echo "Running: $cmd"
        eval "$cmd"
    done
done
