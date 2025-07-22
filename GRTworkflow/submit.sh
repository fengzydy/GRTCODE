#!/bin/bash


runscript="/ncrc/home1/Jing.Feng/scripts/grtcode/GRTworkflow/run-era5.sh"

# experiment list
exps=(PI) 
# full list of experiment name:
#  ls /gpfs/f5/gfdl_m/world-shared/Jing.Feng/GHG
# cfc12eq_4x	  cfc12eq_PI  ch4_4x    ch4_PI    co2_4x	co2_PI hfc134aeq_4x    hfc134aeq_PI  n2o_4x	  n2o_PI
# cfc12eq_2xPI		 cfc12eq_4xPI  ch4_2xPI	 ch4_4xPI  co2_2xPI  co2_4xPI	control   hfc134aeq_2xPI hfc134aeq_4xPI  n2o_2xPI      n2o_4xPI  PI
# cfc12eq_3xPI		 cfc12eq_PD	  ch4_3xPI	 ch4_PD    co2_3xPI  co2_PD	hfc134aeq_3xPI	hfc134aeq_PD    n2o_3xPI      n2o_PD
#
# path to greenhouse gas input files.
# the filename of each input file is $data_dir/${exp}.nc
#
data_dir="/gpfs/f5/gfdl_m/world-shared/Jing.Feng/GHG"

era5_data="/gpfs/f5/gfdl_m/scratch/Jing.Feng/line-by-line/run/era5_noh2o" 

# loop for years
for year in $(seq 2010 2010); do
    # loop for experiment names
    for exp in "${exps[@]}"; do
	# generate command line for job submission.
	#   <runscript>       <path of executable>            <year>   <era5 inputs>  <ghg inputs>   <experiment name>  
        cmd="$runscript -p /ncrc/home1/Jing.Feng/scripts/grtcode $year $era5_data $data_dir/${exp}.nc ${exp}_noh2o"
        echo "Running: $cmd"
        eval "$cmd"
    done
done
