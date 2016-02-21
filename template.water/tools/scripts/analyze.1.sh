#!/bin/bash

source parameters.sh

echo "1. $run_dt"
echo "2. $run_nstep"
echo "3. $run_corr_dt"
echo "4. $analyze_mode_step"
echo "5. $analyze_corr_end_t"

analyze_corr_end_n=`echo "($analyze_corr_end_t + 0.5 * $run_corr_dt) / $run_corr_dt" | bc`
echo "analyze_corr_end_n=  $analyze_corr_end_n"

command="./1/velo.mode.corr \
-c ../2016Jan17out/conf.gro \
-f ../2016Jan17out/traj.trr \
--mode-start $1 \
--mode-step $analyze_mode_step \
--mode-end $2 \
--n-corr $analyze_corr_end_n"

echo "# command is $command"
echo "# command is $command" &> analyze.$1.$2.log
echo "# begin analyze"
$command &> analyze.$1.$2.log

