#!/bin/bash

car_file=car.gro
sol_file=sol.gro

car_low=0.045
car_upp=2.055
sol_upp=`tail $sol_file -n 1 | awk '{print $3}'`
sol_x=`tail $sol_file -n 1 | awk '{print $1}'`
sol_y=`tail $sol_file -n 1 | awk '{print $2}'`

gap=0.22
vacu=3

