#!/bin/bash

source parameters.sh

out_file=conf.gro

if test -f $out_file; then
    echo "file $out_file exists, do nothing"
    exit
fi

sol_trans=`echo "$car_upp + $gap" | bc -l`
car_trans=`echo "$sol_trans + $sol_upp + $gap - $car_low" | bc -l`
editconf -translate 0 0 $sol_trans -f $sol_file
mv -f out.gro tmp_sol.gro.tmp
editconf -translate 0 0 $car_trans -f $car_file
mv -f out.gro tmp_car.gro.tmp
new_box_z=`echo "$car_trans + $car_upp - $car_low + $vacu" | bc -l`

car_natom=`wc $car_file | awk '{print $1}'`
car_natom=`echo "$car_natom -3" | bc`

sol_natom=`wc $sol_file | awk '{print $1}'`
sol_natom=`echo "$sol_natom -3" | bc`

tot_natom=`echo "$sol_natom + $car_natom * 2" | bc`

echo ""					> $out_file
echo "$tot_natom"			>> $out_file
head -n $(($car_natom+2)) $car_file	  | tail -n $car_natom >> $out_file
head -n $(($sol_natom+2)) tmp_sol.gro.tmp | tail -n $sol_natom >> $out_file
head -n $(($car_natom+2)) tmp_car.gro.tmp | tail -n $car_natom >> $out_file
echo "$sol_x $sol_y $new_box_z"		>> $out_file

rm -f tmp_*.gro.tmp
