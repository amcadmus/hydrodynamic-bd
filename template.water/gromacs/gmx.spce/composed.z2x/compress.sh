#!/bin/bash

source parameters.sh

if test $# -ne 2; then
    echo "usage:"
    echo "compress.sh	input	ratio"
    exit
fi

inp_file=$1
out_file=out.gro
ratio=$2

sol_trans=`echo "$car_upp + $gap" | bc -l`
car_old_posi=`echo "$sol_trans + $sol_upp + $gap" | bc -l`
car_new_posi=`echo "$car_old_posi * $ratio" | bc -l`
car_trans=`echo "$car_new_posi - $car_low" | bc -l`
editconf -translate 0 0 $car_trans -f $car_file
mv -f out.gro tmp_car.gro.tmp

car_natom=`wc $car_file | awk '{print $1}'`
car_natom=`echo "$car_natom -3" | bc`

sol_natom=`wc $sol_file | awk '{print $1}'`
sol_natom=`echo "$sol_natom -3" | bc`

tot_natom=`echo "$sol_natom + $car_natom * 2" | bc`

grep -v "CAR" $inp_file > tmp_sol.gro.tmp.1
head -n $(($sol_natom+2)) tmp_sol.gro.tmp.1 | tail -n $sol_natom > tmp_sol.gro.tmp.2
mv -f tmp_sol.gro.tmp.2 tmp_sol.gro.tmp.1
head -n 2 $sol_file > tmp_sol.gro.tmp
cat tmp_sol.gro.tmp.1 >> tmp_sol.gro.tmp
rm -f  tmp_sol.gro.tmp.1
tail -n 1 $sol_file >> tmp_sol.gro.tmp
ln -s tmp_sol.gro.tmp tmp_sol.gro
editconf -scale 1 1 $ratio -f tmp_sol.gro -o tmp.gro
mv -f tmp.gro tmp_sol.gro.tmp
rm -f tmp_sol.gro

new_box_z=`echo "$car_trans + $car_upp + $vacu" | bc -l`

new_car_upp=`echo "$car_upp * $ratio" | bc -l`
new_car_low=`echo "$new_car_upp - $car_upp + $car_low" | bc -l`
car_trans_1=`echo "$new_car_low - $car_low" | bc -l | awk '{printf "%f", $1}'`
editconf -translate 0 0 $car_trans_1 -f $car_file
mv -f out.gro tmp_car_1.gro.tmp

echo ""					> $out_file
echo "$tot_natom"			>> $out_file
head -n $(($car_natom+2)) tmp_car_1.gro.tmp| tail -n $car_natom >> $out_file
head -n $(($sol_natom+2)) tmp_sol.gro.tmp  | tail -n $sol_natom >> $out_file
head -n $(($car_natom+2)) tmp_car.gro.tmp  | tail -n $car_natom >> $out_file
echo "$sol_x $sol_y $new_box_z"		>> $out_file

rm -f tmp_*.gro.tmp

