#!/bin/bash

for i in {1..2500}
do
   mkdir set$i
done

for i in {1..2500}
do
   cp lammps_input_file.in set$i
   cp interactions set$i
   cp initial.py set$i
   cp input_data_generate_f1.py set$i
   cd set$i/
   python3 initial.py > initial_conformation.txt
   python3 input_data_generate_f1.py > param
   #nohup lmp_serial -in lammps_input_file.in &> sim$i.txt &
   #nohup python Pf_lammps.py &> hu_$i.txt &
   sbatch -c 1 --qos=normal --job-name=LF3_$i  --wrap='lmp_serial -in lammps_input_file.in &> chain.txt'
   cd ../
done
