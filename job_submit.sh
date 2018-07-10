#!/usr/bin/env zsh
 
### Job name
#BSUB -J 2x1v_moments
 
### File / path where STDOUT & STDERR will be written
###    %J is the job ID, %I is the array ID
#BSUB -o log_files/2x1v_moments_Inflow_M4
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 1:00
 
### Request memory you need for your job per PROCESS in MB
#BSUB -M 2024
 
### Request the number of compute slots you want to use
#BSUB -n 64
 
### Use esub for Open MPI
#BSUB -a openmpi
 

switch intel/16.0 gcc/system-default
module load openmpi/1.10.2

cd /home/ns179556/DG_MPI
 
export FLAGS_MPI_BATCH="-np 8"
### Execute your application
mpirun $FLAGS_MPI_BATCH ./2x1v_Inflow.out 8 500 4
