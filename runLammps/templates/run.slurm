#!/bin/csh

# simple template job script for submission of an MPI job with
# SLURM directives

# In this script the minimum requirements have been set for SLURM 
# except for these two changes that have to be made
# 
# 1. replacing <number_of_tasks> for the -n option by the actual 
#    number of MPI tasks required
# 2. replacing my_mpi_program by the actual executable name

# The job can be submitted with the command

# sbatch -p parallel_queue  name_of_job_script

# or with overriding the number of tasks as option

# sbatch -p parallel_queue -n number_of_tasks  name_of_job_script

# If successful SLURM will return a jobID for this job which can be
# used to query its status.

#############################################################################

## All lines that start with #SBATCH will be processed by SLURM.
## Lines in this template script that have white space between # and SBATCH 
## will be ignored. They provide examples of further options that can be
## activated by deleting the white space and replacing any text after the 
## option.

## By default SLURM uses as working directory the directory from where the
## job script is submitted. To change this the standard Linux cd command has
## to be used.

## Name of the job as it will appear when querying jobs with squeue (the
## default is the name of the job script)

# SBATCH  -J  job_name

## By default SLURM combines the standard output and error streams in a single
## file based on the jobID and with extension .out
## These streams can be directed to separate files with these two directives

# SBATCH  -o  out_file_name.o%j
# SBATCH  -e  err_file_name.e%j

## where SLURM will expand %j to the jobID of the job.

## Request email to the user when certain type of events occur to 
## the job

#SBATCH  --mail-type=ALL

## where <type> can be one of BEGIN, END, FAIL, REQUEUE or ALL,
## and send to email address

#SBATCH  --mail-user  m.l.jones@durham.ac.uk

## The default email name is that of the submitting user as known to the system.

## Specify project or account name (currently not required).
##SBATCH -A ITSall

#############################################################################

## This job requests number_of_tasks MPI tasks (without OpenMP)

#SBATCH  -n  8

# Request submission to a queue (partition) for parallel jobs

##SBATCH  -p  parallel_queue


module purge
module load sge/6.1u6
module load blas/gcc/64/1
module load lapack/gcc/3.11
module load python/2.7.1
module load slurm/14.03.7
module load intel/14.0.2
module load intelmpi/intel/4.1.0
module load lammps/intel/1Feb14

module li

# execute command
mpirun -np 8 ./lmp_ham < run.in
