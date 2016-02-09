#!/bin/csh

# simple template job script for submission of a serial job with
# SLURM directives

# In this script the minimum requirements have been set for SLURM and 
# it only requires editing the last line and replace 

# ./my_serial_program 

# with the name of the serial program to execute.
# The job can be submitted with the command

# sbatch -p  name_of_serial_queue  name_of_job_script

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

# SBATCH  -o  remove.out
# SBATCH  -e  remove.err

## where SLURM will expand %j to the jobID of the job.

## Request email to the user when certain type of events occur to 
## the job

# SBATCH  --mail-type=<type>

## where <type> can be one of BEGIN, END, FAIL, REQUEUE or ALL,
## and send to email address

# SBATCH  --mail-user  email_name@durham.ac.uk

## The default email name is that of the submitting user as known to the system.

## Specify project or account name (currently not required).
##SBATCH -A ITSall

#############################################################################

## The number of requested processes is specified with the -n option and 
## is 1 for a serial job

##SBATCH  -n  1

## Request submission to a queue (partition) for serial jobs

##SBATCH  -p  seq6.q

module purge
module load slurm/current
module load python/2.7.1

module li
## Load any other modules required here

## Execute serial program

python parallelRunOrca.py -m Mw_57_PDI_1.9_con_40_Equil_290K -c 0 -e None -q None


