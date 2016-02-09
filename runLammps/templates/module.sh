#!bin/sh

module purge
module load sge/current
echo "Loaded sge/current"
module load pgi/current
echo "Loaded pgi/current"
module load lammps/intel/8Aug12
echo "Loaded lammps/intel8Aug12"
echo "Loaded fftw/intel/2.1.5"
module unload openmpi/intel_12.1.4/1.6
module load openmpi/intel/1.6.1
echo "Loaded openmpi/intel/1.6.1"