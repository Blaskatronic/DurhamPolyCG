# Summary #

Repository for the coarse-graining and analysis software used to obtain the charge transport properties for P3HT.

This code was used to generate the data in the following peer-reviewed article, which should be cited if this code is used in any future publications:

Jones, M. L.; Huang, D. M.; Chakrabarti, B.; Groves, C., 2016, ``Relating Molecular Morphology to Charge Mobility in Semicrystalline Conjugated Polymers'', *The Journal of Physical Chemistry C*, **120**, 4240-4250 (dx.doi.org/10.1021/acs.jpcc.5b11511)

# Description #

The code can generate systems with an arbitrary number of chains and molecular weights, coarse-graining them using the force-field described in the above article. The molecular dynamics suite LAMMPS is then used to determine an equilibrated conformation of the morphology, which can then be analyzed using the included scripts.
Of particular note are the morphology characterization scripts (which consider potential energy, proportion of crystallinity and the paracrystallinity), and the kinetic Monte Carlo simulation code to determine the charge transport characteristics.

# License #

This code is released under GPL 3.0. See LICENSE.txt for more details. Code mostly deprecated (superceded by MorphCT, for more details see: https://bitbucket.org/cmelab/morphct), but other queries should be directed to the author, Matthew Jones (mattyl.jones@dunelm.org.uk)
