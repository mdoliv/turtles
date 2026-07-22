import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Cm", "Lo"]

fs = dadi.Spectrum.from_file("../dadi/All.cm_lo.folded.downproj.fs")
ns = fs.sample_sizes

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum:\n")
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#================================================================================
# Calling external 2D models from the Models_2D.py script
#================================================================================

#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
prefix = "_".join(pop_ids)

#**************
#make sure to define your extrapolation grid size (based on your projections)
pts_init = sum(fs.sample_sizes)
pts = [pts_init, pts_init + 10, pts_init + 20]

#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,10,5,5]
maxiters = [3,5,10,15]
folds = [3,2,2,1]

#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

for _ in range(3):
    # Split into two populations, no migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded,
                                            reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")

    # Split into two populations, with continuous symmetric migration.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded,
                                            reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")

    # Split with no gene flow, followed by period of continuous symmetrical gene flow.
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sec_contact_sym_mig", Models_2D.sec_contact_sym_mig, rounds, 5, fs_folded=fs_folded,
                                            reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")


