import dadi
from itertools import product
from datetime import datetime

#===========================================================================
# Function to generate and test all projection combinations, prints to screen
#===========================================================================

def run_projections(dd, pop_ids, maxproj, min_frac):
    # get minprojections for all maxproj items using fraction
    minproj = [int(x * min_frac) for x in maxproj]
    # initiate empty list to store projection ranges
    sizes = []
    # iterate over each item in the maxprojections list
    for i in range(0, len(maxproj)):
        # create a list of values, max to min by decreasing increments of 2
        sizes.append(list(range(maxproj[i], minproj[i], -2)))
    # create all possible combinations of projections using itertools product module
    sizecombinations = list(product(*sizes))
    # print info for this jsfs and projections
    print("\n\n\n{}\nRunning test projections for {}D JSFS containing: {}\n".format("-"*80, len(pop_ids), ", ".join(pop_ids)))
    print("Maximum projection sizes = {}\nMinimum projection sizes = {}".format(maxproj, minproj))
    print("\n\nFound {} projection combinations...".format(len(sizecombinations)))

    # list to store projection sizes and segregating sites
    results = []
    # iterate over all projection combinations
    for combo in sizecombinations:
        # create a spectrum using projection sizes
        fs = dadi.Spectrum.from_file(dd)
        fs = fs.fold()
        #print("projection: {},\tsites: {}".format(combo, int(fs.S())))
        # add sizes and segregating sites in a sub-list to results list
        results.append([combo, int(fs.S())])
        
    # show results by descending order of projections
    print("\n\nResults sorted by combination order:\n")
    for i in range(0, len(results)):
        if i != 0:
            if results[i][0][0] != results[i-1][0][0]:
                print("\n")
        print("\t{1:,} segregating sites with projection sizes of {0}".format(results[i][0], results[i][1]))

    # show results by highest to lowest number segregating sites
    results.sort(key=lambda x: x[1], reverse=True)
    print("\n\nResults sorted by highest number of segregating sites:\n")
    for r in results:
        print("\t{1:,} segregating sites with projection sizes of {0}".format(r[0], r[1]))
    print("\n{}\n\n".format("-"*80))


#===========================================================================
# Example for a 2D joint-site frequency spectrum
#===========================================================================

#**************
snps = "../realsfs/rDerCor1.pri.v4.mean_depth_5.cc_x_ei.unfolded.sfs"

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["cc", "ei"]

#**************
# Maximum projection sizes. This is the maximum possible ALLELES (not individuals) that can be
# found in each of the populations. 
maxproj = [142, 86]

# Choose the minimum fraction of alleles required for both populations.
# For example, if you want to include a minimum of 50% of the alleles, use min_frac = 0.5.
# For a minimum of 75% of alleles, use min_frac = 0.75, etc.
min_frac = 0.5

run_projections(snps, pop_ids, maxproj, min_frac)
