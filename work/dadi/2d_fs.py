import dadi

fs = dadi.Spectrum.from_file("All.cc_ei.folded.full.fs")

dadi.Plotting.plot_single_2d_sfs(fs, vmin=1e-2, vmax=1e4)
