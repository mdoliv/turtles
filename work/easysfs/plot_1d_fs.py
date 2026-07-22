import dadi

sfs = dadi.Spectrum.from_file("sfs_all/dadi/Lo-10.sfs")

dadi.Plotting.plot_1d_fs(sfs)
