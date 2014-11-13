from yt.mods import * # set up our namespace

#dirname = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128/"
dirname = "/scratch/cerberus/d4/mepa/data/RadCosmo_res128/stampede/"
filename = "radCosmoLW_hdf5_chk_0377"
fn = dirname + filename # parameter file to load

pf = load(fn) # load data
pc = PlotCollection(pf) # defaults to center at most dense point
pc.add_slice("Density", 0) # 0 = x-axis
pc.add_slice("Density", 1) # 1 = y-axis
pc.add_slice("Density", 2) # 2 = z-axis
pc.set_width(1.5, 'mpc') # change width of all plots in pc
pc.save(fn) # save all plots
