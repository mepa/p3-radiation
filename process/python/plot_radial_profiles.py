# TODO: implement command-line chk file numbers
dirprefix = "/scratch/cerberus/d4/mepa/"
#dirname = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128/"
dirname = dirprefix + "data/RadCosmo_res128/stampede/"
filename = "radCosmoLW_hdf5_chk_0381"
fn = dirname + filename # parameter file to load

# Load the dataset.
ds = yt.load(fn)

v, c = pf.h.find_max("Density")
print "Found highest density of %0.3e at %s" % (v, c)

# Create a sphere of radius 100 kpc around the cell of max density.
my_sphere = ds.sphere(c, (100.0, "kpc"))

# Create a profile of the average density vs. radius.
plot = yt.ProfilePlot(my_sphere, "radius", "density",
                      weight_field="cell_mass")

# Change the units of the radius into kpc (and not the default in cgs)
plot.set_unit('radius', 'kpc')

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save(dirprefix + 'process/python/plots/' + filename)
