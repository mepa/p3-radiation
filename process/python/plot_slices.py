import yt
from yt.units import kpc


dirprefix = "/scratch/cerberus/d4/mepa/"
#dirname = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128/"
dirname = dirprefix + "data/RadCosmo_res128/stampede/"
filename = "radCosmoLW_hdf5_chk_0381"
fn = dirname + filename # parameter file to load

ds = yt.load(fn) # load data
slc = yt.SlicePlot(ds, "x", "density")
slc.zoom(10)
slc.save(dirprefix + 'process/python/plots/zoom_' + filename)
slc = yt.SlicePlot(ds, "y", "density", center='m')
slc.zoom(10)
#slc.set_width(6*kpc)
slc.save(dirprefix + 'process/python/plots/zoom_' + filename)
slc = yt.SlicePlot(ds, "z", "density")
slc.zoom(10)
slc.save(dirprefix + 'process/python/plots/zoom_' + filename)

