import yt
from yt.units import mp

print mp

# TODO: implement command-line chk file numbers
dirprefix = "/scratch/cerberus/d4/mepa/"
#dirname = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128/"
dirname = dirprefix + "data/RadCosmo_res128/stampede/"
filename = "radCosmoLW_hdf5_chk_0381"
fn = dirname + filename

@derived_field(name='number density', units="1/cm**3")
def num_dens(field, data):
    return data['density'] / mp
    

ds = yt.load(fn)
sp = ds.sphere("max", (100.0, "kpc"))

plot = yt.PhasePlot(sp, "number density", "temperature", "cell_mass",
                    weight_field=None)

plot.set_unit('cell_mass', 'Msun')

plot.save()
