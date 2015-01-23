import yt
from yt.units import kpc
import matplotlib.pyplot as plt

dirprefix = "/scratch/cerberus/d4/mepa/"
#dirname = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128/"
dirname = dirprefix + "data/RadCosmo_res128/stampede/"
filename = "radCosmoLW_hdf5_chk_0381"
fn = dirname + filename

ds = yt.load(fn)
sp = ds.sphere("max", (100.0, "kpc"))

prof1 = yt.create_profile(sp, 'radius', ('flash', 'h   '),
                         units = {'radius': 'kpc'},
                         extrema = {'radius': ((0.1, 'kpc'), (100.0, 'kpc'))})
prof2 = yt.create_profile(sp, 'radius', ('flash', 'hplu'),
                         units = {'radius': 'kpc'},
                         extrema = {'radius': ((0.1, 'kpc'), (100.0, 'kpc'))})
prof3 = yt.create_profile(sp, 'radius', ('flash', 'hel '),
                         units = {'radius': 'kpc'},
                         extrema = {'radius': ((0.1, 'kpc'), (100.0, 'kpc'))})
prof4 = yt.create_profile(sp, 'radius', ('flash', 'hep '),
                         units = {'radius': 'kpc'},
                         extrema = {'radius': ((0.1, 'kpc'), (100.0, 'kpc'))})
prof5 = yt.create_profile(sp, 'radius', ('flash', 'hepp'),
                         units = {'radius': 'kpc'},
                         extrema = {'radius': ((0.1, 'kpc'), (100.0, 'kpc'))})
radius = prof1.x.value
h = prof1['flash', 'h   '].value
hplu = prof2['flash', 'hplu'].value
hel = prof3['flash', 'hel '].value
hep = prof4['flash', 'hep '].value
hepp = prof5['flash', 'hepp'].value

plt.loglog(radius, h, label='H')
plt.loglog(radius, hplu, label='H$^+$')
plt.loglog(radius, hplu, label='He')
plt.loglog(radius, hplu, label='He$^+$')
plt.loglog(radius, hplu, label='He$^{++}$')
plt.xlabel('Radius [kpc]')
plt.ylabel('Abundances')
plt.legend()

plt.savefig('test.png')
