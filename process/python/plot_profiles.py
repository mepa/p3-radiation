
# coding: utf-8

# In[1]:

get_ipython().magic(u'matplotlib inline')
import yt
import numpy as np
import matplotlib.pyplot as plt
from yt import derived_field
from yt.units import g, kpc
from yt.utilities.physical_constants import mp, me, Na


# In[2]:

print 1 / Na
print mp


# In[3]:

e_abar = me / mp
Ah = 1.0 + e_abar
Ahplu = 1.0
Ahmin = 1.0 + 2.0 * e_abar
Ahel = 4.0 + 2.0 * e_abar
Ahep = 4.0 + 1.0 * e_abar
Ahepp = 4.0
Ahtwo = 2.0 + 2.0 * e_abar
Ahtwp = 2.0 + e_abar
Adeut = 2.0 + e_abar
Adplu = 2.0
Ahd = 3.0 + 2.0 * e_abar
Aelec = e_abar


# In[4]:

print Ah
print Ahplu
print Ahmin
print Ahel
print Ahep
print Ahepp
print Ahtwo
print Ahtwp
print Adeut
print Adplu
print Ahd
print Aelec


# In[5]:

@derived_field(name = "number_density", units = "1/cm**3", force_override=True)
def _numdens(field, data):
    return data["density"] / mp
@derived_field(name = "abundance_h", units = "1", force_override=True)
def _numdens(field, data):
    return data["h   "] / Ah
@derived_field(name = "abundance_hplu", units = "1", force_override=True)
def _numdens(field, data):
    return data["hplu"] / Ahplu
@derived_field(name = "abundance_hel", units = "1", force_override=True)
def _numdens(field, data):
    return data["hel "] / Ahel
@derived_field(name = "abundance_hep", units = "1", force_override=True)
def _numdens(field, data):
    return data["hep "] / Ahep
@derived_field(name = "abundance_hepp", units = "1", force_override=True)
def _numdens(field, data):
    return data["hepp"] / Ahepp
@derived_field(name = "abundance_htwo", units = "1", force_override=True)
def _numdens(field, data):
    return data["htwo"] / Ahtwo
@derived_field(name = "abundance_htwp", units = "1", force_override=True)
def _numdens(field, data):
    return data["htwp"] / Ahtwp
@derived_field(name = "abundance_deut", units = "1", force_override=True)
def _numdens(field, data):
    return data["deut"] / Adeut
@derived_field(name = "abundance_dplu", units = "1", force_override=True)
def _numdens(field, data):
    return data["dplu"] / Adplu
@derived_field(name = "abundance_hd", units = "1", force_override=True)
def _numdens(field, data):
    return data["hd  "] / Ahd
@derived_field(name = "abundance_hmin", units = "1", force_override=True)
def _numdens(field, data):
    return data["hmin"] / Ahmin
@derived_field(name = "abundance_elec", units = "1", force_override=True)
def _numdens(field, data):
    return data["elec"] / Aelec


# In[6]:

chk = 381
if chk >= 1: prefix = '000'
if chk >= 10: prefix = '00'
if chk >= 100: prefix = '0'
if chk >= 1000: prefix = ''
outdir = "/scratch/cerberus/d4/mepa/process/python/plots/"
indir = "/scratch/cerberus/d4/mepa/data/RadCosmo_res128/stampede/"
infile = "radCosmoLW_hdf5_chk_" + prefix + str(chk)


# In[7]:

ds = yt.load(indir + infile)
dd = ds.all_data()


# In[8]:

print dd["gas", "abundance_h"][0] + dd["gas", "abundance_hplu"][0] + dd["gas", "abundance_hel"][0] + dd["gas", "abundance_hep"][0] + dd["gas", "abundance_hepp"][0] + dd["gas", "abundance_htwo"][0] + dd["gas", "abundance_htwp"][0] + dd["gas", "abundance_deut"][0] + dd["gas", "abundance_dplu"][0] + dd["gas", "abundance_hmin"][0] + dd["gas", "abundance_hd"][0] + dd["gas", "abundance_elec"][0]
print dd["flash", "h   "][0] + dd["flash", "hplu"][0] + dd["flash", "hel "][0] + dd["flash", "hep "][0] + dd["flash", "hepp"][0] + dd["flash", "htwo"][0] + dd["flash", "htwp"][0] + dd["flash", "deut"][0] + dd["flash", "dplu"][0] + dd["flash", "hmin"][0] + dd["flash", "hd  "][0] + dd["flash", "elec"][0]  


# In[9]:

redshift = "%.4f" % ds.current_redshift
print ds.current_redshift, redshift
#print sorted(ds.field_list)
print ds.print_stats()


# In[10]:

sp = ds.sphere("c", (100.0, "kpc"))


## Radius vs. Ionizing Species Abundances

# In[11]:

rmin = 1.0
rmax = 100.0
prof_h = yt.create_profile(sp, 'radius', ('flash', 'h   '),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})
prof_hel = yt.create_profile(sp, 'radius', ('flash', 'hel '),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})
prof_hplu = yt.create_profile(sp, 'radius', ('flash', 'hplu'),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})
prof_hep = yt.create_profile(sp, 'radius', ('flash', 'hep '),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})                                                              
prof_hepp = yt.create_profile(sp, 'radius', ('flash', 'hepp'),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})                                                             
                                                                


# In[12]:

radius = prof_h.x
h = prof_h['flash','h   '].value
hel = prof_hel['flash','hel '].value
hplu = prof_hplu['flash','hplu'].value
hep = prof_hep['flash','hep '].value
hepp = prof_hepp['flash','hepp'].value


# In[13]:

plt.loglog(radius, h, label="H")
plt.loglog(radius, hel, label="He")
plt.loglog(radius, hplu, label="H$^+$")
plt.loglog(radius, hep, label="He$^+$")
plt.loglog(radius, hepp, label="He$^{++}$")
plt.xlabel("Radius (kpc)")
plt.ylabel("Mass Fractions")
plt.legend()
plt.savefig(outdir + 'profile_abundance_ioniz_' + prefix + str(chk) + '_z' + redshift + '.png')


## Radius vs. Dissociating Species Abundances

# In[14]:

rmin = 1.0
rmax = 100.0
prof_h = yt.create_profile(sp, 'radius', ('flash', 'h   '),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})
prof_htwo = yt.create_profile(sp, 'radius', ('flash', 'htwo'),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})
prof_hmin = yt.create_profile(sp, 'radius', ('flash', 'hmin'),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})
prof_elec = yt.create_profile(sp, 'radius', ('flash', 'elec'),
                          units = {'radius': 'kpc'},
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})                                                                                                                         


# In[15]:

radius = prof_h.x
h = prof_h['flash','h   '].value
htwo = prof_htwo['flash','htwo'].value
hmin = prof_hmin['flash','hmin'].value
elec = prof_elec['flash','elec'].value


# In[16]:

plt.loglog(radius, h, label="H")
plt.loglog(radius, elec, label="e$^-$")
plt.loglog(radius, htwo, label="H_2")
plt.loglog(radius, hmin, label="H$^-$")
plt.xlabel("Radius (kpc)")
plt.ylabel("Mass Fractions")
plt.legend()
plt.savefig(outdir + 'profile_abundance_dissoc_' + prefix + str(chk) + '_z' + redshift + '.png')


## Radius vs. Number Density

# In[17]:

rmin = 1.0e-4
rmax = 100.0
prof_n = yt.create_profile(sp, 'radius', 'number_density',
                          units = {'radius': 'kpc'},
                          n_bins = 512,
                          extrema = {'radius': ((rmin,'kpc'),(rmax,'kpc'))})


# In[18]:

radius = prof_n.x
n = prof_n['number_density'].value


# In[19]:

plt.loglog(radius, n)
plt.xlabel("Radius (kpc)")
plt.ylabel("Number Density (cm$^{-3}$)")
plt.savefig(outdir + 'profile_numdens_' + prefix + str(chk) + '_z' + redshift + '.png')


## Number Density vs. Temperature

# In[20]:

plot = yt.PhasePlot(dd, 'number_density', 'temperature', 'cell_mass', x_bins=512, y_bins=512, weight_field=None)
plot.set_unit('cell_mass', 'Msun')


# In[21]:

plot.save(outdir + 'phase_' + prefix + str(chk) + '_z' + redshift + '.png')


## Number Density vs. Abundance

# In[22]:

plot = yt.PhasePlot(dd, 'number_density', 'htwo', 'cell_mass', x_bins=512, y_bins=512, weight_field=None)
plot.set_unit('cell_mass', 'Msun')


# In[23]:

plot.set_xlim(1.0e0, 1.0e8)
plot.set_ylim(1.0e-9, 1.0e-2)


# In[24]:

plot_zoom = yt.PhasePlot(dd, 'number_density', 'abundance_htwo', 'cell_mass', x_bins=512, y_bins=512, weight_field=None)
plot_zoom.set_unit('cell_mass', 'Msun')


# In[26]:

plot_zoom.set_xlim(1.0e0, 1.0e8)
plot_zoom.set_ylim(1.0e-9, 1.0e-2)


# In[25]:

plot.save(outdir + 'htwo_' + prefix + str(chk) + '_z' + redshift + '.png')


# In[26]:

plot = yt.PhasePlot(dd, 'number_density', 'hd  ', 'cell_mass', x_bins=512, y_bins=512, weight_field=None)
plot.set_unit('cell_mass', 'Msun')


# In[27]:

plot.save(outdir + 'hd_' + prefix + str(chk) + '_z' + redshift + '.png')


# In[28]:

plot = yt.PhasePlot(dd, 'number_density', 'elec', 'cell_mass', x_bins=512, y_bins=512, weight_field=None)
plot.set_unit('cell_mass', 'Msun')


# In[7]:

plot.save(outdir + 'elec_' + prefix + str(chk) + '_z' + redshift + '.png')


# In[ ]:



