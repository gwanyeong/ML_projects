# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 13:41:46 2020
@author: gyjung
"""


import os
import warnings
warnings.filterwarnings('ignore')

"""
Make plot directory
"""
try:
    if not(os.path.isdir('plots')):
        os.makedirs(os.path.join('plots'))
except OSError as e:
    if e.errno != errno.EEXIST:
        print("WARNING: plots directory already exists")
        raise

"""
Plot DOS
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

v = Vasprun('./vasprun.xml', parse_dos = True)
cdos = v.complete_dos
element_dos = cdos.get_element_dos()
plotter = DosPlotter()
plotter.add_dos_dict(element_dos)
plotter.save_plot('plots/dos.pdf', img_format='pdf', xlim= None, ylim = None)
# plotter.save_plot('spin-up_dos.pdf', img_format='pdf', xlim= None, ylim = [0,None])   # up-spin dos


"""
Plot Band
"""

from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter


v = BSVasprun('./vasprun.xml', parse_projected_eigen = True)
bs = v.get_band_structure(kpoints_filename = './KPOINTS', line_mode = True)
bsplot = BSPlotter(bs)

bsplot.get_plot(zero_to_efermi = True, ylim = [-5,5]).savefig('plots/band.pdf')

# add some features
ax = plt.gca()
#ax.set_title("Bandstructure", fontsize = 40) # use this for setting title
xlim = ax.get_xlim()
ax.hlines(0, xlim[0], xlim[1], linestyles="dashed", color="black")

# add legend
ax.plot((), (), "b-", label = "spin up")
ax.plot((), (), "r--", label = "spin-down")
ax.legend(fontsize = 16, loc = "upper left")

bs.as_dict()['vbm']
bs.as_dict()['cbm']

vbm = bs.as_dict()['vbm']['energy']
efermi = bs.as_dict()['efermi']
cbm = bs.as_dict()['cbm']['energy']


if bs.get_band_gap()['direct'] == False:
    bgtype = 'Indirect'
else:
    bgtype = 'Direct'

bg = bs.get_band_gap()['energy']
bgpath = bs.get_band_gap()['transition']

"""
Plot Brillouin Zone
"""
bsplot.plot_brillouin()
plt.savefig('plots/brillouin_fig.pdf')



"""
# Another plot type
bsplot_proj = BSPlotterProjected(bs)
bsplot_proj.get_elt_projected_plots(vbm_cbm_marker = True)
# bsplot_proj.get_projected_plots_dots({"Ti":["s","p","d"], "O":["s", "p"]}, ylim = [-5,5], vbm_cbm_marker = True)
"""


"""
# Getting raw data for specific bands & kpts
from pymatgen import Spin
data = bsplot.bs_plot_data()
ibands = 9
spin = str(Spin.up)
for xpath, epath in zip(data['distances'], data['energy']):
    print(20 * "-")
    for x, bands in zip(xpath, epath[spin][ibands]):
        print("%8.4f %8.4f" % (x,bands))
"""

"""
Plot DOS & Band
"""

from pymatgen.electronic_structure.plotter import BSDOSPlotter

bsdosplot = BSDOSPlotter(
        bs_projection = "elements",
        dos_projection = "elements",
        vb_energy_range = 2,
        cb_energy_range = 2,
        egrid_interval = 2,
        font = 'DejaVu Sans'
        )

bsdosplot.get_plot(bs, cdos).savefig('plots/band_dos.pdf')


"""
Summary
"""

with open('plots/summary', 'w') as f:
    if (type(vbm) and type(cbm)) == float:
        f.writelines("VBM = %4.3f eV,  E_fermi = %4.3f eV, CBM = %4.3f eV\n" % (vbm, efermi, cbm))
        f.writelines("%s gap = %4.3f eV (transition = %s)" % (bgtype, bg, bgpath))
    elif (vbm and cbm) == None:
        f.writelines("VBM = None,  E_fermi = %4.3f eV, CBM = None\n" % (efermi))
        f.writelines("%s gap = %4.3f eV (transition = %s)" % (bgtype, bg, bgpath))


"""
Plot COHP
"""

from pymatgen.electronic_structure.plotter import CohpPlotter
from pymatgen.electronic_structure.cohp import  CompleteCohp
from pymatgen.io.lobster import Icohplist

file_list = [file for file in os.listdir() if file.endswith('.lobster')]

if file_list != []:

    COHPCAR_path = os.getcwd() + "/COHPCAR.lobster"
    POSCAR_path = os.getcwd() + "/POSCAR"
    completecohp = CompleteCohp.from_file(fmt="LOBSTER",filename=COHPCAR_path,structure_file=POSCAR_path)


    # Read in ICOHPLIST.lobster and get lcohpcollection object
    icohplist = Icohplist(filename = os.getcwd() + '/ICOHPLIST.lobster')

    idx_list = [idx for idx in list(icohplist.icohplist)]

    # search for the number of the COHP you would like to plot in ICOHPLIST.lobster
    # (the numbers in COHPCAR.lobster are different!)

    cp = CohpPlotter()

    # make COHPs directory
    try:
        if not(os.path.isdir('COHPs')):
            os.makedirs(os.path.join('COHPs'))
    except OSError as e:
        if e.errno != errno.EEXIST:
            print("WARNING: COHPs directory already exists")
            raise


    for idx in idx_list:
        species1 = completecohp.bonds[idx]['sites'][0]
        species2 = completecohp.bonds[idx]['sites'][1]
        plotlabel = str(species1.species_string) + '-' + str(species2.species_string)

        cp.add_cohp(plotlabel, completecohp.get_cohp_by_label(label = idx))

        # check which COHP you are plotting
        print("This is a COHP between the following sites: " + str(species1) + ' and' + str(species2))

        cohpplot = cp.get_plot(integrated = False)
        cohpplot.ylim([-10, 6])
        cohpplot.savefig('COHPs/cohp_no_%s_%s-%s.pdf' % (idx, species1.species_string, species2.species_string))

else:
    print('lobster files are not found !!')
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
~                                                              
