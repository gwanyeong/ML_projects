from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter, BoltztrapPlotter

from pymatgen.electronic_structure.boltztrap2 import *
from monty.serialization import loadfn


# data = VasprunLoader().from_file('vasprun.xml')

# New unique loader
v = Vasprun('./vasprun.xml', parse_projected_eigen = True)
data = VasprunBSLoader(v)

"""
bs = v.get_band_structure(kpoints_filename = './KPOINTS', line_mode = True)
bs_loader = BandstructureLoader(bs)
"""

bztInterp = BztInterpolator(data, lpfac=10, energy_range=1.5, curvature=True)

# set fname argument to specify a different file name
bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(100,350,50), doping=10.**np.arange(15,18.5,0.5),
                                   save_bztTranspProps=True, fname='bztTranspProps.json.gz')

# bztTransp = BztTransportProperties(bztInterp,load_bztTranspProps=True,fname='bztTranspProps.json.gz')

print('\t'.join(['Temp', '\mu', 'rows', 'columns tensor']))
for p in bztTransp.Conductivity_mu, bztTransp.Seebeck_mu, bztTransp.Kappa_mu, \
         bztTransp.Effective_mass_mu, bztTransp.Power_Factor_mu, bztTransp.Carrier_conc_mu:
    print('\t'.join([str(i) for i in p.shape]))

bztPlotter = BztPlotter(bztTransp, bztInterp)
plt1 = bztPlotter.plot_props('C','temp', 'doping', dop_type = 'n')
plt1.ylim(0,1000)
plt1.legend(['%.0e' % 1e15, '%.0e' % 5e15, 1e16, 5e16, 1e17, 5e17, 1e18])
plt1.savefig('conductivity_vs_T.jpg', dpi = 400)

plt2 = bztPlotter.plot_props('C','doping','temp')
# plt = bztPlotter.plot_props('C','mu','temp')# , ylim=[None,1.2])
# plt.ylim(0*1e6,1.2*1e6)
plt2.ylim(0,1000)
plt2.savefig('conductivity_vs_n.jpg', dpi = 400)

plt3 = bztPlotter.plot_props('E', 'temp', 'doping', dop_type = 'n')
plt3.ylim(0,1)
plt3.legend(['%.0e' % 1e15, '%.0e' % 5e15, 1e16, 5e16, 1e17, 5e17, 1e18])
plt3.savefig('effective_mass_vs_T.jpg', dpi = 400)

plt4 = bztPlotter.plot_props('E','doping','temp')
plt4.ylim(0,1)
plt4.savefig('effective_mass_vs_n.jpg', dpi = 400)



# plt2 = bztPlotter.plot_props('S','mu','temp')
# plt2.ylim(-1500,1500)
#plt2 = bztPlotter.plot_props('S','temp','doping', doping=[1e17, 1e18, 1e19], dop_type = 'n')
#plt2.savefig('Seebeck.pdf')

# plt3 = bztPlotter.plot_props('K','mu','temp')
#plt3 = bztPlotter.plot_props('K','temp','doping', doping=[1e17, 1e18, 1e19], dop_type = 'n')
#plt3.savefig('K_el.pdf')

# plt4 = bztPlotter.plot_props('H','mu','temp')
# plt4 = bztPlotter.plot_props('H', 'temp', 'doping', doping = [1e17, 1e18, 1e19], dop_type = 'n')
# plt4.ylim(1e11, 1e29)
# plt4.savefig('Hall_carrier_conc.pdf')

# plt5 = bztPlotter.plot_props('Ca','mu','temp')
# plt5 = bztPlotter.plot_props('Ca', 'temp', 'doping', doping = [1e17, 1e18, 1e19], dop_type = 'n')
# plt5.savefig('power_factor.pdf')

# plt6 = bztPlotter.plot_props('E', 'mu', 'temp')

# plt7 = bztPlotter.plot_bands()
#plt7.savefig('bands.pdf')

# plt8 = bztPlotter.plot_dos(T=130)
# plt8.savefig('DOS.pdf')

