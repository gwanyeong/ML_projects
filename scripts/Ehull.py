import os
from ase.io import read

from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.io.vasp import Oszicar
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram, PhaseDiagram, GrandPotentialPhaseDiagram, PDPlotter

# import chart_studio
# from chart_studio.plotly import plot, iplot
##############################################################################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def computed_entry(path):
    atoms = read(path + '/' + 'CONTCAR')
    comp = Composition(str(atoms.symbols))

    oszicar = Oszicar(path + '/' +'OSZICAR')
    energy  = oszicar.ionic_steps[-1]['E0']
#   tot_mag = oszicar.ionic_steps[-1]['mag']

    comp_dict = {}
    for i in comp:
        comp_dict[str(i)] = int(comp[i])

    # Calculate corrections
    element_list = [i for i in comp_dict.keys()]
    correction = 0
    for i in element_list:
        if i in list(U_corr_dict):
            correction += U_corr_dict[i] * comp_dict[i]

    return ComputedEntry(composition = comp, energy = energy, correction = correction, entry_id = 'gyjung')


def get_phase_diagram(path):
    atoms = read(path + '/' + 'CONTCAR')
    comp = Composition(str(atoms.symbols))

    comp_dict = {}
    for i in comp:
        comp_dict[str(i)] = int(comp[i])

    element_list = [i for i in comp_dict.keys()]

    entries = mpr.get_entries_in_chemsys(element_list)
    new_entries = entries + [computed_entry(path)]

    createFolder('plots')
    phase_diagram = PhaseDiagram(new_entries)
#    plotter = PDPlotter(phase_diagram, show_unstable = False)
#    plt = plotter.get_contour_pd_plot()
#    plt.tight_layout()
#    plt.savefig('plots/contour_phase_diagram.pdf', dpi = 300, format = 'pdf')

#    fig = plotter.get_plot(phase_diagram, label_unstable = False)
#   fig.show()
#    fig.write_image(file= 'plots/phase_diagram.pdf', format = 'pdf')

    Ehull = phase_diagram.get_e_above_hull(computed_entry(path))
    with open('plots/Ehull', 'w') as f:
        f.writelines("%4.6f eV" % Ehull)
    print(Ehull)

    return Ehull

##############################################################################

mpr = MPRester('MPI_key')

U_corr_dict = {'V':-1.700, 'Cr':-1.999, 'Mn':-1.668, 'Fe':-2.256, 'Co':-1.638, 'Ni':-2.541, 'Mo':-3.202, 'W':-4.438, 
               'O' : -0.687, 'F':-0.462, 'I':-0.379}

# chart_studio.tools.set_credentials_file(username='gwanyeong', api_key = 'api_key')

path = os.getcwd()
get_phase_diagram(path)
