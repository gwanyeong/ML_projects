# my_projects

## 1. KRICT_ML
: Photoelectrochemical water oxidation reaction for highly active and selective H<sub>2</sub>O<sub>2</sub> generation

#### Target data
- Bulk models : ~1,000
- Surface models : ~3,000

#### Target materials
- Perovskite oxides : ABO3
- Layered transition metal (oxy)(hydro)oxides : MOOH, LDH, LTMO
- Transition metal oxides : MO, MO2, M2O3

#### Target properties
- Bulk data : formation enthalpy, band gap, magnetic property, volume
- Surface data : surface energy, adsorption energy (OOH*, O*, OH*)
-> _faradaic efficiency, thermodynamic overpotentials_

### - Workflow for DFT calculations
##### 1) Bulk modeling (from MP or directly)
##### 2) Bulk model analysis (bandgap, DH<sub>f</sub>, magnetization, lattice ..)
##### 3) Slab modeling (most stable or representatively active surface)
##### 4) Adsorbate modeling (OOH*, O* OH*)
##### 5) Adsorbate analysis : calculation of adsorption free energies -> volcanos(2e-,4e-) for h<sub>TD</sub>, f(%)
##### 6) ML-based design of active, selective catalyst using our constructed databse (matminer, automatminer)

## 2. Zwitterionic
: zwitterionic COF as solid electrolyte
- construction of PES (vertical, lateral shift)
- band structure, dos, etc..


## 3. scripts
: useful scripts for analyzing VASP results or any other perl scripts 
