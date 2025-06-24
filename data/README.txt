Standard Test Cases for DL_POLY_5
----------------------------------

Because of the size of the data files for the DL_POLY_5
standard test cases, they are not shipped in the standard
download of the DL_POLY_5 source. Test files are downloaded
automatically when building/ running CMake with the CMake variable BUILD_TESTING=ON.
This can be done as follows:

folder="build-mpi-testing"
rm -rf $folder && mkdir $folder && pushd $folder
cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON 

Unpack the files in the `data' subdirectory using firstly
`gunzip' to uncompress them and then `tar -xf' to create the
`TEST_X' directory.

Version controlled test cases may be downloaded from the ccp5UK
Github: https://github.com/ccp5UK/dl-poly-data

The following is a list of DL_POLY_5 test cases currently used:

TEST001 -  Sodium Chloride (27,000 ions) interacting via Born-Huggins-Meyer potentials
TEST002 -  DMPC in Water (51,737 atoms)
TEST003 -  KNaSi2O5 - Potassium/Sodium Disilicate Glass (69,120 ions)
TEST004 -  Gramicidin A Molecules in Water (99,120 ions)
TEST005 -  SiC with Tersoff Potential (74,088 atoms)
TEST006 -  Cu3Au alloy with Sutton-Chen (metal) Potentials (32,000 atoms)
TEST007 -  Lipid Bilayer in Water (12,428 ions)
TEST008 -  MgO with Adiabatic Shell Model (8,000 ions + 4,000 shells)
TEST009 -  MgO with Relaxed Shell Model (8,000 ions + 4,000 shells)
TEST010 -  Potential of Mean Force on K+ (500 ions) in Water (13,000 ions)
TEST011 -  Cu3Au Alloy with Gupta (metal) Potentials (32,000 atoms)
TEST012 -  Cu with EAM (metal) Potential (32,000 atoms)
TEST013 -  Al with Analytic Sutton-Chen (metal) Potential (32,000 atoms)
TEST014 -  Al with EAM Sutton-Chen (metal) Potential (32,000 atoms)
TEST015 -  NiAl Alloy with EAM (metal) Potentials (27,648 atoms)
TEST016 -  Fe with Finnis-Sinclair (metal) Potential (31,250 atoms)
TEST017 -  Ni with EAM (metal) Potential (32,000 atoms)
TEST018 -  SPC IceVII Water with CBs (34,992 ions)
TEST019 -  SPC IceVII Water with RBs (34,992 ions)
TEST020 -  NaCl Molecules in SPC Water Represented as CBs+RBs (26,816 ions)
TEST021 -  TIP4P water: RBs with a Massless Charged Site (29,052 particles)
TEST022 -  Ionic Liquid Dimethylimidazolium Chloride as RBs (44,352 ions)
TEST023 -  Calcite Nano-Particles (600) in TIP3P Water (6,904 molecules, 23,712 ions)
TEST024 -  Iron/Carbon Alloy with EEAM (metal) Potentials (36,803 particles)
TEST025 -  Iron/Chromium Alloy with 2BEAM (metal) Potentials (32,000 particles)
TEST026 -  Butane in CCl4 Solution with Umbrella Sampling via PLUMED
TEST027 -  Iron with tabulated EAM (metal) Potential (54,000 atoms), Two-Temperature Model, Cascade
TEST028 -  Silicon with original Tersoff (T3) Potential (200,000 atoms), Two-Temperature Model, Swift heavy ion irradiation
TEST029 -  Tungsten with extended Finnis-Sinclair Potential (722,672 atoms), Two-Temperature Model, laser-induced relaxation of thin film
TEST030 - Sodium Chloride (27,000 ions), NVE ensemble
TEST031 - Sodium Chloride (27,000 ions), NVT ensemble, Evans thermostat
TEST032 - Sodium Chloride (27,000 ions), NVT ensemble, Langevin thermostat
TEST033 - Sodium Chloride (27,000 ions), NVT ensemble, Andersen thermostat
TEST034 - Sodium Chloride (27,000 ions), NVT ensemble, Nose-Hoover thermostat
TEST035 - Sodium Chloride (27,000 ions), NVT ensemble, Gentle-Stochastic thermostat
TEST036 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of first order, no drag coefficient
TEST037 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of first order, with drag coefficient
TEST038 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of second order, no drag coefficient
TEST039 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of second order, with drag coefficient
TEST040 - Sodium Chloride (27,000 ions), NPT ensemble, Langevin thermostat-barostat
TEST041 - Sodium Chloride (27,000 ions), NPT ensemble, Berendsen thermostat-barostat
TEST042 - Sodium Chloride (27,000 ions), NPT ensemble, Nose-Hoover thermostat-barostat
TEST043 - Sodium Chloride (27,000 ions), NPT ensemble, Martyna-Tuckerman-Klein thermostat-barostat
TEST044 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Langevin thermostat-barostat
TEST045 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Berendsen thermostat-barostat
TEST046 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Nose-Hoover thermostat-barostat
TEST047 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Martyna-Tuckerman-Klein thermostat-barostat
TEST048 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Langevin thermostat /semi-isotropic barostat : constant normal pressure (Pn)
TEST049 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Berendsen thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST050 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Nose-Hoover thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST051 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST052 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Langevin thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST053 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Berendsen thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST054 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Nose-Hoover thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST055 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST056 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Langevin thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST057 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Berendsen thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST058 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Nose-Hoover thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST059 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST060 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Langevin thermostat/semi-isotropic barostat    : orthorhombic MD cell constraints
TEST061 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Berendsen thermostat/semi-isotropic barostat   : orthorhombic MD cell constraints
TEST062 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Nose-Hoover thermostat/semi-isotropic barostat : orthorhombic MD cell constraints
TEST063 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : orthorhombic MD cell constraints
TEST064 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Langevin thermostat/semi-isotropic barostat    : semi-orthorhombic MD cell constraints
TEST065 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Berendsen thermostat/semi-isotropic barostat   : semi-orthorhombic MD cell constraints
TEST066 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Nose-Hoover thermostat/semi-isotropic barostat : semi-orthorhombic MD cell constraints
TEST067 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : semi-orthorhombic MD cell constraints
TEST068 - Impact and defects analysis for rdiation damage in Sodium Chloride (27,000 ions)
TEST069 - Sodium Chloride (27,000 ions) interacting via 12-6 potentials
TEST070 - Sodium Chloride (27,000 ions) interacting via Lenhard-Jones potentials
TEST071 - Sodium Chloride (27,000 ions) interacting via Lenhard-Jones cohesive potentials
TEST072 - Sodium Chloride (27,000 ions) interacting via n-m potentials
TEST073 - Sodium Chloride (27,000 ions) interacting via Buckingham potentials
TEST074 - Sodium Chloride (27,000 ions) interacting via Hydrogen-bond potentials
TEST075 - Sodium Chloride (27,000 ions) interacting via Shifted force n-m potential potentials
TEST076 - Sodium Chloride (27,000 ions) interacting via Morse potentials
TEST077 - Sodium Chloride (27,000 ions) interacting via Shifted Weeks-Chandler-Andersen potentials
TEST078 - Sodium Chloride (27,000 ions) interacting via Standard DPD potentials
TEST079 - Sodium Chloride (27,000 ions) interacting via 14-7 pair potentials
TEST080 - Sodium Chloride (27,000 ions) interacting via Morse modified potentials
TEST081 - Sodium Chloride (27,000 ions) interacting via Rydberg potentials
TEST082 - Sodium Chloride (27,000 ions) interacting via Ziegler-Biersack-Littmark (ZBL) potentials
TEST083 - Sodium Chloride (27,000 ions) interacting via ZBL mixed with Morse potentials
TEST084 - Sodium Chloride (27,000 ions) interacting via ZBL mixed with Buckingham potentials
TEST085 - Sodium Chloride (27,000 ions) interacting via tabulated potentials
TEST086 - Sodium Chloride (27,000 ions) interacting via Lennard-Jones tapered with Mei-Davenport-Fernando (MDF) taper potentials
TEST087 - Sodium Chloride (27,000 ions) interacting via Buckingham tapered with MDF potentials
TEST088 - Sodium Chloride (27,000 ions) interacting via 12-6 tapered with MDF potentials
TEST089 -  Restart of TEST004
TEST090 -  Formic Acid (970 molecules)-Test for bond potentials from tabulated values
TEST091 -  Formic Acid (970 molecules)-Test for bond potentials of tab+vdW+coul (-tab) type
TEST092 -  Formic Acid (970 molecules)-Test for bond potentials of harmonic (harm) type
TEST093 -  Formic Acid (970 molecules)-Test for bond potentials of harmonic+vdW+coul (-hrm) type
TEST094 -  Formic Acid (970 molecules)-Test for bond potentials of Morse (mors) type
TEST095 -  Formic Acid (970 molecules)-Test for bond potentials of Morse+vdW+coul (-mrs) type
TEST096 -  Formic Acid (970 molecules)-Test for bond potentials of 12-6 type
TEST097 -  Formic Acid (970 molecules)-Test for bond potentials of (12-6)+vdW+coul (-126) type
TEST098 -  Formic Acid (970 molecules)-Test for bond potentials of Lennard-Jones (lj) type
TEST099 -  Formic Acid (970 molecules)-Test for bond potentials of (lj)+vdW+coul (-lj) type
TEST100 -  Formic Acid (970 molecules)-Test for bond potentials of restrained harmonic (rhrm) type
TEST101 -  Formic Acid (970 molecules)-Test for bond potentials of restricted harmonic+vdW+coul (-rhm) type
TEST102 -  Formic Acid (970 molecules)-Test for bond potentials of quar type
TEST103 -  Formic Acid (970 molecules)-Test for bond potentials of quar+vdW+coul (-qur) type
TEST104 -  Formic Acid (970 molecules)-Test for bond potentials of Buckingham (buck) type
TEST105 -  Formic Acid (970 molecules)-Test for bond potentials of Buckingham+vdW+coul (-bck) type
TEST106 -  Formic Acid (970 molecules)-Test for bond potentials of Coulomb (coul) type
TEST107 -  Formic Acid (970 molecules)-Test for bond potentials of Electroctatic+vdW+coul (-cul) type
TEST108 -  Formic Acid (970 molecules)-Test for bond potentials of fene type
TEST109 -  Formic Acid (970 molecules)-Test for bond potentials of FENE+vdW+coul (-fne) type
TEST110 -  Formic Acid (970 molecules)-Test for bond potentials of mmts type
TEST111 -  Formic Acid (970 molecules)-Test for bond potentials of mmst+vdW+coul (-mst) type
TEST112 -  Formic Acid (970 molecules)-Test for bond tabulated (tab) potentials
TEST113 -  Formic Acid (970 molecules)-Test for tabulated angle potentials of tab+coul+vdW (-tab) type
TEST114 -  Formic Acid (970 molecules)-Test for angle potentials of harmonic (harm) type
TEST115 -  Formic Acid (970 molecules)-Test for angle potentials of harmonic+coul+vdW (-hrm) type
TEST116 -  Formic Acid (970 molecules)-Test for angle potentials of quartic (quartic) order
TEST117 -  Formic Acid (970 molecules)-Test for angle potentials of quartic+coul+vdW (-qur) type
TEST118 -  Formic Acid (970 molecules)-Test for angle potentials of Truncated harmonic (thrm) type
TEST119 -  Formic Acid (970 molecules)-Test for angle potentials of Truncate harmonic+coul+vdW (-thm) type
TEST120 -  Formic Acid (970 molecules)-Test for angle potentials of Screened harmonic (shrm) type
TEST121 -  Formic Acid (970 molecules)-Test for angle potentials of Screnned harmonic+coul+vdW (-shm) type
TEST122 -  Formic Acid (970 molecules)-Test for angle potentials of Screened Vessal (bvs1) type
TEST123 -  Formic Acid (970 molecules)-Test for angle potentials of Screnned Vessal+coul+vdW (-bv1) type
TEST124 -  Formic Acid (970 molecules)-Test for angle potentials of Truncated Vessal (bvs2) type
TEST125 -  Formic Acid (970 molecules)-Test for angle potentials of Truncate Vessal+coul+vdW (-bv2) type
TEST126 -  Formic Acid (970 molecules)-Test for angle potentials of harmonic cosine (hacos) type
TEST127 -  Formic Acid (970 molecules)-Test for angle potentials of harmonic cosine+coul+vdW (-hcs) type
TEST128 -  Formic Acid (970 molecules)-Test for angle potentials of cosine (cos) type
TEST129 -  Formic Acid (970 molecules)-Test for angle potentials of cosine+coul+vdW (-hcs) type
TEST130 -  Formic Acid (970 molecules)-Test for angle potentials of MM3 stretch-bend (mmsb) type
TEST131 -  Formic Acid (970 molecules)-Test for angle potentials of MM3 stretch-bend+coul+vdW (-msb) type
TEST132 -  Formic Acid (970 molecules)-Test for angle potentials of compass stretch-stretch (stst) type
TEST133 -  Formic Acid (970 molecules)-Test for angle potentials of Compass stretch-stretch +coul+vdW (-sts) type
TEST134 -  Formic Acid (970 molecules)-Test for angle potentials of compass stretch-bend (stbe) type
TEST135 -  Formic Acid (970 molecules)-Test for angle potentials of Compass stretch-bend +coul+vdW (-stb) type
TEST136 -  Formic Acid (970 molecules)-Test for angle potentials of compass all-terms (cmps) type
TEST137 -  Formic Acid (970 molecules)-Test for angle potentials of Compass all +coul+vdW (-cmp) type
TEST138 -  Formic Acid (970 molecules)-Test for angle potentials of mmbd type
TEST139 -  Formic Acid (970 molecules)-Test for angle potentials of mmbd+coul+vdW (-mbd) type
TEST140 -  Formic Acid (970 molecules)-Test for tabulated potentials for dihedrals
TEST141 -  Formic Acid (970 molecules)-Test for cosine potentials for dihedrals
TEST142 -  Formic Acid (970 molecules)-Test for dihedrals potentials of harmonic (harm) type
TEST143 -  Formic Acid (970 molecules)-Test for dihedrals of harmonic cosine (hcos) type
TEST144 -  Formic Acid (970 molecules)-Test for dihedrals of harmonic cosine (cos3) type
TEST145 -  Formic Acid (970 molecules)-Test for dihedrals potentials of ryck type
TEST146 -  Formic Acid (970 molecules)-Test for dihedrals potentials of rbf type
TEST147 -  Formic Acid (970 molecules)-Test for dihedrals potentials of opls type
TEST148 -  DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for tabulated potentials for inversions
TEST149 -  DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for harmonic potentials for inversions
TEST150 -  DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of harmonic cosine (hcos) type
TEST151 -  DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of planar (plan) type
TEST152 -  DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of xpln type
TEST153 -  Calcite bulk (4320 atoms) -Test for calcite potential
TEST154 -  Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of harmonic (harm) type
TEST155 -  Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of truncated harmonic (thrm) type
TEST156 -  Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of screened harmonic (shrm) type
TEST157 -  Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of screened vessal (bvs1) type
TEST158 -  Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of truncated vessal (bvs2) type
TEST159 -  Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of H-bond (hbnd) type
TEST160 -  Chain polymer with 6 units and Carboxy Head Group-Test for tether potentials of harmonic (harm) type
TEST161 -  Chain polymer with 6 units and Carboxy Head Group-Test for tether potentials of restained harmonic (rhrm) type
TEST162 -  Chain polymer with 6 units and Carboxy Head Group-Test for tether potentials of quartic (quar) type
TEST163 -  Slab of Hf-O-Pt-Zr in vacuum - Test for external electric field, now in units of V/A
TEST164 -  Ficticious gas of Potassium Chloride - Test for external magnetic field
TEST165 -  Slab of Hf-O-Pt-Zr in vacuum - Test for an oscillatory external electric field of w=0.001 ps^-1
TEST 166 - Amorphous silicon with tersoff potential
TEST 167 - Sodium Chloride (27,000 ions) - Direct calculation of vdW interactions instead of evaluation by splining over tabulated values in memory
TEST 168 - Sodium Chloride (27,000 ions) - Use of a force-shifting procedure to vdW interactions
TEST 169 - Sodium Chloride (27,000 ions) - Direct calculation of vdW interactions (as in TEST177) plus use of a force-shifting procedure (as TEST178)

TEST170 -  NaCl Pseudo thermostat 
TEST171 -  LJ liquid with generalised lj potential ala Frenkel
TEST172 -  EVB simulation of a single reactive Malonaldehyde in gas phase, 3 Force fields
TEST173 -  EVB simulation of reactive Malonaldehyde in non-reactive water (599 water molecules) constant NPT
TEST174 -  EVB simulation of reactive Malonaldehyde in non-reactive water (599 water molecules) gaussian NVT
TEST175 -  EVB simulation of reactive Malonaldehyde in non-reactive water (599 water molecules) minimize
TEST176 -  LiF System, velocity autocorrelation with revive
TEST177 -  Argon System, stress autocorrelation
TEST178 -  Argon System, heat_flux autocorrelation
TEST179 -  Argon System, mean square displacement
TEST180 -  Argon System, radial distribution function
TEST181 -  LiF System, currents calculation
TEST182 -  many-body DPD liquid density contrast - Test for mDPD potentials 
TEST183 -  AlkylSuphate DPD test for Slater type charge smearing
TEST184 -  Rigid body positions, velocity, and orientation correlations for SF6
TEST185 -  TIP6P, test of all vdw_mix_methods for LJ
TEST186 -  Mo, OpenKIM portable model SNAP_ChenDengTran_2017_Mo__MO_698578166685_000
