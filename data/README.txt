Standard Test Cases for DL_POLY_4
----------------------------------

Because of the size of the data files for the DL_POLY_4
standard test cases, they are not shipped in the standard
download of the DL_POLY_4 source.  Instead, users are requested
to download them from the CCP5 FTP server as follows:

FTP site : ftp.dl.ac.uk
Username : anonymous
Password : your email address
Directory: ccp5/DL_POLY/DL_POLY_4.0/DATA/4.09
Files    : test_X.tar.xz

where `_X' stands for the test case number.

Remember to use the BINARY data option when transferring these
files.

Unpack the files in the `data' subdirectory using firstly
`gunzip' to uncompress them and then `tar -xf' to create the
`TEST_X' directory.

TEST 01 - Sodium Chloride (27,000 ions) interacting via Born-Huggins-Meyer potentials
TEST 02 - DMPC in Water (51,737 atoms)
TEST 03 - KNaSi2O5 - Potassium/Sodium Disilicate Glass (69,120 ions)
TEST 04 - Gramicidin A Molecules in Water (99,120 ions)
TEST 05 - SiC with Tersoff Potential (74,088 atoms)
TEST 06 - Cu3Au alloy with Sutton-Chen (metal) Potentials (32,000 atoms)
TEST 07 - Lipid Bilayer in Water (12,428 ions)
TEST 08 - MgO with Adiabatic Shell Model (8,000 ions + 4,000 shells)
TEST 09 - MgO with Relaxed Shell Model (8,000 ions + 4,000 shells)
TEST 10 - Potential of Mean Force on K+ (500 ions) in Water (13,000 ions)
TEST 11 - Cu3Au Alloy with Gupta (metal) Potentials (32,000 atoms)
TEST 12 - Cu with EAM (metal) Potential (32,000 atoms)
TEST 13 - Al with Analytic Sutton-Chen (metal) Potential (32,000 atoms)
TEST 14 - Al with EAM Sutton-Chen (metal) Potential (32,000 atoms)
TEST 15 - NiAl Alloy with EAM (metal) Potentials (27,648 atoms)
TEST 16 - Fe with Finnis-Sinclair (metal) Potential (31,250 atoms)
TEST 17 - Ni with EAM (metal) Potential (32,000 atoms)
TEST 18 - SPC IceVII Water with CBs (34,992 ions)
TEST 19 - SPC IceVII Water with RBs (34,992 ions)
TEST 20 - NaCl Molecules in SPC Water Represented as CBs+RBs (26,816 ions)
TEST 21 - TIP4P water: RBs with a Massless Charged Site (29,052 particles)
TEST 22 - Ionic Liquid Dimethylimidazolium Chloride as RBs (44,352 ions)
TEST 23 - Calcite Nano-Particles (600) in TIP3P Water (6,904 molecules, 23,712 ions)
TEST 24 - Iron/Carbon Alloy with EEAM (metal) Potentials (36,803 particles)
TEST 25 - Iron/Chromium Alloy with 2BEAM (metal) Potentials (32,000 particles)
TEST 28 - Butane in CCl4 Solution with Umbrella Sampling via PLUMED
TEST 29 - Iron with tabulated EAM (metal) Potential (54,000 atoms), Two-Temperature Model, Cascade
TEST 30 - Silicon with original Tersoff (T3) Potential (200,000 atoms), Two-Temperature Model, Swift heavy ion irradiation
TEST 31 - Tungsten with extended Finnis-Sinclair Potential (722,672 atoms), Two-Temperature Model, laser-induced relaxation of thin film
TEST 32 - Sodium Chloride (27,000 ions), NVE ensemble
TEST 33 - Sodium Chloride (27,000 ions), NVT ensemble, Evans thermostat
TEST 34 - Sodium Chloride (27,000 ions), NVT ensemble, Langevin thermostat
TEST 35 - Sodium Chloride (27,000 ions), NVT ensemble, Andersen thermostat
TEST 36 - Sodium Chloride (27,000 ions), NVT ensemble, Nose-Hoover thermostat
TEST 37 - Sodium Chloride (27,000 ions), NVT ensemble, Gentle-Stochastic thermostat
TEST 38 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of first order, no drag coefficient
TEST 39 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of first order, with drag coefficient
TEST 40 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of second order, no drag coefficient
TEST 41 - Sodium Chloride (27,000 ions), NVT ensemble, DPD thermostat of second order, with drag coefficient
TEST 42 - Sodium Chloride (27,000 ions), NPT ensemble, Langevin thermostat-barostat
TEST 43 - Sodium Chloride (27,000 ions), NPT ensemble, Berendsen thermostat-barostat
TEST 44 - Sodium Chloride (27,000 ions), NPT ensemble, Nose-Hoover thermostat-barostat
TEST 45 - Sodium Chloride (27,000 ions), NPT ensemble, Martyna-Tuckerman-Klein thermostat-barostat
TEST 46 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Langevin thermostat-barostat
TEST 47 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Berendsen thermostat-barostat
TEST 48 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Nose-Hoover thermostat-barostat
TEST 49 - Sodium Chloride (27,000 ions), NsT (NPT anisotropic ensemble), Martyna-Tuckerman-Klein thermostat-barostat
TEST 50 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Langevin thermostat /semi-isotropic barostat : constant normal pressure (Pn)
TEST 51 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Berendsen thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 52 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Nose-Hoover thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 53 - Wetted POPC membrane (32,107 ions), N-Pn-AT ensemble, Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 54 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Langevin thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 55 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Berendsen thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 56 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Nose-Hoover thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 57 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : constant normal pressure (Pn)
TEST 58 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Langevin thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST 59 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Berendsen thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST 60 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Nose-Hoover thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST 61 - Wetted POPC membrane (32,107 ions), N-Pn-gamma-T ensemble, constant surface tension (gamma), Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : constant normal pressure (Pn)/semi-orthorhombic MD cell constraints
TEST 62 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Langevin thermostat/semi-isotropic barostat    : orthorhombic MD cell constraints
TEST 63 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Berendsen thermostat/semi-isotropic barostat   : orthorhombic MD cell constraints
TEST 64 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Nose-Hoover thermostat/semi-isotropic barostat : orthorhombic MD cell constraints
TEST 65 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : orthorhombic MD cell constraints
TEST 66 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Langevin thermostat/semi-isotropic barostat    : semi-orthorhombic MD cell constraints
TEST 67 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Berendsen thermostat/semi-isotropic barostat   : semi-orthorhombic MD cell constraints
TEST 68 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Nose-Hoover thermostat/semi-isotropic barostat : semi-orthorhombic MD cell constraints
TEST 69 - Wetted POPC membrane (32,107 ions), NPT anisotropic ensemble, Martyna-Tuckerman-Klein thermostat/semi-isotropic barostat : semi-orthorhombic MD cell constraints
TEST 70 - Impact and defects analysis for rdiation damage in Sodium Chloride (27,000 ions)
TEST 71 - Sodium Chloride (27,000 ions) interacting via 12-6 potentials
TEST 72 - Sodium Chloride (27,000 ions) interacting via Lenhard-Jones potentials
TEST 73 - Sodium Chloride (27,000 ions) interacting via Lenhard-Jones cohesive potentials
TEST 74 - Sodium Chloride (27,000 ions) interacting via n-m potentials
TEST 75 - Sodium Chloride (27,000 ions) interacting via Buckingham potentials
TEST 76 - Sodium Chloride (27,000 ions) interacting via Hydrogen-bond potentials
TEST 77 - Sodium Chloride (27,000 ions) interacting via Shifted force n-m potential potentials
TEST 78 - Sodium Chloride (27,000 ions) interacting via Morse potentials
TEST 79 - Sodium Chloride (27,000 ions) interacting via Shifted Weeks-Chandler-Andersen potentials
TEST 80 - Sodium Chloride (27,000 ions) interacting via Standard DPD potentials
TEST 81 - Sodium Chloride (27,000 ions) interacting via 14-7 pair potentials
TEST 82 - Sodium Chloride (27,000 ions) interacting via Morse modified potentials
TEST 83 - Sodium Chloride (27,000 ions) interacting via Rydberg potentials
TEST 84 - Sodium Chloride (27,000 ions) interacting via Ziegler-Biersack-Littmark (ZBL) potentials
TEST 85 - Sodium Chloride (27,000 ions) interacting via ZBL mixed with Morse potentials
TEST 86 - Sodium Chloride (27,000 ions) interacting via ZBL mixed with Buckingham potentials
TEST 87 - Sodium Chloride (27,000 ions) interacting via tabulated potentials
TEST 88 - Sodium Chloride (27,000 ions) interacting via Lennard-Jones tapered with Mei-Davenport-Fernando (MDF) taper potentials
TEST 89 - Sodium Chloride (27,000 ions) interacting via Buckingham tapered with MDF potentials
TEST 90 - Sodium Chloride (27,000 ions) interacting via 12-6 tapered with MDF potentials
TEST 91 - Restart of TEST 04
TEST 92 - Formic Acid (970 molecules)-Test for bond potentials from tabulated values
TEST 93 - Formic Acid (970 molecules)-Test for bond potentials of tab+vdW+coul (-tab) type
TEST 94 - Formic Acid (970 molecules)-Test for bond potentials of harmonic (harm) type
TEST 95 - Formic Acid (970 molecules)-Test for bond potentials of harmonic+vdW+coul (-hrm) type
TEST 96 - Formic Acid (970 molecules)-Test for bond potentials of Morse (mors) type
TEST 97 - Formic Acid (970 molecules)-Test for bond potentials of Morse+vdW+coul (-mrs) type
TEST 98 - Formic Acid (970 molecules)-Test for bond potentials of 12-6 type
TEST 99 - Formic Acid (970 molecules)-Test for bond potentials of (12-6)+vdW+coul (-126) type
TEST 100 -Formic Acid (970 molecules)-Test for bond potentials of Lennard-Jones (lj) type
TEST 101 -Formic Acid (970 molecules)-Test for bond potentials of (lj)+vdW+coul (-lj) type
TEST 102 -Formic Acid (970 molecules)-Test for bond potentials of restrained harmonic (rhrm) type
TEST 103 -Formic Acid (970 molecules)-Test for bond potentials of restricted harmonic+vdW+coul (-rhm) type
TEST 104 -Formic Acid (970 molecules)-Test for bond potentials of quar type
TEST 105 -Formic Acid (970 molecules)-Test for bond potentials of quar+vdW+coul (-qur) type
TEST 106 -Formic Acid (970 molecules)-Test for bond potentials of Buckingham (buck) type
TEST 107 -Formic Acid (970 molecules)-Test for bond potentials of Buckingham+vdW+coul (-bck) type
TEST 108 -Formic Acid (970 molecules)-Test for bond potentials of Coulomb (coul) type
TEST 109 -Formic Acid (970 molecules)-Test for bond potentials of Electroctatic+vdW+coul (-cul) type
TEST 110 -Formic Acid (970 molecules)-Test for bond potentials of fene type
TEST 111 -Formic Acid (970 molecules)-Test for bond potentials of FENE+vdW+coul (-fne) type
TEST 112 -Formic Acid (970 molecules)-Test for bond potentials of mmts type
TEST 113 -Formic Acid (970 molecules)-Test for bond potentials of mmst+vdW+coul (-mst) type
TEST 114 -Formic Acid (970 molecules)-Test for bond tabulated (tab) potentials
TEST 115 -Formic Acid (970 molecules)-Test for tabulated angle potentials of tab+coul+vdW (-tab) type
TEST 116 -Formic Acid (970 molecules)-Test for angle potentials of harmonic (harm) type
TEST 117 -Formic Acid (970 molecules)-Test for angle potentials of harmonic+coul+vdW (-hrm) type
TEST 118 -Formic Acid (970 molecules)-Test for angle potentials of quartic (quartic) order
TEST 119 -Formic Acid (970 molecules)-Test for angle potentials of quartic+coul+vdW (-qur) type
TEST 120 -Formic Acid (970 molecules)-Test for angle potentials of Truncated harmonic (thrm) type
TEST 121 -Formic Acid (970 molecules)-Test for angle potentials of Truncate harmonic+coul+vdW (-thm) type
TEST 122 -Formic Acid (970 molecules)-Test for angle potentials of Screened harmonic (shrm) type
TEST 123 -Formic Acid (970 molecules)-Test for angle potentials of Screnned harmonic+coul+vdW (-shm) type
TEST 124 -Formic Acid (970 molecules)-Test for angle potentials of Screened Vessal (bvs1) type
TEST 125 -Formic Acid (970 molecules)-Test for angle potentials of Screnned Vessal+coul+vdW (-bv1) type
TEST 126 -Formic Acid (970 molecules)-Test for angle potentials of Truncated Vessal (bvs2) type
TEST 127 -Formic Acid (970 molecules)-Test for angle potentials of Truncate Vessal+coul+vdW (-bv2) type
TEST 128 -Formic Acid (970 molecules)-Test for angle potentials of harmonic cosine (hacos) type
TEST 129 -Formic Acid (970 molecules)-Test for angle potentials of harmonic cosine+coul+vdW (-hcs) type
TEST 130 -Formic Acid (970 molecules)-Test for angle potentials of cosine (cos) type
TEST 131 -Formic Acid (970 molecules)-Test for angle potentials of cosine+coul+vdW (-hcs) type
TEST 132 -Formic Acid (970 molecules)-Test for angle potentials of MM3 stretch-bend (mmsb) type
TEST 133 -Formic Acid (970 molecules)-Test for angle potentials of MM3 stretch-bend+coul+vdW (-msb) type
TEST 134 -Formic Acid (970 molecules)-Test for angle potentials of compass stretch-stretch (stst) type
TEST 135 -Formic Acid (970 molecules)-Test for angle potentials of Compass stretch-stretch +coul+vdW (-sts) type
TEST 136 -Formic Acid (970 molecules)-Test for angle potentials of compass stretch-bend (stbe) type
TEST 137 -Formic Acid (970 molecules)-Test for angle potentials of Compass stretch-bend +coul+vdW (-stb) type
TEST 138 -Formic Acid (970 molecules)-Test for angle potentials of compass all-terms (cmps) type
TEST 139 -Formic Acid (970 molecules)-Test for angle potentials of Compass all +coul+vdW (-cmp) type
TEST 140 -Formic Acid (970 molecules)-Test for angle potentials of mmbd type
TEST 141 -Formic Acid (970 molecules)-Test for angle potentials of mmbd+coul+vdW (-mbd) type
TEST 142 -Formic Acid (970 molecules)-Test for angle potentials of KKY type
TEST 143 -Formic Acid (970 molecules)-Test for angle potentials of KKY +coul+vdW (-kky) type
TEST 144 -Formic Acid (970 molecules)-Test for tabulated potentials for dihedrals
TEST 145 -Formic Acid (970 molecules)-Test for cosine potentials for dihedrals
TEST 146 -Formic Acid (970 molecules)-Test for dihedrals potentials of harmonic (harm) type
TEST 147 -Formic Acid (970 molecules)-Test for dihedrals of harmonic cosine (hcos) type
TEST 148 -Formic Acid (970 molecules)-Test for dihedrals of harmonic cosine (cos3) type
TEST 149 -Formic Acid (970 molecules)-Test for dihedrals potentials of ryck type
TEST 150 -Formic Acid (970 molecules)-Test for dihedrals potentials of rbf type
TEST 151 -Formic Acid (970 molecules)-Test for dihedrals potentials of opls type
TEST 152 -DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for tabulated potentials for inversions
TEST 153 -DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for harmonic potentials for inversions
TEST 154 -DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of harmonic cosine (hcos) type
TEST 155 -DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of planar (plan) type
TEST 156 -DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of xpln type
TEST 157 -Calcite bulk (4320 atoms) -Test for calcite potential
TEST 158 -Sauer model sans shell + Benzene-Test for three-body potentials of harmonic (harm) type
TEST 159 -Sauer model sans shell + Benzene-Test for three-body potentials of truncated harmonic (thrm) type
TEST 160 -Sauer model sans shell + Benzene-Test for three-body potentials of screened harmonic (shrm) type
TEST 161 -Sauer model sans shell + Benzene-Test for three-body potentials of screened vessal (bvs1) type
TEST 162 -Sauer model sans shell + Benzene-Test for three-body potentials of truncated vessal (bvs2) type
TEST 163 -Sauer model sans shell + Benzene-Test for three-body potentials of H-bond (hbnd) type

TEST 169 -Chain polymer with   6 units and Carboxy Head Group-Test for tether potentials of harmonic (harm) type
TEST 170 -Chain polymer with   6 units and Carboxy Head Group-Test for tether potentials of restained harmonic (rhrm) type
TEST 171 -Chain polymer with   6 units and Carboxy Head Group-Test for tether potentials of quartic (quar) type

TEST 173 - Slab of Hf-O-Pt-Zr in vacuum (55872 ions) - Test for external electric field, now in units of V/A
TEST 174 - Ficticious gas of Potassium Chloride at 1000K (only 1372 ions) - Test for external magnetic field, now in units of Tesla
TEST 175 - Slab of Hf-O-Pt-Zr in vacuum (55872 ions) - Test for external oscillatory electric field, now in units of V/A
TEST 176 - Amorphous silicon with tersoff potential
TEST 177 - Sodium Chloride (27,000 ions) - Direct calculation of vdW interactions instead of evaluation by splining over tabulated values in memory
TEST 178 - Sodium Chloride (27,000 ions) - Use of a force-shifting procedure to vdW interactions
TEST 179 - Sodium Chloride (27,000 ions) - Direct calculation of vdW interactions (as in TEST177) plus use of a force-shifting procedure (as TEST178)

TEST 181 - Pseudo thermostat
