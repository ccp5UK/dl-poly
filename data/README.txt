Standard Test Cases for DL_POLY_4
----------------------------------

Because of the size of the data files for the DL_POLY_4
standard test cases, they are not shipped in the standard
download of the DL_POLY_4 source. Test input files are downloaded
automatically when building/ running CMake with the CMake variable BUILD_TESTING=ON.
Reference and configuration data for tests is stored in the reference folder.


This can be done as follows:

folder="build-mpi-testing"
rm -rf $folder && mkdir $folder && pushd $folder
cmake ../ -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON 

Unpack the files in the `data' subdirectory using firstly
`gunzip' to uncompress them and then `tar -xf' to create the
`TEST_X' directory.

Version controlled test cases may be downloaded from the ccp5UK 
Github: https://github.com/ccp5UK/dl-poly-data

The following is a list of DL_POLY_4 test cases currently used:

TEST001 SODIUM CHLORIDE WITH (27000 IONS)
TEST002 DMPC in WATER (51737 ATOMS)
TEST003 SODIUM/POTASSIUM DISILICATE GLASS WITH 3 BODY FORCES (69120 IONS)
TEST004 8xGRAMICIDIN A WITH WATER SOLVATING (99120 ATOMS)
TEST005 SiC with Tersoff Potential (74088 ATOMS)
TEST006 Cu3Au ALLOY (32000 ATOMS)
TEST007 DPPC in water
TEST008 Shell Model MgO
TEST009 Shell Model MgO
TEST010 Potential of mean force on K+ in water
TEST011 Cu3Au ALLOY (32000 ATOMS)
TEST012 COPPER WITH TABULATED EAM POTENTIAL
TEST013 ALUMINIUM WITH SUTTON-CHEN POTENTIAL
TEST014 ALUMINIUM WITH TABULATED EAM POTENTIAL
TEST015 ALUMINIUM/NICKEL ALLOY WITH TABULATED EAM POTENTIAL
TEST016 DL_POLY TEST CASE: Iron with Finnis Sinclair Potential
TEST017 NICKEL WITH TABULATED EAM POTENTIAL
TEST018 SPC IceVII Minimisation
TEST019 SPC IceVII Minimisation
TEST020 SODIUM CHLORIDE MOLECULES IN SPC WATER(CB+RB)
TEST021 TIP4P WATER
TEST022 dmimcl rigid
TEST023 Calcite Nanoparticles in Water
TEST024 IRON/CARBON ALLOY WITH TABULATED EEAM POTENTIAL
TEST025 FeCr WITH TABULATED EAM POTENTIAL
TEST026 Butane
TEST027 IRON WITH TABULATED EAM POTENTIAL, TTM AND CASCADE
TEST028 SWIFT HEAVY ION IRRADIATION OF SILICON WITH TERSOFF T3 POTENTIAL
TEST029 RELAXATION OF A TUNGSTEN THIN FILM
TEST030 SODIUM CHLORIDE WITH (27000 IONS)
TEST031 SODIUM CHLORIDE WITH (27000 IONS)
TEST032 SODIUM CHLORIDE WITH (27000 IONS)
TEST033 SODIUM CHLORIDE WITH (27000 IONS)
TEST034 SODIUM CHLORIDE WITH (27000 IONS)
TEST035 SODIUM CHLORIDE WITH (27000 IONS)
TEST036 SODIUM CHLORIDE WITH (27000 IONS)
TEST037 SODIUM CHLORIDE WITH (27000 IONS)
TEST038 SODIUM CHLORIDE WITH (27000 IONS)
TEST039 SODIUM CHLORIDE WITH (27000 IONS)
TEST040 SODIUM CHLORIDE WITH (27000 IONS)
TEST041 SODIUM CHLORIDE WITH (27000 IONS)
TEST042 SODIUM CHLORIDE WITH (27000 IONS)
TEST043 SODIUM CHLORIDE WITH (27000 IONS)
TEST044 SODIUM CHLORIDE WITH (27000 IONS)
TEST045 SODIUM CHLORIDE WITH (27000 IONS)
TEST046 SODIUM CHLORIDE WITH (27000 IONS)
TEST047 SODIUM CHLORIDE WITH (27000 IONS)
TEST048 POPC membrane in contact with water (32107 IONS)
TEST049 POPC membrane in contact with water (32107 IONS)
TEST050 POPC membrane in contact with water (32107 IONS)
TEST051 POPC membrane in contact with water (32107 IONS)
TEST052 POPC membrane in contact with water (32107 IONS)
TEST053 POPC membrane in contact with water (32107 IONS)
TEST054 POPC membrane in contact with water (32107 IONS)
TEST055 POPC membrane in contact with water (32107 IONS)
TEST056 POPC membrane in contact with water (32107 IONS)
TEST057 POPC membrane in contact with water (32107 IONS)
TEST058 POPC membrane in contact with water (32107 IONS)
TEST059 POPC membrane in contact with water (32107 IONS)
TEST060 POPC membrane in contact with water (32107 IONS)
TEST061 POPC membrane in contact with water (32107 IONS)
TEST062 POPC membrane in contact with water (32107 IONS)
TEST063 POPC membrane in contact with water (32107 IONS)
TEST064 POPC membrane in contact with water (32107 IONS)
TEST065 POPC membrane in contact with water (32107 IONS)
TEST066 POPC membrane in contact with water (32107 IONS)
TEST067 POPC membrane in contact with water (32107 IONS)
TEST068 SODIUM CHLORIDE WITH (27000 IONS)
TEST069 SODIUM CHLORIDE WITH (27000 IONS)
TEST070 SODIUM CHLORIDE WITH (27000 IONS)
TEST071 SODIUM CHLORIDE WITH (27000 IONS)
TEST072 SODIUM CHLORIDE WITH (27000 IONS)
TEST073 SODIUM CHLORIDE WITH (27000 IONS)
TEST074 SODIUM CHLORIDE WITH (27000 IONS)
TEST075 SODIUM CHLORIDE WITH (27000 IONS)
TEST076 SODIUM CHLORIDE WITH (27000 IONS)
TEST077 SODIUM CHLORIDE WITH (27000 IONS)
TEST078 SODIUM CHLORIDE WITH (27000 IONS)
TEST079 SODIUM CHLORIDE WITH (27000 IONS)
TEST080 SODIUM CHLORIDE WITH (27000 IONS)
TEST081 SODIUM CHLORIDE WITH (27000 IONS)
TEST082 SODIUM CHLORIDE WITH (27000 IONS)
TEST083 SODIUM CHLORIDE WITH (27000 IONS)
TEST084 SODIUM CHLORIDE WITH (27000 IONS)
TEST085 SODIUM CHLORIDE WITH (27000 IONS)
TEST086 SODIUM CHLORIDE WITH (27000 IONS)
TEST087 SODIUM CHLORIDE WITH (27000 IONS)
TEST088 SODIUM CHLORIDE WITH (27000 IONS)
TEST089 8xGRAMICIDIN A WITH WATER SOLVATING (99120 ATOMS)
TEST090 Formic Acid (970 molecules)-Test for bond potentials from tabulated values
TEST091 Formic Acid (970 molecules)-Test for bond potentials of tab+vdW+coul (-tab) type
TEST092 Formic Acid (970 molecules)-Test for bond potentials of harmonic (harm) type
TEST093 Formic Acid (970 molecules)-Test for bond potentials of harmonic+vdW+coul (-hrm) type
TEST094 Formic Acid (970 molecules)-Test for bond potentials of Morse (mors) type
TEST095 Formic Acid (970 molecules)-Test for bond potentials of Morse+vdW+coul (-mrs) type
TEST096 Formic Acid (970 molecules)-Test for bond potentials of 12-6 type
TEST097 Formic Acid (970 molecules)-Test for bond potentials of (12-6)+vdW+coul (-126) type
TEST098 Formic Acid (970 molecules)-Test for bond potentials of Lennard-Jones (lj) type
TEST099 Formic Acid (970 molecules)-Test for bond potentials of (lj)+vdW+coul (-lj) type
TEST100 Formic Acid (970 molecules)-Test for bond potentials of restrained harmonic (rhrm) type
TEST101 Formic Acid (970 molecules)-Test for bond potentials of restricted harmonic+vdW+coul (-rhm) type
TEST102 Formic Acid (970 molecules)-Test for bond potentials of quar type
TEST103 Formic Acid (970 molecules)-Test for bond potentials of quar+vdW+coul (-qur) type
TEST104 Formic Acid (970 molecules)-Test for bond potentials of Buckingham (buck) type
TEST105 Formic Acid (970 molecules)-Test for bond potentials of Buckingham+vdW+coul (-bck) type
TEST106 Formic Acid (970 molecules)-Test for bond potentials of Coulomb (coul) type
TEST107 Formic Acid (970 molecules)-Test for bond potentials of Electroctatic+vdW+coul (-cul) type
TEST108 Formic Acid (970 molecules)-Test for bond potentials of fene type
TEST109 Formic Acid (970 molecules)-Test for bond potentials of FENE+vdW+coul (-fne) type
TEST110 Formic Acid (970 molecules)-Test for bond potentials of mmts type
TEST111 Formic Acid (970 molecules)-Test for bond potentials of mmst+vdW+coul (-mst) type
TEST112 Formic Acid (970 molecules)-Test for bond tabulated (tab) potentials  
TEST113 Formic Acid (970 molecules)-Test for tabulated angle potentials of tab+coul+vdW (-tab) type
TEST114 Formic Acid (970 molecules)-Test for angle potentials of harmonic (harm) type
TEST115 Formic Acid (970 molecules)-Test for angle potentials of harmonic+coul+vdW (-hrm) type
TEST116 Formic Acid (970 molecules)-Test for angle potentials of quartic (quartic) order
TEST117 Formic Acid (970 molecules)-Test for angle potentials of quartic+coul+vdW (-qur) type
TEST118 Formic Acid (970 molecules)-Test for angle potentials of Truncated harmonic (thrm) type 
TEST119 Formic Acid (970 molecules)-Test for angle potentials of Truncate harmonic+coul+vdW (-thm) type
TEST120 Formic Acid (970 molecules)-Test for angle potentials of Screened harmonic (shrm) type
TEST121 Formic Acid (970 molecules)-Test for angle potentials of Screnned harmonic+coul+vdW (-shm) type
TEST122 Formic Acid (970 molecules)-Test for angle potentials of Screened Vessal (bvs1) type 
TEST123 Formic Acid (970 molecules)-Test for angle potentials of Screnned Vessal+coul+vdW (-bv1) type
TEST124 Formic Acid (970 molecules)-Test for angle potentials of Truncated Vessal (bvs2) type 
TEST125 Formic Acid (970 molecules)-Test for angle potentials of Truncate Vessal+coul+vdW (-bv2) type
TEST126 Formic Acid (970 molecules)-Test for angle potentials of harmonic cosine (hacos) type
TEST127 Formic Acid (970 molecules)-Test for angle potentials of harmonic cosine+coul+vdW (-hcs) type
TEST128 Formic Acid (970 molecules)-Test for angle potentials of cosine (cos) type
TEST129 Formic Acid (970 molecules)-Test for angle potentials of cosine+coul+vdW (-hcs) type
TEST130 Formic Acid (970 molecules)-Test for angle potentials of MM3 stretch-bend (mmsb) type 
TEST131 Formic Acid (970 molecules)-Test for angle potentials of MM3 stretch-bend+coul+vdW (-msb) type
TEST132 Formic Acid (970 molecules)-Test for angle potentials of compass stretch-stretch (stst) type 
TEST133 Formic Acid (970 molecules)-Test for angle potentials of Compass stretch-stretch +coul+vdW (-sts) type
TEST134 Formic Acid (970 molecules)-Test for angle potentials of compass stretch-bend (stbe) type 
TEST135 Formic Acid (970 molecules)-Test for angle potentials of Compass stretch-bend +coul+vdW (-stb) type
TEST136 Formic Acid (970 molecules)-Test for angle potentials of compass all-terms (cmps) type 
TEST137 Formic Acid (970 molecules)-Test for angle potentials of Compass all +coul+vdW (-cmp) type
TEST138 Formic Acid (970 molecules)-Test for angle potentials of mmbd type
TEST139 Formic Acid (970 molecules)-Test for angle potentials of mmbd+coul+vdW (-mbd) type
TEST140 Formic Acid (970 molecules)-Test for tabulated potentials for dihedrals
TEST141 Formic Acid (970 molecules)-Test for cosine potentials for dihedrals
TEST142 Formic Acid (970 molecules)-Test for dihedrals potentials of harmonic (harm) type
TEST143 Formic Acid (970 molecules)-Test for dihedrals of harmonic cosine (hcos) type
TEST144 Formic Acid (970 molecules)-Test for dihedrals of harmonic cosine (cos3) type
TEST145 Formic Acid (970 molecules)-Test for dihedrals potentials of ryck type
TEST146 Formic Acid (970 molecules)-Test for dihedrals potentials of rbf type
TEST147 Formic Acid (970 molecules)-Test for dihedrals potentials of opls type
TEST148 DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for tabulated potentials for inversions 
TEST149 DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for harmonic potentials for inversions 
TEST150 DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of harmonic cosine (hcos) type
TEST151 DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of planar (plan) type
TEST152 DPPC (200 molecules) in Water (32626 molecules) Lipid Bilayer -Test for inversions potentials of xpln type
TEST153 Calcite bulk (4320 atoms) -Test for calcite potential
TEST154 Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of harmonic (harm) type                        
TEST155 Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of truncated harmonic (thrm) type                        
TEST156 Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of screened harmonic (shrm) type                        
TEST157 Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of screened vessal (bvs1) type                        
TEST158 Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of truncated vessal (bvs2) type                        
TEST159 Silicalite: Sauer model sans shell + Benzene-Test for three-body potentials of H-bond (hbnd) type                        
TEST160 Chain polymer with   6 units and Carboxy Head Group-Test for tether potentials of harmonic (harm) type                        
TEST161 Chain polymer with   6 units and Carboxy Head Group-Test for tether potentials of restained harmonic (rhrm) type                        
TEST162 Chain polymer with   6 units and Carboxy Head Group-Test for tether potentials of quartic (quar) type                        
TEST163 Slab of Hf-O-Pt-Zr in vacuum - Test for external electric field, now in units of V/A 
TEST164 Ficticious gas of Potassium Chloride - Test for external magnetic field
TEST165 Slab of Hf-O-Pt-Zr in vacuum - Test for an oscillatory external electric field of w=0.001 ps^-1
TEST166 SWIFT HEAVY ION IRRADIATION OF SILICON WITH TERSOFF T3 POTENTIAL
TEST167 SODIUM CHLORIDE WITH (27000 IONS)
TEST168 SODIUM CHLORIDE WITH (27000 IONS)
TEST169 SODIUM CHLORIDE WITH (27000 IONS)
TEST170 SODIUM CHLORIDE WITH (27000 IONS)
TEST171 DL_POLY: ar
TEST172 EVB simulation of a single reactive Malonaldehyde in gas phase
TEST173 EVB simulation of reactive Malonaldehyde in non-reactive water (599 water molecules)
TEST174 EVB simulation of reactive Malonaldehyde in non-reactive water (599 water molecules)
TEST175 EVB simulation of reactive Malonaldehyde in non-reactive water (599 water molecules)
TEST176 Argon, Velocity auto correlation function
TEST177 Argon, Stress auto correlation function
TEST178 Argon, heat flux auto correlation function
TEST179 Argon, MSD
TEST180 Argon, RDF
