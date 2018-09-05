Release notes for DL_POLY_4
===========================

Version 4.09: September 2018
----------------------------

+ New features, functionality, (re)implementations, etc.
  + new VDW potentials
    + modified morse 
    + Ziegler-Biersack-Littmark, ZBL, ZBL mixed with Morse and Buckingham
    + Mei-Devenport-Fernando taper for Lennard-Jones, 12-6 and Buckingham
    + Rydberg
    + modified LJ
  + new RDF error bars calculation functionality
  + new optional I/O file-naming functionality
    + CONTROL file can be passed as command line argument under any name
    + main input/output filenames can be changed via directives in CONTROL file
  + two-temperature model (TTM)
    + electronic temperature evolution for metal and non-metal systems
    + initial electronic energy depositions for lasers and Swift heavy ion irradiation
    + inhomogeneous Langevin thermostat (also available without TTM)
  + new shell polarisation handling on multipolar level, build\_tplg\_intra - extended
  + new CHARMM shells enabled (buildi\_chrm\_intra, coul\_chrm\_forces)
  + new zero\_k smoothing performance handling with a new frequency option
  + new extra reporting recommendations (and fixes) on densvar for passing in read\_config\*, build\_book\_intra and relocate\_particles
  + new, better exclusion array length scanning & reporting
  + new, better sanity check for Ewald conditions
  + new, better sanity check for array sizes non-PBC cells
  + new implementation of dry run 0K velocity field and stress reporting
  + new umbrella sampling restraint external field and RDFs
  + new quaternions unit measure rescale for VV
  + new scale\_temperature frozen rotational degrees exception
  + new set\_temperature and scale\_temperature exception handling for 0K field for non-0K point restart
  + new faster shell filtering in all routines
  + new zero fire option, optional CGM stepping
  + new INSTALL and CITATION.cff

+ New and changed defaults
  + default changed for l\_dis - applied at every step instead at step 0 only
  + default changed for VNL application - enabled for all runs if no explicit zero rpad/delr
  + default change for core\_shell\_on\_top - apply to new non-restarted models only
  + default changed for core-shell distances - excluded when l\_dis is applied

+ Improvements, updates, etc.
  + improved in conditional, adaptive VNL setup and checking
  + improved core\_shell\_relax and minimise\_relax information extension and interplay
  + improved metal\_forces handling
  + improved replay\_histor\* (EoF for handling last frame preservation), dl\_poly and read\_history
  + improved (doubled) number of tests
  + improved STATIS writing (thanks to Mike Pacey, University of Lancaster,UK)
  + improved windows cross-compilation from linux, both serial and parallel
  + updated README.txt
  + improved build/Makefile\* to make SERIAL build easier
  + improved, cleaner NL/EL characters parsing
  + updated energy unit constants in setup\_module
  + improved end of frequency options reporting in w\_statistsics\_report
  + updated MANUAL - two-temperature model chapter, outline surface thickness/energy; core-shell rules and references; tabulated interactions; input files information fresh up for new and old options; minimise information; MPOLES information and notation wrt CHARMM functionality; full CHARMM druder consequences description.

+ Bug-fixes:
  + HISTORY writing info line when using 'io write master' in trajectory\_write
  + PDF sampling setting in CONTROL (scan\_control and read\_control routines)
  + restart for velocity autocorrelation functions
  + mixing potential reading
  + elec bonds
  + xpiston, zres and zrs+/- external fields
  + molecule COM calculation in kinetics\_module
  + two\_body\_forces for placement of frozen list of interactions
  + reading from tables intra\_table\_read routines, thanks to Tom Potter @ University of Durham (UK)
  + TABDIH zero element reading in dihedrals\_forces
  + TABBND force calculation in bonds\_forces
  + initialisation of nst\_b\* semi-isotropic integrators
  + scaling factors in nst\_h\*\_scl routines
  + mmst potential in bonds\_forces
  + non-isotropic system expansion of cubic boxes in system\_expan\*
  + units for TABLE file

Version 4.08: March 2016
------------------------

+ new compilation framework with CMake (old legacy Makefiles kept)
+ implementation of multipolar electrostatics interactions, up to hexadecapole moment (order 4)
+ PLUMED integrability
+ new subcelling keyword to improve performance of VNL build-up
+ improved handling of array bounds
+ improved handling of image condition indicator
+ improved handling of per-run counters for constraints and CVNL skips
+ improved coordinates DD passing, wrapping up only when crossing MD cell boundaries
+ improved RDF accuracy for degenerated systems when averaging over very few species
+ removal of images calls for all non-bonded calculations
+ improved handling of residual halo for VNL build-up
+ case sensitivity of atom labels
+ capture of <PV> and extended degree of freedom energy (cons) in STATIS and OUTPUT
+ improved instantaneous per-atom MSD calculation in MSDTMP for NP/sT ensembles
+ new connectivity restore procedure for replay HISTORY options
+ new failsafes for dynamics calculations in replay HISTORY
+ long range corrections for 14-7 buffered VDW
+ general code clean-up, refactoring and restructuring
+ extra metadata collection at build and runtime reported
+ README, Makefiles, CMake, macros and user manual updates
+ fixes to:
  + OpenKIM integrability routines
  + VNL subcelling safety for particle on grid line
  + global safety for CVNL when no particles occupy a domain
  + global safety for SHAKE (CB and PMF) when no constraints present on a domain
  + correct density calculation for extended Finnis-Sinclair (metal) potential
  + deport\_atomic\_data re-stacking safety
  + replay HISTORY for generation of MSDTMP
  + setting halo size w.r.t. Ewald requirements
  + 3D diffusion averaging estimate
  + parallel reading without MPI-IO of very large CONFIGs
  + 14-7 buffered VDW definition
  + MM3 bond stretch potential definition

Version 4.07: January 2015
--------------------------

+ 'no strict' CONTROL option reset of rcut now depends on RBs' existence
+ improved SHAKE and RATTLE/QUENCH statistics and optimisation, improved RATTLE/QUENCH convergence criteria
+ improved CVNL statistics, exceptions handling and error reporting
+ improved CGM minimisation stability w.r.t. processor count
+ new torsional and inversional restraints
+ NVT DPD thermostat with enhanced defaults, potential mixing and integration with Shardlow's VV first and second order splittings
+ 'densvar' treatment and application to link\_cell algorithms generalised by removing 3D link-cell size limit for non-periodic models (imcon=0 or imcon=6), with exceeded limit error check and condition removed where applicable
+ 'densvar' application extended to affect mxexcl (exclusion list length) and mxfdih (dihedrals list length) values
+ mxdih default value increased
+ optimisation of two\_body\_forces vector length usage and deport\_atomic\_data reordering (re-stacking)
+ improved reporting for build\_book\_intra
+ improved handling for ewald\_real zero distance exception
+ introduced handling of NVT Langevin impulse VV implementation when entering algorithm's undefined solution zone by timestep size adjustment
+ redefinition of MBPC metal-like interaction
+ OpenKIM integration, thanks to Henry Boateng at STFC DL (UK) and Ryan Elliott at University of Minnesota (USA)
+ general code clean-up, refactoring and restructuring
+ DL\_POLY GUI code moved to Java 1.7 by default (Java 1.6 now obsolete)
+ README, Makefiles, macros and user manual updates
+ fixes to:
  + scale\_temperature RB COM velocity generation, thanks to Laurence Ellison at STFC DL (UK)
  + revert to use of integers in n-m powers types of VDW/metal potentials to avoid underflow-driven cumulative numerical inaccuracy
  + metal density table key calculation for FS-type of metal potentials when 'metal direct' CONTROL option specified
  + optimise two\_body\_forces and deport\_atomic\_data routines, thanks to Alin Elena at ICHEC (Ireland) and Victor Gamayunov at INTEL (UK)

Version 4.06.1: September 2014
------------------------------

+ implementation of many-body perturbation component potential for actinide oxides (metal potentials), improvements to metal routines
+ inclusion of run statistics reporting on bond and PMF constraints, VNL conditional update skips in OUTPUT
+ backwards compatibility for DL\_POLY\_Classic's delr option with DL\_POLY\_4's rpad (VNL conditional update)
+ relaxation of efficiency mode warning, advice given on largest numbers of nodes for balanced performance and possible simulation
+ new test cases for coarse-grained functionality
+ fixes to:
  + link\_cell\_pairs for large cutoffs leading to self-haloing
  + change delth\_max in setup\_module
  + change defaults for intramolecular PDF grids and smoothing sizes, and PDF printing
  + rename intermolecular PDF analysis output files (based on RDFs) to VDWPMF/VDWTAB
  + arrays bound for erfc+ grid for Ewald related electrostatics and zero point added to handle Coulomb collapse
  + VNL conditional set-up auto-switch with 'no strict' option in CONTROL and skip run statistics in OUTPUT
  + routines dealing with halo and halo-driven transmission data for cases with one link cell per domain
  + contention between SPME and LC/DD requirements on halo width with one link cell per domain (if k-vector size too small to accommodate maximum off-domain particles' blurring due to preset padding cutoff)
  + incorrect reading of collection intervals for intramolecular PDF analysis options in CONTROL
  + failure to collect data for PDF analysis in real run
  + extend spacing for large simulation times and/or large number of timesteps in OUTPUT
  + Tersoff KIHS potential: grids extended to add zero point for easy particle collapse handling with force capping and/or zero temperature optimisation
  + rules for 1-4 overlapping dihedral interactions with different scaling factors, to use first non-zero scaling factors from any dihedral of same 1-4 kind rule over all remaining ones (which are zeroed)
+ user manual updates

Version 4.06: August 2014
-------------------------

+ optional quartic term for core-shell interaction
+ full feedback on breaking of intra-molecular and like array sizes in build\_book\_intra
+ full estimates of array bound violations in export/deport routines
+ full system and topology matching when 'replay history' options triggered
  + new restarting capabilities
  + correct statistics for calculating dynamical properties if HISTORY frames are equi-temporal
+ improved coordinate reading
+ improved fall-back parallel reading
+ improved handling of extended NsT LFV integrators
+ improved iterative coupling of thermostat and barostat for LFV implementations of NPT Nose-Hoover, NsT Nose-Hoover and MTK integrators
+ 'impact' not affected by 'no vom' enforcement
+ Verlet neighbour list conditional update functionality by defining system cutoff padding in CONTROL
  + speed-ups from 5% to 120% depending on force-field complexity and user choice of rpad
  + also triggered with 'no strict' option
+ four new external fields: zres, zrs+, zrs- and osel
+ new potentials
  + AMOEBA-specific VDW force-field (buffered 14-7), bond and angle potential interactions
  + DPD soft repulsive (Groot-Warren) VDW interaction potential
+ five new mixing rule types for Lennard-Jones' like interactions for mixed component interactions
  + Lorentz-Berthelot
  + Fender-Halsey
  + Hogervorst (good hope)
  + Halgren HHG
  + Kong
+ failure reports and suggested densvar to pass a stage
+ speed improvements of SHAKE (divisions) and linked cells (reducing numbers of loops and IF conditions)
+ increasing array sizes for dihedrals and inversions to reduce densvar expansion for overloaded dihedrals/inversions in bio-chemical systems
+ inclusion of frozen-frozen contributions to RDFs
+ introduction of SARU random number generator to reduce dependence on DD for velocity field generation, Langevin, Andersen and GST thermostats, thanks to Michael Seaton at STFC DL (UK)
+ error check on intramolecular sites being within molecule
+ improved handling of rigid bodies and reporting algorithmic failures
+ inclusion of instructions for deploying Breathing Shell Model
+ systematic coarse-graining features: intramolecular potentials via table files and intramolecular potentials' distribution function analysis, thanks to Andrey Brukhno at University of Bath (UK)
+ velocity autocorrelation function collection and output feature, thanks to Michael Seaton at STFC DL (UK)
+ 'nfold' option check on system contiguity changed to depend on two cutoffs' criteria
+ reporting of run statistics on conditional VNL skips, configurational and RSM CGM minimisation cycles
+ user manual improvements, thanks to Michael Seaton at STFC DL (UK)
+ fixes to:
  + force-shifted n-m potential, thanks to Hari Yadav at IIT Delhi (India)
  + bug in Langevin random forces by removing COM motion constraint to avoid artificial thermalisation effects for small systems with disparate masses, thanks to Jake Stinson at UCL (UK)
  + correct bin numbering for RDFs and their contribution sums in 'replay history' option, thanks to Andrey Brukhno at University of Bath (UK)
  + correct handling of stress tensor's conjugate variable for extended NsT VV integrators, thanks to Mohammad Sedghi at University of Wyoming (USA)
  + correct buffer sizes for export and reporting exceeded buffers
  + correct handling of dihedrals with core-shell units on their 1-4 members
  + improve speed of SHAKE, thanks to Alin Elena at ICHEC (Ireland)
  + correct definition of x-piston external field handling, thanks to Andrey Brukhno at University of Bath (UK)
  + refresh of 'mock' eta for NPT VV integrators to properly update particle tracking in xscale
  + RB consistency and contiguity checks
+ discontinuation of CUDA port

Version 4.05.1: September 2013
------------------------------

+ big fix in force-shifting of van der Waals (vdw\_forces), thanks to Joost van den Ende at Radboud University Nijmegen (Netherlands)
+ optimisation of cache reuse in link\_cell\_pairs, thanks to CRAY engineers
+ CUDA port update, thanks to Buket Gursoy at ICHEC (Ireland)
+ user manual updates

Version 4.05: July 2013
-----------------------

+ new two band (2B) EAM and EEAM metal potentials, thanks to Ben Palmer at University of Birmingham (UK)
+ new boundary thermostat 'gauss' (pseudo)
+ 'no vom' option to stop removing centre-of-mass motion during integration and velocity rescaling
+ new Xpiston external force field
+ 'replay history force' option, using HISTROF positions instead of integration and calculate forces
+ core-shell interaction exclusions now apply for 1-4 dihedral interaction electrostatics scaling
+ special CONTROL options given in manual:
  + 'l\_scr' - redirects OUTPUT contents to OS default output channel (screen, interactive)
  + 'l\_fast' - abandons global safety checks, speeds up simulations with too few link cells per domain, but uncontrolled termination in case of parallel failure
  + 'l\_rout' - writes REVIVE in ASCII
  + 'l\_rin' - reads REVOLD in ASCII
  + 'l\_his' - generates one-frame HISTORY from CONFIG and terminates immediately
  + 'l\_tor' - abandons production of REVIVE and REVCON
  + 'l\_dis' - checks and reports minimum separation distance between all Verlet neighbour list pairs at (re)start
+ default array and buffer size changed to minimise memory footprint
+ parallel reading default buffer changed to increase speed, thanks to Asmina Maniopoulou at NAG Ltd. (UK)
+ extensive code clean up and refactoring
+ user manual updates
+ fixes to:
  + correct NVT Langevin VV integrator, thanks to Jake Stinson at UCL (UK)
  + boundary thermostats (pseudo), thanks to Eva Zarkadoula at QMUL (UK)

Version 4.04.2: January 2013
----------------------------

+ changes to scan\_field for maximum extents of site numbers and grid calculations
+ user manual updates
+ file format changes from UTF-8 to ASCII

Version 4.04.1: December 2012
-----------------------------

+ CUDA port update, thanks to Michael Lysaght at ICHEC (Dublin, Ireland)
+ better estimates of list sizes for exclusions and angles
+ larger upper limit for defects cutoff
+ all domains mapped onto vacuum reporting to OUTPUT file
+ user manual improvements

Version 4.04: November 2012
---------------------------

+ extended EAM (EEAM) for metals, thanks to Ruslan Davidchak at Leicester University (UK)
+ optional reading of TABEAM embedding functions, interpolated over square root of density
+ improved memory management when handling metal interactions
+ KIHS Tersoff
+ shared units compression
+ metals tables referencing synchronised with DL\_POLY\_Classic
+ extra verification and error checks
+ user manual updates
+ fixes to:
  + forces on adiabatic shells in gravitational field, thanks to Sherwin Singer at Ohio State University (USA)
  + metals potentials
  + RDF calculation for massless particles
  + reporting topology (report\_topology)
  + reading angle potential parameters (and to user manual), thanks to Özgür Yazaydin at University of Surrey (UK)
  + shared RBs' communications
  + frozen RBs' setup
  + CGM minimisation of RBs
  + RBs' COM stress and virial estimators at start
  + Velocity Verlet NPT/NsT Berendsen iteration cycles for systems with constraints
  + Leapfrog Verlet NVT Berendsen/Nose-Hoover/Langevin iteration cycles for systems with both RBs and constraints
  + Velocity Verlet NsT Langevin RBs' COM stress
  + Velocity Verlet NVT Langevin RBs' quaternion momentum
  + Leapfrog Verlet NVT Langevin RBs' rotational motion
  + safety check when quenching adiabatic core-shell units, thanks to Ian Bush at NAG Ltd. (UK)

Version 4.03.4: June 2012
-------------------------

+ critical updates of CUDA port, thanks to Michael Lysaght at ICHEC (Dublin, Ireland)
+ minor fixes and improvements to eliminate syntax and concept inconsistencies, thanks to Yaser Afshar at JGU (Mainz, Germany) and Michael Seaton at STFC DL (UK)
+ modification to handling of metal interactions in line with DL\_POLY\_Classic
+ essential revisions to topology verification in 'nfold' option in line with DL\_FIELD standards, thanks to Chin Yong at STFC DL (UK)

Version 4.03.3: May 2012
------------------------

+ fix in set\_halo\_particles to correct electrostatics in parallel
+ non-critical change to system\_expand
+ improved dealing with EAM densities in metal\_ld\_collect\_eam with no smooth decent to or ascent from zero
+ increase in maximum number of variables in stack arrays in set\_bounds
+ backport of Java GUI to Java 1.6 (to replace syntax only available in Java 1.7)

Version 4.03.2: March 2012
--------------------------

+ changes in core algorithms due to exposure of ill-defined conditions, thanks to Eva Zarkadoula at Queen Mary University of London (UK)
+ fixes to Ewald frozen-frozen corrections, thanks to Henry Foxall at University of Sheffield (UK)
+ fixes to freedom distribution problem for fully frozen constraints, thanks to Sahar Mirshamsi at University of Florida (USA)
+ code refactoring
+ updates to user manual

Version 4.03.1: February 2012
-----------------------------

+ essential fixes to CUDA port
+ corrections to user manual

Version 4.03: January 2012
--------------------------

+ corrections and extensions to extended NsT ensembles
+ recasting of 'nfold' to better deal with intra-molecularly/topologically rich systems, thanks to Dr Gontrani Lorenzo at University of Rome 'La Sapienza' (Italy)
+ extra FIELD verification procedures
+ new 'mxstep' CONTROL option for variable timestep algorithm
+ speed-ups of constraint solvers, thanks to Valene Pellissier at NAG Ltd. (UK)
+ PMF units now allowed to be larger than major cutoff (but still smaller than half smallest MD cell width)
+ new Java GUI with manual and new features, thanks to Bill Smith:
  + cluster search in CONFIG file
  + comparison of two CONFIG files
  + new graph plotter with added functionality
  + in-built text editing capability
+ FENE bond breaking check
+ newly redefined Ewald intra-molecular and frozen-frozen exclusion routines
+ redefinition/reuse of linked list, leading to extension of RDF routine to full RDFs for intra-molecular species and performance boost for bio-chemical systems
+ LFV NPT and NsT Hoover, Langevin and MTK refinements
+ new Gentle Stochastic thermostat
+ better handling of TABLE reading and problem reporting
+ fixes to:
  + cosine angle potential bug
  + calcite potentials bug, thanks to Özgür Yazaydin at University of Surrey (UK)
  + bug reading large RBs, thanks to Adam Rigby at University of Manchester (UK)
  + free communicators in io\_module, thanks to Michael Seaton and Stephen Pickles at STFC DL (UK)
+ many user manual corrections and improvements

Version 4.02: July 2011
-----------------------

+ addition of residual CGM core-shell minimisation forces, thanks to Omololu Akin-Ojo at ICTP (Trieste, Italy)
+ new 'no topology' option in CONTROL to minimise reporting for bio-chemical simulations
+ new 4-body intra-atomic calcite potential as inversion angle
+ new option to generate displacements tracking file RSDDAT (similar to HISTORY and DEFECTS)
+ netCDF I/O improvements to avoid reprinting constant information and use single head node for single/serial operations
+ I/O buffer and batch size defaults modified
+ RSD available by replaying HISTORY without available velocity data
+ new Java GUI (beta release) available
+ fixes to:
  + typos in read\_field PMF parsing, thanks to Marco Molinari at University of Bath (UK)
  + typos in writing DEFECTS and reading REFERENCE, thanks to Eva Zarkadoula at Queen Mary University of London (UK)
  + Domain Decomposition exception, thanks to Emma Cai (project pupil) at STFC DL (UK)
  + typo in angle\_forces, thanks to Hamid Mosaddeghi at Isfahan University of Technology (Iran)
  + typo in read\_field when parsing long RBs, thanks to Lorenzo Gontrani at University of Rome (Italy)
  + typos and problematic definitions in shifted VDW interactions and long-range contributions, thanks to Marcelo Sepliarsky at IFIR CONICET (Argentina)
  + typos in read\_field for excluded/non-excluded interactions, thanks to Volodymyr Babin at North Carolina State University (USA)
  + typos in reporting correct energy units, thanks to Chin Yong at STFC DL (UK)
  + potential problem calculating restrained harmonic force in bonds\_forces, thanks to Yaser Afshar at Isfahan University of Technology (Iran)
  + inconsistency in system\_init reading NST Langevin data, thanks to Yaser Afshar at Isfahan University of Technology (Iran)
  + improve allocatable arrays' upper limit in two\_body\_forces and elsewhere, thanks to Yaser Afshar at Isfahan University of Technology (Iran)
  + bug in system\_expand affecting nfold option, thanks to Stephen Pickles at STFC DL (UK)

Version 4.01.1: November 2010
-----------------------------

+ fixes to:
  + issue in io\_module, thanks to Dr Andres De Virgilis at Max Planck Institute for Polymer Research (Mainz, Germany)
  + inconsistencies in set\_bounds and read\_control, thanks to Yaser Afshar at Isfahan University of Technology (Iran)
  + old style formatting sequences in CUDA f90 files, thanks to Yaser Afshar at Isfahan University of Technology (Iran)

Version 4.01: October 2010
--------------------------

+ rigid body (RB) dynamics from DL\_POLY\_2
  + not part of constraint bond, PMF constraint or contain a shell
  + available for all integrators
  + fully and semi-frozen RBs available
  + testing thanks to Laurence Ellison at STFC DL (UK)
+ velocity fields generated by Box-Mueller method in iso-kinetic manner for free particles (FPs) and RBs
+ parallel reading of CONFIG, thanks to Ian Bush at NAG Ltd. (UK) with funding from EPSRC's dCSE programme
+ reading and writing of CONFIG-like files available in netCDF binary (Amber-like format)
+ CGM timestep for relaxed shell model independent of simulation timestep (related to stiffest core-shell spring coefficient)
+ direct numerical evaluation option of untabulated VDW and metal interactions available
+ option to force-shift VDW interactions to zero at rvdw
+ LJ bonds for Amber force fields
+ Coulomb bonds and electrostatic exclusions available for all electrostatic options
+ zero Kelvin optimisation for LFV integrators similar to DL\_POLY\_2
+ iterative solver-related changes to velocities and forces (CBs, PMFs) in all integrators
+ velocity-dependent external force fields applied on half-timestep velocity
+ efficient domain decomposition of SPME electrostatics for processor numbers in each direction as multiples of 2, 3 and/or 5 (upgrade to 1-D FFT solver used in DaFT), thanks to Ian Bush at NAG Ltd. (UK) with funding from EPSRC's dCSE programme
+ CUDA+OpenMP port, thanks to Christos Kartsaklis and Ruairi Nestor at ICHEC (Dublin, Ireland)
+ Microsoft port with Microsoft self-installers (MSI), thanks to Igor Kozin at STFC DL (UK)
