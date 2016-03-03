Test & Verify KIM:
------------------

1) Obtain DL_POLY_4 TEST21 following the README.txt advice from the DL_POLY_4.Z/data directory.

	% cd ../data
	% ftp ftp.dl.ac.uk
	ftp username: anonymous
	ftp password: email@google.com        # enter your email
	ftp> cd ccp5/DL_POLY/DL_POLY_4.0/DATA
	ftp> binary                           # transfer mode
	ftp> get TEST21.tar.gz
	ftp> quit
	% tar xzvf TEST21.tar.gz

2) Move DL_POLY input files to the DL_POLY_4.Z/execute directory.

	% cd ../execute
	% cp ../data/TEST21/CONFIG .
	% cp ../data/TEST21/CONTROL .
	% cp ../data/TEST21/FIELD .

3) Edit FIELD and change

	metal 1 table
	Cu      Cu      eam
	close

   to

	vdw 1
	Cu      Cu      mors 3.42900e-01 2.86600e+00 1.35880e+00
	close

4) Edit CONTROL and:

   remove

	restart scale

   change

	steps               6000     steps
	equilibration       5000     steps

   to

	steps                100     steps
	equilibration         10     steps

   and

	cutoff                 5.0   Angstroms

   to

	cutoff                 6.76342e+00   Angstroms

5) Run and save the "native" simulation.

	% ./DLPOLY.Z
	% mv OUTPUT OUTPUT.native

6) Edit FILED and change

	vdw 1
	Cu      Cu      mors 3.42900e-01 2.86600e+00 1.35880e+00
	close

   to

	kim  Pair_Morse_GirifalcoWeizer_LowCutoff_Cu__MO

7) Run and save the "kim" simulation.

	% ./DLPOLY.Z
	% mv OUTPUT OUTPUT.kim

8) Compare simulation outputs form the native and KIM interaction models.

	% diff OUTPUT.native OUTPUT.kim

Good luck!

March 2015
Ryan Elliott, Ilian Todorov and Henry Boateng
