KIM API Installation Instructions:
----------------------------------

0) Presumptions

	Linux OS (for wget downloads)
	internet
	starting in the DL_POLY_4.ZZ main directory

1) Download and unpack the kim-api-v1.X.X.tgz (1.6.0 or later) tarball.

	% wget http://s3.openkim.org/kim-api/kim-api-v1.X.X.tgz
	% mv kim-api-v1.6.X.tgz DL_POLY_4.ZZ/source/
	% cd DL_POLY_4.ZZ/source
	% tar xzvf kim-api-v1.6.X

2) Download and unpack the testing Model and its Model Driver.

	% cd kim-api-v1.X.X/src/model_drivers
	% wget http://www.aem.umn.edu/~elliott/Pair_Morse__MD.tgz
	% tar xzvf Pair_Morse__MD.tgz
	% cd ../models
	% wget http://www.aem.umn.edu/~elliott/Pair_Morse_GirifalcoWeizer_LowCutoff_Cu__MO.tgz
	% tar xzfv Pair_Morse_GirifalcoWeizer_LowCutoff_Cu__MO.tgz
	% cd ../../

3) Configure, build and install the KIM API package according to the
   instructions in the INSTALL file.  Also "make add-" any/all OpenKIM
   Models you would like to have available.

	% cp Makefile.KIM_Config.example Makefile.KIM_Config
	% vi Makefile.KIM_Config  # Edit file as appropriate.
	                          # In particular, change the
	                          # KIM_DIR variable.
	% make add-Pair_Morse_GirifalcoWeizer_LowCutoff_Cu__MO
	% make
	% make install
	% make install-set-default-to-v1

4) Copy the DL_POLY_4.ZZ makefile and compile.

	% cd ..
	% cp ../build/Makefile_SRL2_KIM Makefile
        % vi Makefile # Edit the KIM_API_DIR variable as appropriate.
	% make gnu
	% cd ../execute
	% ls -haltr   # The last listed file should be DLPOLY.Z
	              # and have the current date stamp

  Similarly, one may build DL_POLY_4 in parallel by copying
  ../build/Makefile_MPI_KIM and compiling accordingly.

Test & Verify KIM:
------------------

1) Obtain DL_POLY_4 TEST21 following the README.txt advice from the DL_POLY_4.ZZ/data directory.

	% cd ../data
	% ftp ftp.dl.ac.uk
	ftp username: anonymous
	ftp password: email@google.com        # enter your email
	ftp> cd ccp5/DL_POLY/DL_POLY_4.0/DATA
	ftp> binary                           # transfer mode
	ftp> get TEST21.tar.gz
	ftp> quit
	% tar xzvf TEST21.tar.gz

2) Move DL_POLY input files to the DL_POLY_4.ZZ/execute directory.

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
