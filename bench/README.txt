Benchmark Test Cases for DL_POLY_4
----------------------------------

Interested users can obtain the benchmark test cases
from the CCP5 FTP server as follows:

FTP site : ftp.dl.ac.uk
Username : anonymous
Password : your email address
Directory: ccp5/DL_POLY/DL_POLY_4.0/BENCH

The user may as well download BENCH.tar.gz containing the
whole distribution.  The directories within contain small
samples of the systems for benchmarking. To enlarge the
systems one can use two approaches:

(i ) Enlarge them by using the F90 source that will only
     apply to the CONFIG or *.inp crystallographic file.
     The FIELD must then be amended respectively by hand
     according to the extent of the enlargement.

(ii) Enlarge the system by using the DL_POLY_4's 'nfold'
     option in CONTROL by carrying out a dry run (0
     timesteps).  Since the starting samples are quite
     small the cutoff value in the corresponding CONTROLs
     must be augmented by hand to 2 A for a successful
     enlargement and then reverted back to its original
     value.  This task can be done in serial as well as
     in parallel.  Consult the user manual for more
     information on the aforementioned options.
