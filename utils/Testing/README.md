#Generate test data ([../../cmake/all.xml.cmake](../../cmake/all.xml.cmake))

*download the tests by using **-DBUILD_TESTING=On** 

*build dlpoly with MPI support 

*change the [rerun.sh](rerun.sh) script to match your **DLPOLY.Z** binary and data

*if all is fine edit in [generateTests.sh](generateTests.sh) data and DLROOT variables and then run it

*[all.xml.cmake](all.xml.cmake) containing all the reference data is in $DLROOT/cmake/ now

*now things are ready and you can test again with _make test_

