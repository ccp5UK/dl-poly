# Building notes
* these notes are for building with [**cmake**](https://cmake.org)
* you can pass options to cmake using **-DOPTION=value**. For a complete list of options inspect [cmake/DLPOLYBuildOptions.cmake](cmake/DLPOLYBuildOptions.cmake)
* cmake -L <path to CMakeLists.txt> will show you a list of all available options.
* explicit compiler specification can be achieved by using environment variable **FC** (eg. using Intel ifort *FC=ifort*)
* compiler flags can be altered via **FFLAGS**, (eg *FFLAGS="-O3 -xHost"*)
* one also can use **cmake-gui** or **ccmake** to setup the build options
* to change the install path use **-DCMAKE_INSTALL_PREFIX=<path>** (*-DCMAKE_INSTALL_PREFIX=$HOME/101/DL_POLY*)
* automatic testing can be done after **DL_POLY_4** is built, using **make test**
* to see all the tests available use **ctest -N**
* to run one specific test use **ctest -R <TESTNAME>**
* for a list of all supported targets **make help**
* TODO check it works on Windows...

## Standard MPI
```sh
mkdir build-mpi-pure
pushd build-mpi-pure
FFLAGS="-O3" cmake ../
make -j10
make install
```
* will use whatever default MPI is found

*Intel Compilers - Intel MPI*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DMPI_Fortran_COMPILER=mpiifort 
```
*Intel Compilers - Some default mpi library, other than Intel MPI*
```sh
FC=ifort FFLAGS="-O3" cmake ../ 
```

## Hybrid MPI and OpenMP
```sh
mkdir build-mpi-openmp
pushd build-mpi-openmp
FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON
make -j10
make install
```
*Intel Compilers - Intel MPI*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON -DMPI_Fortran_COMPILER=mpiifort
```

## Serial
```sh
mkdir build-serial
pushd build-serial
FFLAGS="-O3" cmake ../ -DWITH_MPI=OFF
```
*Intel Compilers*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DWITH_MPI=OFF
```

## Serial with OpenMP threads
```sh
mkdir build-openmp
pushd build-openmp
FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON -DWITH_MPI=OFF
```
*Intel Compilers*
```sh
FC=ifort FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON -DWITH_MPI=OFF
```
## Xeon Phi
```sh
mkdir build-xeonphi
pushd build-xeonphi
FC=ifort FFLAGS="-O3 " cmake ../ -DWITH_PHI=ON -DWITH_MPI=ON
```

## Optimisation flags
* gfortran

```sh
FFLAGS="-O3 -mtune=native"
```

* Intel

```sh
FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15"
```

* If you plan to run the binary on a different type of a machine than you build it, check the manual of your compiler
for the flags matching the _running machine_

## Debugging, or when things go merdre
* gfortran

```sh
FFLAGS="-g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=35 -ffpe-trap=invalid,zero,overflow -fdump-core"
```
* Intel

```sh
FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv"
```

## Building with NETCDF support
```sh
mkdir build-mpi-netcdf
pushd build-mpi-netcdf
FFLAGS="-O3" cmake ../ -DWITH_NETCDF=ON
make -j10
make install
```

## Building with KIM support
```
mkdir build-mpi-kim
pushd build-mpi-kim
FFLAGS="-O3" cmake ../ -DWITH_KIM=ON
make -j10
make install
```

## Building with PLUMED support
```sh
mkdir build-mpi-plumed
pushd build-mpi-plumed
FFLAGS="-O3" cmake ../ -DWITH_PLUMED=ON
make -j10
make install
```

# Using the git for development

* The development model and structure

We have few concepts.

**Project git**: the main project where all the code and tests live
```sh
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git dl-poly-internal
```
* use ssh rather than https (self signed key issues)

**User fork**: a copy of the project git, which is your working copy.
create the fork on the web ui
Assuming the user is *alin*

```sh
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin
```

**Rationale**: While the _Project_ git and the _Fork_ git contain the same info, consider the fork
as your playground from where once you are happy you can make _Merge Requests_ to the Project.
You shall be able to play as much as you like with the code but keep the Project git clean and tidy.

## Branch, fix, merge  model:
Let us assume you have an issue with yout code which needs to be fixed.

* **Step 1**: Branch from your fork (We assume the fork is up to date with Project git)
- from the dashboard of your project you can branch from the + icon (issueXYZ)
clone the branch

```sh
git clone -b issueXYZ --single-branch ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin-issueXYZ
```

or
- create branch from cli

```sh
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin-issueXYZ
pushd dl-poly-alin-issueXYZ
git checkout -b issueXYZ
git push -u origin issueXYZ
```

* **Step 2**: fix
fix whatever is wrong.
use git status to see what is changed

```sh
git add <filename|folder> to add the new things
git commit -m "[tag] my cool message"
git push
```

Now you are ready for
* **Step 3a**: merge
simplest way is to go to webui and merge the branch
look for your Forked project, list the branches. Compare and if happy click Merge Request.
now all shall be created and one shall be able to accept the merge.

One can merge in two places.
1. Fork git
2. Project git

Each of them has its own merits. If you choose 2 be sure you sync your fork with the latest changes from Project.
If you choose 1 be sure you create another merge request when happy to the Project git

* **Step 3b**. Accept the merge
If all is ok merging is one click business. If not,
you will need to merge the things by hand
clone your project git and then follow the instructions gitlab indicates.
in short, if you tried to merge to the master branch of Project git:

```sh
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git dl-poly-internal
pushd dl-poly-internal
git fetch ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git issueXYZ
git checkout -b alin/dl-poly-issueXYZ FETCH_HEAD
git checkout master
git merge --no-ff alin/dl-poly-issueXYZ
fix the resulting conflicts if any
commit with a [git] tag
git add <whatever you fixed>
git commit -m "[git] fix commits"
git push origin master
```

## Advanced merdre
* Keep your fork in sync with the Project

```sh
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin
pushd dl-poly-alin
git remote add project ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git
```

these commands need to be done only once git remote -v shall show you the origin and project fetch and push links

```sh
[10:07:22 alin@baphomet: ...dl-poly-alin-master]: git remote -v
origin  ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git (fetch)
origin  ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git (push)
project ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git (fetch)
project ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git (push)
```

* when you need to sync

```sh
git pull
git fetch project
git checkout master
git merge project/master
git push
```

of course one can try to merge any other branch or available projects.

* rebasing commits
you worked on your issue and you have pushed few commits, eg. 5 , in the branch
squash them

```sh
git rebase -i HEAD~5
follow the instructions. Pick the first commit then s or f the rest.
git push origin issueXYZ --force
```

* cleaning stale branches. Deleting branches from the interface will get rid of the remotes and not of your
local copies. Command *git branch -a* will list remotes which are long gone. These are called stale branches. To get rid of them

```sh
git remote prune origin
```

* delete a local and remote branch. When you have merged a request forgot to delete the origin. You can delete the branch from the
web ui or from command line as:

```sh
git push origin :issueXYZ
```

to delete a local branch

```sh
git branch -d localBranch
```

if unmerged stuff exists but you still want to delete

```sh
git branch -D localBranch
```

# Code Coverage
if one builds DL_POLY_4 with **-DWITH_COVERAGE=ON** two targets will be available *make coverage* and *make runcoverage*.
First will run the code coverage on all tests from *make test*.

*make runcoverage* will run on the inputs which are put by user in **CodeAnalysis**. If one uses MPI **-DMPI_NPROCS**, default 4,
controls on how many processes the job will run.

# Workflow

Simple principles:
 * No commit without an issue
 * No merge request without a review
 * No new feature without a test case
 * Open issues for design discussions too. Things like the questions Aidan had about link cells/domain decomposition. The answers may be important for newer members. (these shall use a
   label called: _design_)
 * New features, eg task parallelism by Aidan, shaped particles by Vlad, shall have an issue too, the comments to the issue shall provide a succint
   progress on. (this shall use a label called: _leader_)

There are three labels, from which one shall be assigned to each issue
 * Testing: like in testing the water, shall be attached to things one may think to make one day.
 * Development: like this is what I am working now on. 
 * Production: like anything in this is critical and shall be fixed asap.

Adapting these principles shall give us the chance to have a quick overview of the state here.
https://ccforge.dl.ac.uk/dl-poly/dl-poly/boards (STFC networks only)
