# Building notes 
* these notes are for building with [**cmake**](https://cmake.org)
* you can pass options to cmake using **-DOPTION=value**. For a complete list of options inspect [cmake/cmake/DLPOLYBuildOptions.cmake](cmake/cmake/DLPOLYBuildOptions.cmake) 
* explicit compiler specification can be achieved by using environment variable **FC** (eg. using Intel ifort *FC=ifort*)
* compiler flags can be altered via **FFLAGS**, (eg *FFLAGS="-O3 -xHost"*)
* one also can use **cmake-gui** to setup the build options
* to change the install path use **-DCMAKE_INSTALL_PREFIX=<path>** (*-DCMAKE_INSTALL_PREFIX=$HOME/101/DL_POLY*)
* automatic testing can be done after **DL_POLY_4** is built using **make test**
* for a list of all supported targets **make help**
* TODO check it works on Windows...

## Standard MPI 
```
mkdir build-mpi-pure
pushd build-mpi-pure
FFLAGS="-O3" cmake ../ -DWITH_MPI=ON
make -j10 
make install
```
* will use whatever default MPI is found

*Intel Compilers*
```
FC=mpiifort FFLAGS="-O3" cmake ../ -DWITH_MPI=ON
```

## Hybrid MPI and OpenMP
```
mkdir build-mpi-openmp
pushd build-mpi-openmp
FFLAGS="-O3" cmake ../ -DWITH_MPI=ON -DWITH_OPENMP=ON
make -j10 
make install
```
*Intel Compilers*
```
FC=mpiifort FFLAGS="-O3" cmake ../ -DWITH_MPI=ON -DWITH_OPENMP=ON
```

## Serial 
```
mkdir build-serial
pushd build-serial
FFLAGS="-O3" cmake ../ 
```
*Intel Compilers*
```
FC=ifort FFLAGS="-O3" cmake ../ 
```

## Serial with OpenMP threads
```
mkdir build-openmp
pushd build-openmp
FFLAGS="-O3" cmake ../ -DWITH_OPENMP=ON
```
*Intel Compilers*
```
FC=ifort FFLAGS="-O3" cmake ../ 
```
## Xeon Phi
```
mkdir build-xeonphi
pushd build-xeonphi
FC=ifort FFLAGS="-O3 " cmake ../ -DWITH_PHI=ON -DWITH_MPI=ON
```

## Optimisation flags
* gfortran 
```
FFLAGS="-O3 -mtune=native"
```
* Intel
```
FFLAGS="-fpp -O3 -xHost -fimf-domain-exclusion=15"
```

* If you plan to run the binary on a different type of a machine than you build it, check the manual of your compiler 
for the flags matching the _running machine_

## Debugging, or when things go merdre
* gfortran
```
FFLAGS="-g -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483648 -finit-real=nan -finit-logical=true -finit-character=35 -ffpe-trap=invalid,zero,overflow -fdump-core"
```
* Intel 
```
FFLAGS="-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv"
```

## known issues: 
1. FindMPI may trip with Intel MPI
```
I_MPI_F90=ifort FC=ifort 
```
shall cure all issues.

# Using the git for development

* The development model and structure

We have few concepts.

**Project git**: the main project where all the code and tests live
```
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git dl-poly-internal
```
* use ssh rather than https (self signed key issues)

**User fork**: a copy of the project git, which is your working copy.
create the fork on the web ui
Assuming the user is *alin*
```
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
```
git clone -b issueXYZ --single-branch ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin-issueXYZ
```

or
- create branch from cli
```
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin-issueXYZ
pushd dl-poly-alin-issueXYZ
git checkout -b issueXYZ
git push -u origin issueXYZ
```

* **Step 2**: fix
fix whatever is wrong.
use git status to see what is changed
```
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
```
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

```
git clone ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git dl-poly-alin
pushd dl-poly-alin 
git remote add project ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git
```

these commands need to be done only once git remote -v shall show you the origin and project fetch and push links

```
[10:07:22 alin@baphomet: ...dl-poly-alin-master]: git remote -v
origin  ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git (fetch)
origin  ssh://gitlab@ccforge.dl.ac.uk:1980/alin/dl-poly.git (push)
project ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git (fetch)
project ssh://gitlab@ccforge.dl.ac.uk:1980/dl-poly/dl-poly.git (push)
```

* when you need to sync
```
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
```
git rebase -i HEAD~5
follow the instructions. Pick the first commit then s or f the rest.
git push origin issueXYZ --force
```

* cleaning stale branches. Deleting branches from the interface will get rid of the remotes and not of your 
local copies. Command *git branch -a* will list remotes which are long gone. These are called stale branches. To get rid of them 
```
git remote prune origin
```

* delete a local and remote branch. When you have merged a request forgot to delete the origin. You can delete the branch from the 
web ui or from command line as:

```
git push origin :issueXYZ 
```
to delete a local branch 
```
git branch -d localBranch
```
if unmerged stuff exists but you still want to delete
```
git branch -D localBranch
```

# Code Coverage
if one builds DL_POLY_4 with **-DWITH_COVERAGE=ON** two targets will be available *make coverage* and *make runcoverage*.
First will run the code coverage on all tests from *make test*.

*make runcoverage* will run on the inputs which are put by user in **CodeAnalysis**. If one uses MPI **-DMPI_NPROCS**, default 4,
controls on how many processes the job will run.
