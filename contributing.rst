Contribution workflow and the review process
============================================

This document outlines the best practice, using git and CI, which
**must** be followed for all contributions to DL_POLY. Also contained
are instructions and tips for managing your fork of the project which
will help keep merges clean and avoid many headaches.

Golden rules
------------

In brief the rules for contribution are as follows,

-  Follow the branch, fix, merge model, from your own fork
-  An issue must be created for every piece of work (bug, feature,
   *etc.*)
-  Merge requests will not be accepted without a review
-  New features must have a test
-  All tests must pass, no regressions may be merged

Issues
------

Using issues
~~~~~~~~~~~~

-  Open an issue for each piece of work done.
-  Open issues for design discussions. For example the questions Aidan
   had about link cells/domain decomposition. The answers may be
   important for newer members.
-  New features, *e.g.* task parallelism by Aidan, shaped particles by
   Vlad, shall have an issue too, the comments shall be used to provide
   succint reports on progress.

Labels
~~~~~~

Labels may be assigned to issues to help classify them. Examples
include,

-  Testing: as in testing the water. This label shall be attached to
   things one may think to make one day.
-  Development: this is what I am working now on.
-  Production: anything in this is critical and shall be fixed asap.
-  Design: queries or suggestions about the structure of the program.
-  Leader: for issues relating to new features.

Review
------

All merge requests will be reviewed to ensure the integrity of the code.

The reviewer/s have the following responsibilities, \* Ensuring all
contribution rules have been followed \* Ensuring the `coding
style <./coding_style.md>`__ is adhered to \* Only accepting a merge if
all tests have passed \* Using the comments system to request changes
for the submittor to make

Using the git for development
-----------------------------

The Gitlab instance hosts a *devel* repository, which we will refer to
as *devel*. Contributors will work on their own personal copies of the
repository by creating *forks*. This allows us to keep *devel* clean
(one commit per merge request, if possible, all commits passing tests)
while users may work on their own *fork*, creating commits and changing
the code as they see fit.

The *devel* repository may be cloned as follows,

.. code:: sh

   git clone git@gitlab.com:ccp5/dl-poly.git dl-poly-devel

A *fork* is created using the web UI. It may then be cloned for a user
called ‘username’ as follows,

.. code:: sh

   git clone git@gitlab.com:username/dl-poly.git dl-poly-fork

Branch, fix, merge model:
~~~~~~~~~~~~~~~~~~~~~~~~~

All work should follow the workflow of branch, fix, merge. Let us assume
you have an issue with your code which needs to be fixed.

Step 1: Branch from your fork
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new branch for the issue on the dashboard of your fork, we will
assume the branch is called ‘issueXYZ’. Clone the branch,

.. code:: sh

   $ git clone -b issueXYZ --single-branch git@gitlab.com:username/dl-poly.git dl-poly-issueXYZ

Alternatively you can create the branch in the cli using

.. code:: sh

   # clone the repository, if you already have a local repository this is not nessecary
   $ git clone git@gitlab.com:username/dl-poly.git dl-poly-issueXYZ
   $ pushd dl-poly-issueXYZ
   # create and checkout a new branch (this is equivilent to git branch followed by git checkout)
   $ git checkout -b issueXYZ
   # create a remote tracking branch for you to push your changes to
   $ git push -u origin issueXYZ

Step 2: Fix the issue and commit your changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fix whatever is wrong. Use git status to see which files you have
changed and prepare a commit.

.. code:: sh

   # stage changes
   $ git add <filename|folder> to add the new things
   # commit the changes with a clear and brief message
   $ git commit -m "<commit message>"
   # push the commit to origin
   $ git push

Step 3a: Merge your branch into devel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On the web interface navigate to *devel* and create a merge request for
your branch on your *fork*. Add any relevant labels or milestones and
assign a reviewer. Compare the code and if you are happy click Submit
Merge Request.

After the merge request has been submitted tests will be run and your
reviewer will be notified.

Step 3b: Finalising the merge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If all is OK with the commit your reviewer may set the request to be
merged once all tests pass. Otherwise the reviewer may open discussions
using the Gitlab comment system to point out issues that may need to be
addressed before the commit can be merged.

If changes need to be made you may make more commits onto your branch.
When you push your branch to your *fork* the merge request will be
automatically updated to use the latest commit. Reply to the discussions
to indicate when and how they have been addressed.

If your branch has become out of sync with *devel* then conflicts may
arise. Sometimes these cannot be automatically resolved and you will
need to resolve them by hand. Gitlab provides instructions for this, or
you can follow this routine,

.. code:: sh

   # add devel as a remote if you have not already
   $ git remote add devel git@gitlab.com:ccp5/dl-poly.git
   # get the changes to devel since you started working on your issue
   $ git fetch devel
   # merge these changes into your branch (assuming you want to merge into the master branch on devel)
   $ git merge devel/devel
   # resolve any conflicts
   # push to your fork
   $ git push

Alternatively you may use rebase which will replay the changes you made
in your branch on top of *devel/devel* however be sure you understand
the differences between merge and rebase

.. code:: sh

   # add devel as a remote if you have not already
   $ git remote add devel git@gitlab.com:ccp5/dl-poly.git
   # get the changes to devel since you started working on your issue
   $ git fetch devel
   # merge these changes into your branch (assuming you want to merge into the master branch on devel)
   $ git rebase devel/devel
   # resolve any conflicts
   # push to your fork
   $ git push

Advanced git
~~~~~~~~~~~~

Keeping your fork in sync with project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By adding two remotes, one for *devel* and one for your *fork* it is
possible to keep your *fork* in sync with *devel*. This will greatly
simplify merge requests.

.. code:: sh

   # clone your fork
   $ git clone git@gitlab.com:username/dl-poly.git dl-poly-fork
   pushd dl-poly-fork
   # add a remote for devel
   $ git remote add devel git@gitlab.com:ccp5/dl-poly.git

These commands need to be done only once. ``git remote -v`` shall show
you the origin and project fetch and push links

.. code:: sh

   $ git remote -v
   origin  git@gitlab.com:username/dl-poly.git (fetch)
   origin  git@gitlab.com:username/dl-poly.git (push)
   devel git@gitlab.com:ccp5/dl-poly.git (fetch)
   devel git@gitlab.com:ccp5/dl-poly.git (push)

When you need to sync your *fork* with *devel*, do the following,

.. code:: sh

   # get the latest commits from devel
   $ git fetch devel
   # ensure you are in the master branch of your fork
   $ git checkout master
   # merge your master branch into the master branch of devel
   $ git merge devel/devel
   # push these changes back to the remote of your fork (origin)
   $ git push

of course one can use a similar process to merge any other branch or
available projects.

Rebasing commits
^^^^^^^^^^^^^^^^

When working on an issue you may use multiple commits. When you are
ready to create a merge request, you should squash your changes into one
commit in order to keep *devel* clean. This is most easily achieved with
an interactive rebase.

Assuming you have made five commits,

.. code:: sh

   # rebase your branch five commits before HEAD i.e. where your branch originally diverged
   $ git rebase -i HEAD~5
   # follow the instructions. 'pick' the first commit then 'sqaush' or 'fixup' the rest.
   # You should now be left with a single commit containing all your changes
   # Push your commmit to the remote, use --force if you have already pushed this branch to
   # 'rewrite history'
   $ git push origin branchname --force

using force is a powerful and dangerous option. use it only if you know
150% nobody else touched that branch.

Cleaning stale branches
^^^^^^^^^^^^^^^^^^^^^^^

Deleting branches from the web interface will get rid of the remotes and
not of your local copies. The local branches left behind are called
stale branches. To get rid of them

.. code:: sh

   $ git remote prune origin

To delete a local branch

.. code:: sh

   $ git branch -d localBranch

if unmerged commits exists but you still want to delete use

.. code:: sh

   $ git branch -D localBranch

To delete a remote branch on the remote *origin* use

.. code:: sh

   $ git push -d origin remoteBranch

Code Coverage
-------------

If one builds DL_POLY_4 with **-DWITH_COVERAGE=ON** two targets will be
available *make coverage* and *make runcoverage*. First will run the
code coverage on all tests from *make test*.

*make runcoverage* will run on the inputs which are put by user in
**CodeAnalysis**. If one uses MPI **-DMPI_NPROCS**, default 4, controls
on how many processes the job will run.
