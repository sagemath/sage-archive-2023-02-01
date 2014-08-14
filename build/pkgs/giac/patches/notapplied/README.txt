!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
The patches in this directory are not applied. 
   We keep them for memory in case they are needed again.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*   the configure.in patch was solving a linking problem.
    NB:  after application we need to rebuilt the scripts with:
     aclocal
     autoconf
     automake --add-missing
     libtoolize
