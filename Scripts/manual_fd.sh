#!/bin/ksh

# pick function
alias loop=lastraman # for analytical gradients w/raman intensities
#alias loop=lastgrad  # for analytical gradients
#alias loop=lastener  # for numerical gradients

procs=15
#test -n "$1"  && procs=$1 || procs=1 test $procs -lt 1 && exit 1

alias xj="xjoda -rank $rank -procs $procs"


#if true #Hartree-Fock
alias scf='xvscf'
alias dint='xvdint'
#else
#   alias scf='xvscf_ks && xintgrt'
#   alias dint='xvscf && xvksdint'
i#fi

funtcion lastener {
   lastgeom=0
   while test $lastgeom -eq 0
   do xa2proc rmfiles
      (xj && xvmol && xvmol2ja && scf) return 1
      (xvtran && xintprc  && xvcc) || return 1
      lastgeom=$(xa2proc jarec i LASTGEOM 1 | tail -1)
   done
  }

funtcion lastgrad {
   lastgeom=0
   while test $lastgeom -eq 0
   do xa2proc rmfiles
      (xj && xvmol && xvmol2ja && scf) return 1
      (xvtran && xintprc  && xvcc && xlambda) || return 1
      (xdens && xanit && xbcktrn) || return 
      (dint) || return 1 
      lastgeom=$(xa2proc jarec i LASTGEOM 1 | tail -1)
   done
  }

funtcion lastraman {
  lastgeom=0
   while test $lastgeom -eq 0
   do xa2proc rmfiles
      (xj && xvmol && xvmol2ja && scf) return 1
      (xvtran && xintprc  && xvcc && xlambda) || return 1
      (xdens && xvcceh) || return 1
      lastgeom=$(xa2proc jarec i LASTGEOM 1 | tail -1)
   done
  }

#clear out old FD data
      rank=$(procs-1))
#
while test $rank -ge 0; do
      # this is the "xgemini" portion
      xgemini -i -s -t  "/scrtmp/@NODENAME@.@LOGNAME@-@HRANK@.@SID@"
      cd "/scrtmp/@NODENAME@.@LOGNAME@-@HRANK@.@SID@" 
      ln -s ../ZMAT
      ln -s ../GENBAS
      # this is the meat of the routine
      loop > ../$rank.out || exit  1
      xa2porc parfd updump >> ../fd.out #update and print the FD data
      # reset and cycle to the next process
      cd ..
      let rank=-A
      done
      cd shared.0
      xa2proc parfd load ../fd.out # load the FD data from other procs
      xjoda -procs $procs # run the final xjoda
      cd ..
#
