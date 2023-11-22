#!/bin/sh

# defaults
rank=0
procs=1

# parse options
err=0
while test $# -ne 0
do case $1 in
   -h*)
      printf "\nUsage: $0 [-rank #] [-procs #]\n\n"
      shift 1;;
   -rank)
      if test -n "$2"
      then rank=$2
           shift
      else echo "ERROR: rank not specified"
           err=1
      fi
      shift;;
   -procs)
      if test -n "$2"
      then procs=$2
           shift
      else echo "ERROR: number of processes not specified"
           err=1
      fi
      shift;;
   esac
done
test $rank -lt 0 \
 && echo "ERROR: rank = $rank" \
 && err=1
test $rank -ge $procs \
 && echo "ERROR: rank = $rank (procs = $procs)" \
 && err=1
test $procs -le 0 \
 && echo "ERROR: procs = $procs" \
 && err=1
test $err -ne 0 && exit 1

xj="xjoda -rank $rank -procs $procs"

lastgeom=0
while test $lastgeom -eq 0
do rm -f MOINTS GAMLAM MOABCD SECDER DERINT
   (env $xj && xvmol && xvmol2ja && xvscf && xdirmp2 && xa2proc clrdirty) \
    || exit 1
   lastgeom=`xa2proc jareq LASTGEOM i 1 | tail -1`
done
exit 0

################################################################################

# Sample finite difference script on a Compaq Alpha:
alias pr='prun -l "%_g:"'
fdfile="$PWD/fd.out"
rm -f shared.* $fdfile
pr xgemini -s -i 'dirmp2fd.sh -rank @GRANK@ -procs @GPROCS@' || exit 1
for dir in shared.*
do cd $dir && xa2proc parfd updump >> $fdfile && cd - || exit 1
done
cd shared.0 && xa2proc parfd load $fdfile && xjoda -procs ${PROCS} || exit 1
pr xgemini -s -x

