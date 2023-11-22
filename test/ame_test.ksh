# typical usage of this script would be in `gmake $(./ame_test.ksh xvcc)`,
# which would run only those tests that call xvcc.

if test $# -eq 0
then echo FORCE
     echo "Usage: $0 vcc lambda ..." 1>&2
fi

list="(null"
for ame in $*
do list="${list}|${ame}"
done
list="${list})"

egrep -l "@ACES2: Executing .x*${list}" OUTPUT/out* | sed 's/OUTPUT.out.//'
