echo "Error summary:"
test -f "$1" && log=$1 || log=gmake.out
for out in out* out.null
do test -f $out || continue
   sed -n "/^ *test value =/ {
               /[DdEe]-0*8$/d;
               s/$/ (${out})/;
               p;
           }" $out
done
grep 'FAILED in [xs]' $log
echo "Finished error summary."
