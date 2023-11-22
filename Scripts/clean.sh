
# Usage: clean.sh [options]
# options: (bracketed text is default)
#    [+]/-out = [leave]/remove *out
#    [+]/-src = [leave]/remove source code (e.g., *.F *.f *.c ...)
#    +/[-]mos = leave/[remove] *NEWMOS *OLDMOS *GUESS
#    +/[-]aos = leave/[remove] *AOBASMOS *OLDAOMOS
#    +/[-]arc = leave/[remove] JOBARC JAINDX

# This script ruthlessly deletes files not needed or recognized by ACES.
# Various flags may be passed in which leave or remove certain files.
# As a safety feature, if ZMAT and GENBAS do not exist, then clean.sh assumes
# this is not an execute directory and exits.
#    PLEASE examine the actual script to make sure nothing gets deleted by
# accident. Some folks like working on regular files and executing ACES in
# the same directory. This behavior is discouraged and clean.sh enforces that.
# - Anthony Yau

if [ -f ZMAT ] && [ -f GENBAS ] ; then
   rm_out=no
   rm_src=no
   rm_mos=yes
   rm_aos=yes
   rm_arc=yes
   for arg in "$@" ; do
      if [ "$arg" = "-out" ] ; then rm_out=yes ; fi
      if [ "$arg" = "-src" ] ; then rm_src=yes ; fi
      if [ "$arg" = "+mos" ] ; then rm_mos=no  ; fi
      if [ "$arg" = "+aos" ] ; then rm_aos=no  ; fi
      if [ "$arg" = "+arc" ] ; then rm_arc=no  ; fi
   done
   rm .molcas_info
   for file in * ; do
      case $file in
      ZMAT|GENBAS|*zmt)
      ;;
      *out)
         if [ "$rm_out" = "yes" ] ; then rm $file ; fi
      ;;
      JOBARC|JAINDX)
         if [ "$rm_arc" = "yes" ] ; then rm $file ; fi
      ;;
      *NEWMOS|*OLDMOS|*GUESS)
         if [ "$rm_mos" = "yes" ] ; then rm $file ; fi
      ;;
      *AOBASMOS|*OLDAOMOS)
         if [ "$rm_aos" = "yes" ] ; then rm $file ; fi
      ;;
      *job|*exe)
      ;;
      *bz2|*gz|*tar|*old|*.csv)
      ;;
      x*|lib*.a|*.h|*.cpp|*.C|*.c|*.F95|*.f95|*.F90|*.f90|*.F|*.f|*.o)
         if [ "$rm_src" = "yes" ] ; then rm $file ; fi
      ;;
      *)
         if [ -f $file ] || [ -L $file ] ; then rm $file ; fi
      ;;
      esac
   done
   echo ""
   /bin/ls -F
else
   printf "\nclean only acts on directories with ZMAT and GENBAS files.\n\n"
fi

