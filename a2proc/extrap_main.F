
      SUBROUTINE EXTRAP_MAIN(ARGS,DIMARGS,Bverbose)
      IMPLICIT NONE

      INTEGER DIMARGS
      CHARACTER*80 ARGS(DIMARGS)
      Character*80 Feature
      Logical bVerbose

      INTEGER NATOMS, zFREE, zSCFGD, zTOTGD
      Integer Ilst, Linblnk
      External Linblnk

#include "machsp.com"
#include "icore.com"
#include "istart.com"
#include "iuhf.com"

      IF (DIMARGS.EQ.0) THEN
         PRINT *, '@EXTRAP_MAIN: unspecified argument'
         CALL ERREX
      END IF

      Ilst = Linblnk(args(1))
      if (bVerbose) print *, 'datatype  = ',args(1)(1:Ilst)

      Feature = Args(1)(1:Ilst)

      CALL GETREC(1, "JOBARC", "NREALATM", 1, NATOMS)

      If (Feature(1:Ilst) .Eq. "gradient" .or. 
     +    Feature(1:Ilst) .Eq. "energy") Then

         zSCFGD = I0
         zTOTGD = zSCFGD + 3*NATOMS*IINTFP
         zFREE  = zTOTGD + 3*NATOMS*IINTFP
         IF (zFREE-I0.GT.ICRSIZ) cALL INSMEM('EXTRAP_MAIN',
     +                                        zFREE-1,ICRSIZ)
         CALL EXTRAP_GRAD(ICORE(zSCFGD),ICORE(zTOTGD),NATOMS,ARGS(1))

      Elseif (Feature(1:Ilst) .Eq. "amplitudes") Then

         Call Extrap_amps(Icore(I0),Icrsiz,Iuhf)

      Endif 

      RETURN
      END

