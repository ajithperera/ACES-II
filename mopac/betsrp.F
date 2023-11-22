C 
C    BETSRP unit (Ivan Rossi - April 94)                               
C
C include in mopcvar.par the two following paramrters:
C
C     MXATSP : MAXIMUM NUMBER OF atomic species in the system
C     MXSRPB : MAXIMUM NUMBER OF SRP special beta allowed
*
*     COMMON BLOCKS FOR SRP (Ivan Rossi - April 94)                    IR0494
*
*     COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
*    *       /SRPL/ ISSRP
*    *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
*     LOGICAL ISSRP
c
      SUBROUTINE INIBET
*
*  Initialize arrays and pointers for special SRP BETAs (Ivan Rossi- April '94)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mopcvar.par"                                                
C
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
*
*     COMMON BLOCKS FOR STANDARD BETAs
*
      COMMON /BETAS / BETAS(107),BETAP(107),BETAD(107)
c
c     Initialize pointers to RSP betas
c
      NATSP=0
      do 5 i=1,107
   5     IBTPTR(i)=0
      DO 10 i=1,NATOMS
        if (LABELS(i) .eq. 99 .or. LABELS(i) .eq.107) goto 10
        if (IBTPTR(LABELS(i)) .eq. 0) then
           NATSP=NATSP+1
           if (NATSP .gt. MXATSP) then
             Write(6,'("***  The parameter MXATSP is too SMALL ! ***")')    
             Write(6,'("*Increase it in mopcvar.par and recompile *")')      
             STOP
           ENDIF
           NATPTR(NATSP)=LABELS(I)
           IBTPTR(LABELS(i))=NATSP
        endif
 10   continue
        do 50 j=1,NATSP
          do 50 i=j,NATSP
          BSS= 0.5d0*(BETAS(NATPTR(i))+BETAS(NATPTR(j)))
          BPP= 0.5d0*(BETAP(NATPTR(i))+BETAP(NATPTR(j)))
          call SETBET( NATPTR(i), NATPTR(j), 'BETSS', BSS )
          call SETBET( NATPTR(i), NATPTR(j), 'BETPP', BPP )
 50     continue
        do 60 j=1,NATSP
          do 60 i=1,NATSP
          BSP= 0.5d0*(BETAS(NATPTR(i))+BETAP(NATPTR(j)))
          call SETBET( NATPTR(i), NATPTR(j), 'BETSP', BSP )
 60     continue
      RETURN
      END
*
      SUBROUTINE SETBET( NATM1, NATM2, BETTYP, BETVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*5 BETTYP
      DOUBLE PRECISION BETVAL
* 
*  Set the value of the special SRP BETAs ( Ivan Rossi - April '94 )
*      WARNING: BETSP(x,y) is DIFFERENT from BETSP(y,x)
*
*   INPUT:
*     NATM1, NATM2 : Atomic numbers of the atom couple
*     BETTYP : Name of the Beta parameter type (UPPERCASE!)
*     BETVAL : value of the BETA to set
*
#include "mopcvar.par"
C
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     *       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
c
      DIMENSION NATM(2)
C
        NATM(1)=NATM1
        NATM(2)=NATM2
C
C   INPUT Error checking
C
      do 10 j=1,2
        if (IBTPTR(NATM(j)) .eq. 0) then 
          WRITE(6,'("SETBET: Atom type",I4," NOT present!")') NATM(j)
cmgc          STOP
          return
        endif
 10   continue
c
c  Calculate address of beta (packed lower simmetric matrix) and set value
c
     
      IROW=IBTPTR(NATM(1))
      JCOL=IBTPTR(NATM(2))
      if( JCOL .gt. IROW ) then
        IROW=IBTPTR(NATM(2))
        JCOL=IBTPTR(NATM(1))
      ENDIF
      j=NATSP*(JCOL-1)+IROW -(JCOL*(JCOL-1))/2
      if( BETTYP .eq. 'BETPP') then
        BETPP(j)=BETVAL
      else if( BETTYP .eq. 'BETSP') then
        BETSP(IBTPTR(NATM1),IBTPTR(NATM2))=BETVAL
      else if( BETTYP .eq. 'BETSS') then
        BETSS(j)=BETVAL
      else
        write(6,'("SETBET: Wrong beta type",A5)') BETTYP
        STOP
      endif
      return
      END
*
      DOUBLE PRECISION FUNCTION GETBET( NATM1, NATM2, BETTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*5 BETTYP
* 
*  Get the value of the special SRP BETAs ( Ivan Rossi - April '94 )
*      WARNING: BETSP(x,y) is DIFFERENT from BETSP(y,x)
*
*   INPUT:
*     NATM1, NATM2 : Atomic numbers of the atom couple
*     BETTYP : Name of the Beta parameter type (UPPERCASE!)
*   OUTPUT:
*     Value of the SRP BETA requested
*
#include "mopcvar.par"
C
      COMMON /SRPI/ IBTPTR(107), NATPTR(MXATSP), NATSP
     1       /SRPR/ BETSS(MAXBET), BETSP(MXATSP,MXATSP), BETPP(MAXBET)
c
      DIMENSION NATM(2)
C
        NATM(1)=NATM1
        NATM(2)=NATM2
C
C   INPUT Error checking
C
      do 10 j=1,2
        if (IBTPTR(NATM(j)) .eq. 0) then 
          WRITE(6,'("GETBET: Atom type",I4," NOT present!")') NATM(j)
          STOP
        endif
 10   continue
c
c  Calculate address of beta and set value
c
      IROW=IBTPTR(NATM(1))
      JCOL=IBTPTR(NATM(2))
      if( JCOL .gt. IROW ) then
        IROW=IBTPTR(NATM(2))
        jCOL=IBTPTR(NATM(1))
      ENDIF
      j=NATSP*(JCOL-1)+IROW -(JCOL*(JCOL-1))/2
      if( BETTYP .eq. 'BETPP') then
        GETBET=BETPP(j)
      else if( BETTYP .eq. 'BETSP') then
        GETBET=BETSP(IBTPTR(NATM1), IBTPTR(NATM2))
      else if( BETTYP .eq. 'BETSS') then
        GETBET=BETSS(j)
      else
        write(6,'("GETBET: Wrong beta type",A5)') BETTYP
          STOP
      endif
      return
      END
