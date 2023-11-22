      SUBROUTINE QRHFERR(A)
C
C  THIS ROUTINE WRITES AN ERROR MESSAGE IF THE GRADIENT CALCULATIONS
C  CAN NOT BE DONE FOR THE GIVEN QRHF PARAMETERS
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*(*) A
      DIMENSION ISET(255)
C
      WRITE(6,*) '@',A,'-E:  SORRY, ORBITAL RELAXATION IS NOT YET',
     1      ' PROGRAMMED FOR THIS QRHF CASE'
C
      IONE=1
      CALL GETREC(20,'JOBARC','QRHFTOT ',IONE,NMOD)
       WRITE(6,*) ' QRHFTOT: ',NMOD
      if (nmod.gt.255) then
         print *, '@QRHFERR: Assertion failed.'
         print *, '          iset dimension = 255'
         print *, '          nmod = ',nmod
         call errex
      end if
      CALL GETREC(20,'JOBARC','QRHFIRR ',NMOD,ISET)
       WRITE(6,*) ' QRHFIRR: ',(ISET(I),I=1,NMOD)
      CALL GETREC(-1,'JOBARC','QRHFSPN ',NMOD,ISET)
       WRITE(6,*) ' QRHFSPN: ',(ISET(I),I=1,NMOD)
      CALL GETREC(-1,'JOBARC','QRHFLOC ',NMOD,ISET)
       WRITE(6,*) ' QRHFLOC: ',(ISET(I),I=1,NMOD)
C
      RETURN
      END
