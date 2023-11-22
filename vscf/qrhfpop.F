
      SUBROUTINE QRHFPOP(IRRP,ILOC,NBAS,EVAL,EVEC,IPOSABS,LDIM2,
     &                   ILOCATE,SCR,ISPN,IOS)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION EVAL(1),EVEC(1),SCR(1),ZJUNK
      DIMENSION ILOCATE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /POPUL/ NOCC(8,2)
#include "symm2.com"
C
      INDX2(I,J,N)=I+(J-1)*N
C
      IF(NOCC(IRRP,ISPN).EQ.0) THEN
       WRITE(LUOUT,9980)IRRP,ISPN
9980   FORMAT(T3,'@QRHFPOP-F, Irrep ',I1,' for spin ',I1,
     &           ' has fewer electrons '
     &          ,'than number to be removed.',/,
     &           T3,'I can''t do this calculation.')
       CALL ERREX
      ENDIF
      WRITE(LUOUT,9987)ILOCATE(IPOSABS),ISPN,EVAL(IPOSABS)
9987  FORMAT(T3,'@QRHFPOP-I, Electron removed from orbital ',I3,
     &          ' of spin ',I2,/,
     &          T14,' with eigenvalue ',F10.5,'.') 
C
C CALCULATE ILOC ASSUMING THAT SOME EIGENVALUES HAVE BEEN REORDERED
C
      IPOSREL=INEWVC(IPOSABS)
      ILOC=1+NOCC(IRRP,ISPN)-IPOSREL
C
C  Now check to see if we need to switch some eigenvectors around,
C  i.e., we're popping out an electron other than in the HOMO.  We
C  need to switch eigenvectors in both the alpha and beta blocks
C  so the transformation and what not are correct. This gets a bit
C  complicated if more than one change is made to a particular irrep.
C
      IF(ILOC.GT.1)THEN
       IF(ILOC.GT.NOCC(IRRP,ISPN))THEN
        WRITE(LUOUT,9982)IRRP,ISPN,NOCC(IRRP,ISPN),ILOC
9982    FORMAT(T3,'@QRHFPOP-F, Irrep ',I1,' of spin ',I1,
     &            ' has only ',I3,
     &            ' occupied orbitals.',/,
     &         T3,'You asked to depopulate a non-existent ',
     &            'orbital (',I3,').')
       CALL ERREX
       ENDIF
       IOFF=(NOCC(IRRP,ISPN)+1)-ILOC
       INOC=NOCC(IRRP,ISPN)
       
       DO 10 ISPIN=1,2
        IF(IOS.NE.0.AND.ISPIN.NE.ISPN) GOTO 10
#ifdef _DEBUG_LVLM
        Write(6,*) "IOFF, INOC", IOFF, INOC
#endif 
        IPOS1=(ISPIN-1)*LDIM2+ISQROF(IRRP)-1+INDX2(1,INOC,NBFIRR(IRRP))
        IPOS2=(ISPIN-1)*LDIM2+ISQROF(IRRP)-1+INDX2(1,IOFF,NBFIRR(IRRP))
c YAU : old
c       CALL ICOPY(NBFIRR(IRRP)*IINTFP,EVEC(IPOS1),1,SCR,1)
c       CALL ICOPY(NBFIRR(IRRP)*IINTFP,EVEC(IPOS2),1,EVEC(IPOS1),1)
c       CALL ICOPY(NBFIRR(IRRP)*IINTFP,SCR,1,EVEC(IPOS2),1)
c YAU : new
        CALL DSWAP(NBFIRR(IRRP),EVEC(IPOS1),1,EVEC(IPOS2),1)
c YAU : end
C
C   Now switch around the eigenvalues.
C
        IOFF1=IRPOFF(IRRP)+NOCC(IRRP,ISPN)+1-ILOC
        IOFF2=IRPOFF(IRRP)+NOCC(IRRP,ISPN)
        ZJUNK=EVAL((ISPIN-1)*NBAS+IOFF1)
        EVAL((ISPIN-1)*NBAS+IOFF1)=EVAL((ISPIN-1)*NBAS+IOFF2)
        EVAL((ISPIN-1)*NBAS+IOFF2)=ZJUNK
C
   10  CONTINUE
       ITMP=ILOCATE(IOFF1)
       ILOCATE(IOFF1)=ILOCATE(IOFF2)
       ILOCATE(IOFF2)=ITMP
      ENDIF
      NOCC(IRRP,ISPN)=NOCC(IRRP,ISPN)-1
C
      RETURN
      END
