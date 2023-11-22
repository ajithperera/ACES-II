      SUBROUTINE SETBAS
C
C This subroutine fills the common block /BASSYM/. This take into
C account for the dorpped orbitals in drop core calcualtions.
C 
C  NBAS   ... Number of basis functions per irrep.
C  NBASIS ... Total number of basis functions
C  NBASSQ ... Size of a symmetry packed square matrix
C  NBASTT ... Size of a symmetry packed triangular matrix
C
C
      IMPLICIT INTEGER(A-Z)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      DIMENSION NDRPOP(8), NDRVRT(8) 
C
C Take care of the dropped orbitals in drop core gradient calculations
C
      call getrec(20, 'JOBARC', 'NUMDROPA', 1, NDROP)
C
      if (ndrop.ne.0) then
C
         call getrec(20, 'JOBARC', 'NDROPPOP', NIRREP, NDRPOP) 
         call getrec(20, 'JOBARC', 'NDROPVRT', NIRREP, NDRVRT) 
C
         DO 5 i=1, nirrep
            vrt(i,1) = vrt(i,1) + ndrvrt(i) 
            vrt(i,2) = vrt(i,2) + ndrvrt(i) 
            pop(i,1) = pop(i,1) + ndrpop(i) 
            pop(i,2) = pop(i,2) + ndrpop(i) 
 5       Continue 
C
      endif
C
      DO 100 IRREP = 1, NIRREP
C     
         NBAS(IRREP) = POP(IRREP,1) + VRT(IRREP,1)
C
 100  CONTINUE
C
      NBASIS=0
      NBASSQ=0
      NBASTT=0
C
      DO 200 IRREP = 1, NIRREP
C
         NBASIS = NBASIS + NBAS(IRREP)
         NBASSQ = NBASSQ + NBAS(IRREP)*NBAS(IRREP)
         NBASTT = NBASTT+(NBAS(IRREP)*(NBAS(IRREP)+1))/2
C
 200  CONTINUE
C
      RETURN
      END

