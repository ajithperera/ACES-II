
      SUBROUTINE SETBAS
C
C  THIS SUBROUTINE FILLS THE COMMON BLOCK /BASSYM/
C
C  ITS PARAMETER ARE
C 
C      NBAS ..... NUMBER OF BASIS FUNCTION PER IRREP
C      NBASIS ... TOTAL NUMBER OF BASIS FUNCTION
C      NBASSQ ... SIZE OF A SYMMETRY PACKED SQUARE MATRIX
C      NBASTT ... SIZE OF A SYMMETRY PACKED TRIANGULAR MATRIX
C
C  AND SET UP FOR CARTESIAN D AND F FUNCTIONS
C
C  THE INFORMATION CONCERNING THE BASIS SET WITH
C  SPHERICAL D AND F FUNCTIONS IS HOLD IN /BASSPH/.
C  ALL VARIABLES ARE IDENTIFIED BY THE IDNETIFIER 5
C  AT THE END
C
CEND            
C
      IMPLICIT INTEGER(A-Z)
      LOGICAL IDFGH,SCF,NONHF

      DIMENSION IAOPOP(8)

      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/BASSPH/NBAS5(8),NBASIS5,NBASSQ5,NBASTT5
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/LSYMSPH/NLENQ5(8),NLENT5(8)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/DFGH/IDFGH
      dimension ndrpop(8),ndrvrt(8) 
C
c---------------------------------------------------------------
c---   for drop-mo  cases   -----  Mar, 94,   ---  KB   --------
c---------------------------------------------------------------
c---   Reset  pop, vrt, nt, nf1, and nf2 -----------------------
c---------------------------------------------------------------
      call getrec(20,'JOBARC','NUMDROPA',1,NDROP)
      if (ndrop.ne.0.and..not.scf) then
       call getrec(20,'JOBARC','NDROPPOP',NIRREP,NDRPOP) 
       call getrec(20,'JOBARC','NDROPVRT',NIRREP,NDRVRT) 
       do 5 i=1,nirrep
        vrt(i,1) = vrt(i,1) + ndrvrt(i) 
        vrt(i,2) = vrt(i,2) + ndrvrt(i) 
        pop(i,1) = pop(i,1) + ndrpop(i) 
        pop(i,2) = pop(i,2) + ndrpop(i) 
  5    continue 
       NT(1) = 0
       NF1(1) = 0
       NF2(1) = 0
       DO 10 I=1,NIRREP
        NT(1) = POP(I,1) * VRT(I,1) + NT(1) 
        NF1(1) = POP(I,1) * POP(I,1) + NF1(1) 
        NF2(1) = VRT(I,1) * VRT(I,1) + NF2(1) 
  10   CONTINUE 
c
       IF (IUHF.EQ.0) THEN 
        NT(2) = NT(1) 
        NF1(2) = NF1(1) 
        NF2(2) = NF2(1) 
       ELSE
c
       NT(2) = 0
       NF1(2) = 0
       NF2(2) = 0
       DO 15 I=1,NIRREP
        NT(2) = POP(I,2) * VRT(I,2) + NT(2) 
        NF1(2) = POP(I,2) * POP(I,2) + NF1(2) 
        NF2(2) = VRT(I,2) * VRT(I,2) + NF2(2) 
  15   CONTINUE 
       ENDIF
c
      endif 
c---------------------------------------------------------------
      IF(IDFGH) THEN
C
        IF (NDROP .NE. 0) THEN
            CALL GETREC(20,'JOBARC','NUMBASI0',NIRREP,IAOPOP)
        ELSE
            CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,IAOPOP)
        ENDIF
C
        DO 90 IRREP=1,NIRREP

           NBAS5(IRREP) = IAOPOP(IRREP)
C
90     CONTINUE
C
       CALL DETBAS(NBAS)
C
      ELSE
C
       DO 100 IRREP=1,NIRREP
C
        NBAS(IRREP)=POP(IRREP,1)+VRT(IRREP,1)
        NBAS5(IRREP)=POP(IRREP,1)+VRT(IRREP,1)
C
100    CONTINUE
C
      ENDIF
C
      NBASIS=0
      NBASSQ=0
      NBASTT=0
      NBASIS5=0
      NBASSQ5=0
      NBASTT5=0
C
      DO 200 IRREP=1,NIRREP
C
       NBASIS=NBASIS+NBAS(IRREP)
C
       NBASSQ=NBASSQ+NBAS(IRREP)*NBAS(IRREP)
C
       NBASTT=NBASTT+(NBAS(IRREP)*(NBAS(IRREP)+1))/2
C
       NBASIS5=NBASIS5+NBAS5(IRREP)
C
       NBASSQ5=NBASSQ5+NBAS5(IRREP)*NBAS5(IRREP)
C
       NBASTT5=NBASTT5+(NBAS5(IRREP)*(NBAS5(IRREP)+1))/2
C
200   CONTINUE
C
      NLENQ(1)=NBASSQ
      NLENT(1)=NBASTT
      NLENQ5(1)=NBASSQ5
      NLENT5(1)=NBASTT5
C
      DO 300 IRREP=2,NIRREP
C
       NLENQ(IRREP)=0
       NLENT(IRREP)=0
       NLENQ5(IRREP)=0
       NLENT5(IRREP)=0
C
       DO 400 IRREP1=1,NIRREP
C
        IRREP2=DIRPRD(IRREP,IRREP1)
C
        IF(IRREP2.GT.IRREP1) THEN
C
         NLENT(IRREP)=NLENT(IRREP)+NBAS(IRREP1)*NBAS(IRREP2)
         NLENT5(IRREP)=NLENT5(IRREP)+NBAS5(IRREP1)*NBAS5(IRREP2)
C
        ENDIF
C
        NLENQ(IRREP)=NLENQ(IRREP)+NBAS(IRREP1)*NBAS(IRREP2)
        NLENQ5(IRREP)=NLENQ5(IRREP)+NBAS5(IRREP1)*NBAS5(IRREP2)
C
400    CONTINUE
C
300   CONTINUE
C
      CALL PUTREC(20,'JOBARC','SYMMLENG',16,NLENQ)
C
      RETURN
      END
