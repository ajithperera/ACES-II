      SUBROUTINE DS25AB(ICORE,MAXCOR,ITYPE,IUHF,NMO,IRRVEC)
      IMPLICIT INTEGER(A-Z)
      INTEGER DIRPRD
      LOGICAL REDUCE
      CHARACTER*4 INTYPE(6)
      CHARACTER*2 SPCASE(3)
C 
      DIMENSION ICORE(MAXCOR),IRRVEC(NMO,2)
      DIMENSION IPW(9),IPDSZ(8),IPDIS(8)
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /BUFLEN/ ILNBUF
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/MINTPC/REDUCE
C
      DATA INTYPE /'HHHH','PHHH','PPHH','PHPH','PPPH','PPPP'/
      DATA SPCASE /'AA','BB','AB'/
C
      ISPIN=2
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      IF(ITYPE.EQ.2) THEN
       NWDSZ=NOCCA*NOCCB
      ELSE
       NWDSZ=NVRTA*NVRTB
      ENDIF
C
C ONLY NEED ONE CALL IF WE ARE GENERATING THE AB INTEGRALS FOR AN RHF
C  REFERENCE.  MUCH OF WHAT IS PASSED IS IGNORED IN THAT INSTANCE.
C
      DO 10 I=1,IUHF+1
       ISPIN=ISPIN+1
       IF(ITYPE.EQ.2)THEN
        ISW1=14
        IF(I.EQ.1) THEN
         NWDIS=NVRTA*NOCCB
         ISW2=11
         IF(IUHF.EQ.0) ISW2=18
        ELSE
         NWDIS=NVRTB*NOCCA
         ISW2=18
        ENDIf
       ELSE IF(ITYPE.EQ.5)THEN
        ISW1=13
        IF(I.EQ.1) THEN
         NWDIS=NVRTB*NOCCA
         ISW2=18
         IF(IUHF.EQ.0) ISW2=11
        ELSE
         NWDIS=NVRTA*NOCCB
         ISW2=11
        ENDIF
       ENDIF
       CALL SETUP3(ICORE,MAXCOR,IR1,IR2,IR3,ISW1,ISW2,
     &             NWDIS,NWDSZ)
C
C   CALCULATE OF W-ARRAY AND OFFSETS IN SYMMETRY PACKED W-ARRAY
C 
       IPW(1)=0
       IPDSZ(1)=0
       IPDIS(1)=0
       DO 2 IRREP=1,NIRREP-1
       IPW(IRREP+1)=IPW(IRREP)+ICORE(IR1+IRREP-1)*
     &              ICORE(IR1+NIRREP+IRREP-1)
       IPDSZ(IRREP+1)=IPDSZ(IRREP)+ICORE(IR1+IRREP+NIRREP-1)
       IPDIS(IRREP+1)=IPDIS(IRREP)+ICORE(IR1+IRREP-1)
2      CONTINUE
       NSIZE=IPW(NIRREP)+ICORE(IR1+NIRREP-1)*ICORE(IR1+2*NIRREP-1)
      
       IR4=IR3+NSIZE*IINTFP
       IR5=IR4+ILNBUF*IINTFP
       IR6=IR5+ILNBUF
       IF(IR6.GT.MAXCOR)THEN
        IF(ITYPE.EQ.2)THEN
         WRITE(6,3333)
3333     FORMAT(T3,'@DS25AB-F, Insufficient core to sort ijka ',
     &             'MO integrals!')
         CALL INSMEM('DS25AA',MAXCOR,ITEST)
        ENDIF
        IF(IFLAGS(1).GE.1)THEN
         WRITE(LUOUT,1005)INTYPE(ITYPE),SPCASE(ISPIN)
        ENDIF
        MEMLFT=MAXCOR-IR3
        IF(REDUCE) MEMLFT=MEMLFT/2
        IF(IUHF.EQ.0) THEN
C
         NVRT=NVRTA
         NOCC=NOCCA
         ISPIN=1
C
C  ONE PROBLEM HERE, WE NEED A DIFFERENT SYMMETRy VECTOR FOR
C  THE OUT CORE VERSION (I,J FULL ARRAY)
C  STARTING POSITION OF THIS SYMMETRY VECTOR
C
         IR2A=IR2+NWDIS
         IR1A=IR1+NIRREP
         DO 5 IRREP=1,NIRREP
          ICORE(IR1A+IRREP-1)=0
5        CONTINUE
         IND=0
         DO 6 IA=1,NVRT
         DO 6 IB=1,NVRT
          IND=IND+1
          IRREPAB=DIRPRD(IRRVEC(IA+NOCC,ISPIN),IRRVEC(IB+NOCC,ISPIN))
          ICORE(IR1A+IRREPAB-1)=ICORE(IR1A+IRREPAB-1)+1
          ICORE(IR2A+IND-1)=ICORE(IR1A+IRREPAB-1)
6        CONTINUE
         NWDSZ=0
         DO 7 IRREP=1,NIRREP
          NWDSZ=NWDSZ+ICORE(IR1A+IRREP-1)
7        CONTINUE
         IR3=IR2+NWDIS+NWDSZ
         IF(MOD(IR3,2).EQ.0) IR3=IR3+1
         IPW(1)=0
         DO 8 IRREP=1,NIRREP-1
          IPW(IRREP+1)=IPW(IRREP)+ICORE(IR1+IRREP-1)*
     &                 ICORE(IR1+IRREP+NIRREP-1)
8        CONTINUE

         CALL DABCI0(ICORE(IR3),MEMLFT,ICORE(IR1),ICORE(IR2),
     &               IPW,IPDIS,NWDIS,1,NMO,IRRVEC)
        ELSE IF(I.EQ.2)THEN
         CALL DABCI1(ICORE(IR3),MEMLFT,ICORE(IR1),ICORE(IR2),
     &               IPW,IPDIS,IPDSZ,NWDIS,NMO,IRRVEC)
        ELSE IF(I.EQ.1) THEN
         CALL DABCI2(ICORE(IR3),MEMLFT,ICORE(IR1),ICORE(IR2),
     &               IPW,IPDIS,IPDSZ,NWDIS,NMO,IRRVEC)
        ENDIF
       ELSE
        IF(IFLAGS(1).GE.1)THEN
         WRITE(LUOUT,1002)INTYPE(ITYPE),SPCASE(ISPIN)
         WRITE(LUOUT,1001)IR6,MAXCOR
        ENDIF
        CALL ST25AB(ICORE(IR3),ICORE(IR4),ICORE(IR5),ICORE(IR1),
     &              ICORE(IR2),IPW,IPDIS,IPDSZ,NSIZE,NWDIS,NOCCA,
     &              NOCCB,NVRTA,NVRTB,ITYPE,ISPIN,IUHF,I,
     &              IRRVEC(1,1),IRRVEC(1,2))
       ENDIF
10    CONTINUE
1005  FORMAT(T3,'@DS25AB-I, Out-of-core sort used for ',A4,A2,
     &             ' integrals.')
1002  FORMAT(T3,'@DS25AB-I, ',A4,A2,' integrals will be sorted in ',
     & 'core.')
1001  FORMAT(T3,'   Words required: ',I8,'  Words available: ',I8)
      RETURN
      END