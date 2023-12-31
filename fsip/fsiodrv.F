   
      SUBROUTINE FSIODRV(ICORE,MAXCOR,IUHF,IDO,IACTIVE,
     &                   NACTIVE,NBAS,I0)
C
C DRIVER FOR CONSTRUCTION OF THE GATHER-SCATTER I/O VECTORS 
C NEEDED FOR CALCULATIONS INVOLVING ACTIVE OR INACTIVE INDICES
C
      IMPLICIT INTEGER (A-Z)
#include "maxbasfn.par"
      DIMENSION ICORE(MAXCOR),IDO(22),SPIN(2,22),NACTIVE(2)
      CHARACTER*3 TYPE(2,22),LFTTYP,RHTTYP
      CHARACTER*4 PACK(22),TYPEX
      DIMENSION IACTIVE(MAXBASFN,2),IVECTOR(255,3,2),POPLFT(8),POPRHT(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM   / POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
C COMMON BLOCK WITH POPULATIONS FOR VARIOUS INDEX TYPES:
C
C
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
C
C POSITIONS IN CORE WHERE THE GATHER-SCATTER VECTORS ARE LOCATED.
C
C        LEFT INDEX    - IRREP
C        CENTRAL INDEX - SYMMETRY TYPE
C        RIGHT INDEX   - ORBITAL TYPE (AI, IA, ETC)
C
      COMMON /FSIOLOC/ IFSPOS(8,22,8),ISCRLC
C
      DATA TYPE /'VRT','VRT' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','OCC' , 'VRT','VRT' , 'VRT','VRT' ,
     &           'OCC','OCC' , 'OCC','OCC' , 'VRT','OCC' ,
     &           'VRT','OCC' , 'VRT','OCC' , 'VRT','OCC' ,
     &           'VRT','VRT' , 'OCC','OCC' , 'VRT','VRT' ,
     &           'OCC','VRT' , 'OCC','VRT' , 'OCC','VRT' ,
     &           'VRT','VRT' , 'VRT','VRT' , 'OCC','OCC' ,
     &           'OCC','OCC'/
      DATA SPIN / 1,1 , 2,2 , 1,1 , 2,2 , 1,1 , 2,2 , 1,1 , 2,2 , 
     &            1,1 , 2,2 ,
     &            1,2 , 2,1 , 1,2 , 1,2 , 1,2 , 1,1 , 2,2 , 1,2 ,
     &            1,1 , 2,2 ,
     &            1,1 , 2,2/
      DATA PACK /'PACK','PACK','PACK','PACK','PCK2',
     &           'PCK2','PCK2','PCK2','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL',
     &           'FULL','FULL','FULL','FULL','FULL',
     &           'FULL','FULL'/
C
      IOFF=I0
C
C FILL VECTOR OF LENGTH NBAS WITH "1" 
C
      CALL IZERO(IFSPOS,8*22*8)
      CALL IZERO(IVECTOR,255*6)
      DO 5 ISPIN=1,2
       DO 10 ITYPE=1,2
        DO 20 IPOS=1,NBAS
         IVECTOR(IPOS,ITYPE,ISPIN)=1
20      CONTINUE
10     CONTINUE
5     CONTINUE
C
C NOW ZERO OUT POSITIONS WHICH CORRESPOND TO ACTIVE ORBITALS (ITYPE=2)
C AND POSITIONS WHICH CORRESPOND TO INACTIVE ORBITALS (ITYPE=3)
C
       DO 25 ISPIN=1,2
        DO 28 I=1,NACTIVE(ISPIN)
         IPOS=IACTIVE(I,ISPIN) 
         IVECTOR(IPOS,2,ISPIN)=0
         IVECTOR(IPOS,3,ISPIN)=1
28      CONTINUE
25     CONTINUE
C
       CALL FSINIPOP(IVECTOR,IUHF)
C
C NOW GO THROUGH THE DIFFERENT DISTRIBUTION TYPES ONE BY ONE
C
      DO 30 ITYPE=1,22
       IF(IDO(ITYPE).EQ.1)THEN
        LFTTYP=TYPE(1,ITYPE)
        LFTSPN=SPIN(1,ITYPE) 
        RHTTYP=TYPE(2,ITYPE)
        RHTSPN=SPIN(2,ITYPE)
        TYPEX=PACK(ITYPE)
        IF(LFTTYP.EQ.'OCC')CALL ICOPY(POP(1,LFTSPN),POPLFT,8)
        IF(LFTTYP.EQ.'VRT')CALL ICOPY(VRT(1,LFTSPN),POPLFT,8)
        IF(RHTTYP.EQ.'OCC')CALL ICOPY(POP(1,RHTSPN),POPRHT,8)
        IF(RHTTYP.EQ.'VRT')CALL ICOPY(VRT(1,RHTSPN),POPRHT,8)
C
C LEFT INDEX ACTIVE, RIGHT INDEX FULL
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,3,LFTSPN),IVECTOR(1,1,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDAN(1,ITYPE),
     &             IFSPOS(1,ITYPE,1),IOFF)
C
C LEFT INDEX FULL  , RIGHT INDEX ACTIVE
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,1,LFTSPN),IVECTOR(1,3,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDNA(1,ITYPE),
     &             IFSPOS(1,ITYPE,2),IOFF)
C
C LEFT INDEX ACTIVE  , RIGHT INDEX ACTIVE
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,3,LFTSPN),IVECTOR(1,3,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDAA(1,ITYPE),
     &             IFSPOS(1,ITYPE,3),IOFF)
C
C LEFT INDEX INACTIVE, RIGHT INDEX FULL
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,2,LFTSPN),IVECTOR(1,1,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDIN(1,ITYPE),
     &             IFSPOS(1,ITYPE,4),IOFF)
C
C LEFT INDEX FULL , RIGHT INDEX INACTIVE
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,1,LFTSPN),IVECTOR(1,2,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDNI(1,ITYPE),
     &             IFSPOS(1,ITYPE,5),IOFF)
C
C LEFT INDEX INACTIVE, RIGHT INDEX INACTIVE
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,2,LFTSPN),IVECTOR(1,2,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDII(1,ITYPE),
     &             IFSPOS(1,ITYPE,6),IOFF)
C
C LEFT INDEX ACTIVE, RIGHT INDEX INACTIVE
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,3,LFTSPN),IVECTOR(1,2,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDAI(1,ITYPE),
     &             IFSPOS(1,ITYPE,7),IOFF)
C
C LEFT INDEX INACTIVE, RIGHT INDEX ACTIVE
C
        CALL FSMAP(LFTTYP,RHTTYP,LFTSPN,RHTSPN,TYPEX,POPLFT,POPRHT,
     &             IVECTOR(1,2,LFTSPN),IVECTOR(1,3,RHTSPN),NBAS,
     &             ISZFUL,ISZOUT,ICORE(IOFF),FSDPDIA(1,ITYPE),
     &             IFSPOS(1,ITYPE,8),IOFF)
       ENDIF
30    CONTINUE
C
      WRITE(6,1000)IOFF-I0
1000  FORMAT(T3,'@FSIODRV-I, I/O vectors require ',I8,' words.')
C
C NOW ALLOCATE SPACE FOR THE LARGEST DISTRIBUTION.  THIS IS USED
C AS A SCRATCH AREA IN FSGET. 
C
      IDIFF=IOFF-I0
      IF(MOD(IDIFF,2).NE.0)IOFF=IOFF+1
      ISCRLC=IOFF
      DO 100 I=1,22
       DO 101 J=1,8
        DSZMAX=MAX(DSZMAX,IRPDPD(J,I))
101    CONTINUE
100   CONTINUE
      I0=ISCRLC+DSZMAX*IINTFP
      RETURN
      END
