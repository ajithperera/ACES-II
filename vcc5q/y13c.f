      SUBROUTINE Y13C(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZD,DISSIZT,DISSIZY
      INTEGER POP,VRT,DIRPRD
      DIMENSION CORE(1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I * (I-1)/2
C
      WRITE(6,1000)
 1000 FORMAT(' @Y13C-I, T2 * Y2 contributions. ')
C
C     Y13(C,E,K,M) = + SUM  T(CE,KJ) * Y2(J,M)
C                       J
C
C     CEKM = CEKJ * JM
C     cekm = cekj * jm
C
      DO  100 ISPIN=1,IUHF+1
C
      I000 = 1
      I010 = I000 + NFMI(ISPIN)
      CALL GETLST(CORE(I000),1,1,1,ISPIN,91)
C
      DO   40 IRPKM=1,NIRREP
C
      IRPCE = IRPKM
      IRPKJ = IRPCE
C
      IF(IRPDPD(IRPCE,18+ISPIN).EQ.0.OR.
     1   IRPDPD(IRPKM,20+ISPIN).EQ.0) GOTO 40
C
      I020 = I010 + IRPDPD(IRPCE,18+ISPIN) * IRPDPD(IRPKJ,20+ISPIN)
      I030 = I020 + IRPDPD(IRPCE,18+ISPIN) * IRPDPD(IRPKM,20+ISPIN)
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
 1020 FORMAT(' @Y13C-I, Insufficient memory. Need ',I15,' . Got ',I15)
      STOP 'Y13C'
      ENDIF
C
      DISSIZY = IRPDPD(IRPCE,18+ISPIN)
      NDISY   = IRPDPD(IRPKM,20+ISPIN)
C
C     T(C<E,K<J) or T(c<e,k<j)
      CALL GETLST(CORE(I010),1,IRPDPD(IRPKJ,2+ISPIN),1,
     1            IRPKJ,43+ISPIN)
C     (C<E,K<J) ---> (C,E,K<J)
      CALL SYMEXP2(IRPCE,VRT(1,ISPIN),DISSIZY,
     1             IRPDPD(IRPCE,  ISPIN),
     1             IRPDPD(IRPKJ,2+ISPIN),CORE(I010),CORE(I010))
C     (C,E,K<J) ---> (C,E,K,J)
      CALL  SYMEXP(IRPKJ,POP(1,ISPIN),DISSIZY,CORE(I010))
C
      NSIZY   = DISSIZY * NDISY
C
      CALL ZERO(CORE(I020),NSIZY)
      DO   30 IRPM=1,NIRREP
C
      IRPK = DIRPRD(IRPKM,IRPM)
      IRPJ = IRPM
C
      IF(POP(IRPK,ISPIN).EQ.0.OR.POP(IRPM,ISPIN).EQ.0) GOTO 30
C
      IOFFD = I000 +                        IOFFOO(IRPM,1,2+ISPIN)
      IOFFT = I010 + IRPDPD(IRPCE,18+ISPIN)*IOFFOO(IRPJ,IRPKJ,2+ISPIN)
      IOFFY = I020 + IRPDPD(IRPCE,18+ISPIN)*IOFFOO(IRPM,IRPKM,2+ISPIN)
C
      DISSIZT = IRPDPD(IRPCE,18+ISPIN) * POP(IRPK,ISPIN)
      NDIST   = POP(IRPJ,ISPIN)
      DISSIZD = POP(IRPJ,ISPIN)
      NDISD   = POP(IRPM,ISPIN)
C
      CALL XGEMM('N','N',DISSIZT,NDISD,NDIST,
     1             1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFD),DISSIZD,1.0D+00,
     1            CORE(IOFFY),DISSIZT)
   30 CONTINUE
C
      CALL GETLIST(CORE(I010),1,NDISY,2,IRPKM,18+ISPIN)
c     CALL    VADD(CORE(I010),CORE(I010),CORE(I020),NSIZY,-1.0D+00)
C     We have KJ !
      CALL    VADD(CORE(I010),CORE(I010),CORE(I020),NSIZY, 1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISY,2,IRPKM,18+ISPIN)
   40 CONTINUE
  100 CONTINUE
C
C     CeKm = CeKj * jm
C     cEkM = cEkJ * JM
C
      DO  200 ISPIN=1,IUHF+1
C
      I000 = 1
      I010 = I000 + NFMI(3-ISPIN-1+IUHF)
      CALL GETLST(CORE(I000),1,1,1,3-ISPIN-1+IUHF,91)
C
      DO  140 IRPKM=1,NIRREP
C
      IRPCE = IRPKM
      IRPKJ = IRPCE
C
      IF(IRPDPD(IRPCE,13).EQ.0.OR.
     1   IRPDPD(IRPKM,14).EQ.0) GOTO 140
C
      I020 = I010 + IRPDPD(IRPCE,13) * IRPDPD(IRPKJ,14)
      I030 = I020 + IRPDPD(IRPCE,13) * IRPDPD(IRPKM,14)
      I040 = I030 + MAX(IRPDPD(IRPCE,13),IRPDPD(IRPKJ,14))
      I050 = I040 + MAX(IRPDPD(IRPCE,13),IRPDPD(IRPKJ,14))
      I060 = I050 + MAX(IRPDPD(IRPCE,13),IRPDPD(IRPKJ,14))
      NEED = IINTFP * I060
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y13C'
      ENDIF
C
      DISSIZY = IRPDPD(IRPCE,13)
      NDISY   = IRPDPD(IRPKM,14)
C
C     T(C,e,K,j) or T(E,c,J,k)
      CALL GETLST(CORE(I010),1,IRPDPD(IRPKJ,14),1,IRPKJ,46)
C
      IF(ISPIN.EQ.2)THEN
C     T(E,c,J,k) ---> T(c,E,J,k)
      CALL SYMTR3(IRPKJ,VRT(1,1),VRT(1,2),IRPDPD(IRPKJ,13),
     1            IRPDPD(IRPKJ,14),CORE(I010),
     1            CORE(I030),CORE(I040),CORE(I050))
C     T(c,E,J,k) ---> T(c,E,k,J)
      CALL SYMTR1(IRPKJ,POP(1,1),POP(1,2),IRPDPD(IRPKJ,13),
     1            CORE(I010),CORE(I030),CORE(I040),CORE(I050))
      ENDIF
C
      NSIZY   = DISSIZY * NDISY
C
      CALL ZERO(CORE(I020),NSIZY)
      DO  130 IRPM=1,NIRREP
C
      IRPK = DIRPRD(IRPKM,IRPM)
      IRPJ = IRPM
C
      IF(POP(IRPK,ISPIN).EQ.0.OR.POP(IRPM,3-ISPIN).EQ.0) GOTO 130
C
      IOFFD = I000 +                  IOFFOO(IRPM,1,5-ISPIN)
      IOFFT = I010 + IRPDPD(IRPCE,13)*IOFFOO(IRPJ,IRPKJ,4+ISPIN)
      IOFFY = I020 + IRPDPD(IRPCE,13)*IOFFOO(IRPM,IRPKM,4+ISPIN)
C
      DISSIZT = IRPDPD(IRPCE,13) * POP(IRPK,ISPIN)
      NDIST   = POP(IRPJ,3-ISPIN)
      DISSIZD = POP(IRPJ,3-ISPIN)
      NDISD   = POP(IRPM,3-ISPIN)
C
      CALL XGEMM('N','N',DISSIZT,NDISD,NDIST,
     1             1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFD),DISSIZD,1.0D+00,
     1            CORE(IOFFY),DISSIZT)
  130 CONTINUE
C
      CALL GETLIST(CORE(I010),1,NDISY,2,IRPKM,20+ISPIN)
c     CALL    VADD(CORE(I010),CORE(I010),CORE(I020),NSIZY,-1.0D+00)
      CALL    VADD(CORE(I010),CORE(I010),CORE(I020),NSIZY, 1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISY,2,IRPKM,20+ISPIN)
  140 CONTINUE
  200 CONTINUE
C
C     CekM = CekJ * JM
C     cEKm = cEKj * jm
C
      DO  300 ISPIN=1,IUHF+1
C
      I000 = 1
      I010 = I000 + NFMI(ISPIN)
      CALL GETLST(CORE(I000),1,1,1,ISPIN,91)
C
      DO  240 IRPKM=1,NIRREP
C
      IRPCE = IRPKM
      IRPKJ = IRPCE
C
      IF(IRPDPD(IRPCE,13).EQ.0.OR.
     1   IRPDPD(IRPKM,14).EQ.0) GOTO 240
C
      I020 = I010 + IRPDPD(IRPCE,13) * IRPDPD(IRPKJ,14)
      I030 = I020 + IRPDPD(IRPCE,13) * IRPDPD(IRPKM,14)
      I040 = I030 + MAX(IRPDPD(IRPCE,13),IRPDPD(IRPKJ,14))
      I050 = I040 + MAX(IRPDPD(IRPCE,13),IRPDPD(IRPKJ,14))
      I060 = I050 + MAX(IRPDPD(IRPCE,13),IRPDPD(IRPKJ,14))
      NEED = IINTFP * I060
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y13C'
      ENDIF
C
      DISSIZY = IRPDPD(IRPCE,13)
      NDISY   = IRPDPD(IRPKM,14)
C
C     T(C,e,J,k) or T(E,c,K,j)
      CALL GETLST(CORE(I010),1,IRPDPD(IRPKJ,14),1,IRPKJ,46)
C
      IF(ISPIN.EQ.1)THEN
C     T(C,e,J,k) ---> T(C,e,k,J)
      CALL SYMTR1(IRPKJ,POP(1,1),POP(1,2),IRPDPD(IRPKJ,13),
     1            CORE(I010),CORE(I030),CORE(I040),CORE(I050))
      ELSE
C     T(E,c,K,j) ---> T(c,E,K,j)
      CALL SYMTR3(IRPKJ,VRT(1,1),VRT(1,2),IRPDPD(IRPKJ,13),
     1            IRPDPD(IRPKJ,14),CORE(I010),
     1            CORE(I030),CORE(I040),CORE(I050))
      ENDIF
C
      NSIZY   = DISSIZY * NDISY
C
      CALL ZERO(CORE(I020),NSIZY)
      DO  230 IRPM=1,NIRREP
C
      IRPK = DIRPRD(IRPKM,IRPM)
      IRPJ = IRPM
C
      IF(POP(IRPK,3-ISPIN).EQ.0.OR.POP(IRPM,ISPIN).EQ.0) GOTO 230
C
      IOFFD = I000 +                  IOFFOO(IRPM,    1,2+ISPIN)
      IOFFT = I010 + IRPDPD(IRPCE,13)*IOFFOO(IRPJ,IRPKJ,7-ISPIN)
      IOFFY = I020 + IRPDPD(IRPCE,13)*IOFFOO(IRPM,IRPKM,7-ISPIN)
C
      DISSIZT = IRPDPD(IRPCE,13) * POP(IRPK,3-ISPIN)
      NDIST   = POP(IRPJ,ISPIN)
      DISSIZD = POP(IRPJ,ISPIN)
      NDISD   = POP(IRPM,ISPIN)
C
      CALL XGEMM('N','N',DISSIZT,NDISD,NDIST,
c    1            -1.0D+00,
     1             1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFD),DISSIZD,1.0D+00,
     1            CORE(IOFFY),DISSIZT)
  230 CONTINUE
C
      CALL GETLIST(CORE(I010),1,NDISY,2,IRPKM,22+ISPIN)
      CALL    VADD(CORE(I010),CORE(I010),CORE(I020),NSIZY,-1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISY,2,IRPKM,22+ISPIN)
  240 CONTINUE
  300 CONTINUE
      RETURN
      END