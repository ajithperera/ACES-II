      SUBROUTINE TPDABCI4(ICORE,MAXCOR,IUHF,LISTL2,LISTR2)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TWO,ONEM
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2),I0G(2)
C
C CALCULATION OF THE FOURTH ABCI CONTRIBUTION TO
C THE EOM-CCSD TWO-PARTICLE DENSITY MATRIX
C
C G(AB,CI) = 1/4 P(AB) [R(MN,EA)*L(EC,MN) 
C
C          = P(AB) * X(AC)*T(I,B) 
C
C  WHERE THE X ARE CALCULATED BY GFORMG.
C
C ALSO, ONLY THE SUM OF G(AB,CI) AND G(AI,BC) IS STORED ON DISK
C
CEND 
C
C CODED JG SEPTEMBER/93
C
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA TWO,ONEM /2.0D0,-1.0D0/
C
C FOR IRREPX=1, PUT R(MN,EA) ON LISTS 61-63.
C
      IF(IRREPX.EQ.1)THEN
       DO 10 ISPIN=3,3-2*IUHF,-1
        LISTT=43+ISPIN
c        LISTR=460+ISPIN
        LISTR=LISTR2-1+ISPIN
        LISTZ=60+ISPIN
        ISIZE=ISYMSZ(ISYTYP(1,LISTT),ISYTYP(2,LISTT))
        I000=1
        I010=I000+ISIZE*IINTFP
        IEND=I010+ISIZE*IINTFP
        CALL GETALL(ICORE(I010),ISIZE,1,LISTR)
        CALL PUTALL(ICORE(I010),ISIZE,1,LISTZ)
10     CONTINUE
       LISTR=61
      ELSE
c       LISTR=461
       LISTR=LISTR2
      ENDIF
C
C CALL GFORMG TO GET THE X INTERMEDIATES
C
c      CALL GFORMG(IRREPX,IRREPX,444,LISTR,0,ICORE,MAXCOR,2,ONEM,IUHF)      
      CALL GFORMG(IRREPX,IRREPX,LISTL2,LISTR,0,ICORE,MAXCOR,2,ONEM,IUHF)

C
C READ IN T1 AND X1 AMPLITUDES AND TRANSPOSE THE FORMER TO I,A 
C ORDERING
C
      I0T(1)=1
      I0T(2)=I0T(1)+IRPDPD(1,9)*IINTFP*IUHF
      I0G(1)=I0T(2)+IRPDPD(1,10)*IINTFP
      I0G(2)=I0G(1)+IRPDPD(1,19)*IINTFP*IUHF
      ISTART=I0G(2)+IRPDPD(1,20)*IINTFP
      I000=ISTART
      MXCOR=MAXCOR-ISTART+1
C
      CALL GETLST(ICORE(ISTART),1,1,1,1,90)
      CALL SYMTRA(1,VRT(1,1),POP(1,1),1,ICORE(ISTART),ICORE(I0T(1)))
      CALL GETLST(ICORE(I0G(1)),1,1,1,3,92)
      IF(IUHF.EQ.1) THEN
       CALL GETLST(ICORE(ISTART),1,1,1,2,90)
       CALL SYMTRA(1,VRT(1,2),POP(1,2),1,ICORE(ISTART),ICORE(I0T(2)))
       CALL GETLST(ICORE(I0G(2)),1,1,1,4,92)
      ENDIF
C
C FORM G*R TERMS
C
       IF(IUHF.EQ.1) THEN
C
C AAAA AND BBBB SPIN CASES
C
        DO 2000 ISPIN=1,2
         LISTG=126+ISPIN
         DO 2100 IRREP=1,NIRREP
          NUMDSG=IRPDPD(IRREP,ISYTYP(2,LISTG))
          DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
          NINCOR=MXCOR/(IINTFP*MAX(1,DISSYG))
          NLEFT =NUMDSG
          IFIRST=1
1         NREAD=MIN(NINCOR,NLEFT)
          CALL GETLST(ICORE(I000),IFIRST,NREAD,1,IRREP,LISTG)
          CALL GTTAU2(ICORE(I000),ICORE(I0G(ISPIN)),ICORE(I0T(ISPIN)),
     &                DISSYG,NUMDSG,VRT(1,ISPIN),POP(1,ISPIN),
     &                VRT(1,ISPIN),1,IRREP,ISPIN,IFIRST,NREAD)
          CALL PUTLST(ICORE(I000),IFIRST,NREAD,1,IRREP,LISTG)
          IFIRST=IFIRST+NREAD
          NLEFT=NLEFT-NREAD
          IF(NLEFT.NE.0)GOTO 1
2100     CONTINUE
2000    CONTINUE
       ENDIF
C
       DO 3000 ISPIN=1,1+IUHF
        LISTG=131-ISPIN
        DO 3100 IRREP=1,NIRREP
         NUMDSG=IRPDPD(IRREP,ISYTYP(2,LISTG))
         DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
         NINCOR=MXCOR/(IINTFP*MAX(1,DISSYG))
         NLEFT =NUMDSG
         IFIRST=1
2        NREAD=MIN(NINCOR,NLEFT)
         CALL GETLST(ICORE(I000),IFIRST,NREAD,1,IRREP,LISTG)
         CALL GTTAU2(ICORE(I000),ICORE(I0G(ISPIN)),
     &               ICORE(I0T(3-ISPIN)),DISSYG,NUMDSG,
     &               VRT(1,ISPIN),POP(1,3-ISPIN),VRT(1,3-ISPIN),
     &               1,IRREP,2+ISPIN,IFIRST,NREAD)
         CALL PUTLST(ICORE(I000),IFIRST,NREAD,1,IRREP,LISTG)
         IFIRST=IFIRST+NREAD
         NLEFT=NLEFT-NREAD
         IF(NLEFT.NE.0)GOTO 2
3100    CONTINUE
3000   CONTINUE
C
C ALL DONE, RETURN
C
      TWO=0.D0
      if(iuhf.eq.0) then
       call checkgam1(icore,30,130,two,iuhf,2,vrt)
      endif
      IF(IUHF.EQ.1) THEN
       CALL CHECKGAM(ICORE,30,130,TWO)
       CALL CHECKGAM(ICORE,27,127,TWO)
       CALL CHECKGAM(ICORE,28,128,TWO)
       CALL CHECKGAM(ICORE,29,129,TWO)
      ENDIF
C
      RETURN
      END
