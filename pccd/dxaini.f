      SUBROUTINE DXAINI(DOO,DVO,DVV,IOO,ICORE,MAXCOR,IUHF,
     &                  NONHF_TERMS_EXIST)
C
C DRIVER FOR THE FOLLOWING CONTRIBUTION TO THE I INTERMEDIATE:
C
C       I(ij) = - (1/2) SUM D(pq) [<pi||qj> + <pj||qi>]
C                       p,q
C
CEND
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION DOO,DVO,DVV,IOO
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,
     &        GABCD,RELAXED,TRULY_NONHF,NONHF_TERMS_EXIST
      DIMENSION ICORE(MAXCOR),DOO(1),DVO(1),DVV(1),IOO(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPP(255,2),DIRPRD(8,8)
      COMMON /DERIV/ DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &               TRIP2,GABCD,RELAXED,TRULY_NONHF
CJDW KKB stuff
C     COMMON /SYM/  POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
CJDW END
C
      CALL OOINOO(IOO,DOO,ICORE,MAXCOR,IUHF)
CSSS      call checksum("IOO-1",Ioo,nfmi(1))
      CALL VVINOO(IOO,DVV,ICORE,MAXCOR,IUHF)
CSSS      call checksum("IOO-2",Ioo,nfmi(1))
      IF(QRHF.OR.ROHF .OR. NONHF_TERMS_EXIST)THEN
       I0AA=1
       I0BB=I0AA+NT(1)*IINTFP
       CALL SYMTRA(1,VRT(1,1),POP(1,1),1,DVO,ICORE(I0AA))
       IF(IUHF.NE.0)THEN
        IBOT=I0BB+NT(2)*IINTFP
        IOFF=1+NT(1)
        CALL SYMTRA(1,VRT(1,2),POP(1,2),1,DVO(IOFF),ICORE(I0BB))
       ELSE
        IBOT=I0BB
       ENDIF
       MXCOR=MAXCOR-IBOT+1
       CALL VOINOO(IOO,ICORE,ICORE(IBOT),MXCOR,IUHF)
CSSS       call checksum("IOO-3",Ioo,nfmi(1))
      ENDIF
      RETURN
      END