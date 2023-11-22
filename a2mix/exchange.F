      SUBROUTINE EXCHANGE(N,ILNBUF,BUF,IBUF,
     & IUHF,EEXCH,EHAR,ESCF,DENSA,DENSB,DENSAA,DENSBB,
     & CDENSA,CDENSB,REPULS,ISYMC1,IMAP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C     BUF - AO integral buffer; dimension is ILNBUF
C
C    IBUF - AO integral index buffer; dimension is ILNBUF
C
      PARAMETER(LUINT=10)
C
      INTEGER AND
      CHARACTER*80 FNAME
      CHARACTER*8 LABEL
      CHARACTER*32 JUNK
      LOGICAL YESNO
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION DENSA(N*N),DENSB(N*N),CDENSA(N*N),
     & CDENSB(N*N),DENSAA(N*(N+1)/2),DENSBB(N*(N+1)/2)
      DIMENSION IMAP(N)
C
      COMMON // ICORE(1)
      COMMON /ISTART/ I0
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /ISYMINF/ NIRREP,NSOIRP(8)
      COMMON /IPAR/ LUOUT
C
C  Statement functions for packing and unpacking indices.
C
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
C
      N2=N*N
C
C     Set up density matrix for rhf case
      IF(IUHF.EQ.0)THEN
      DO 810 I=1,N2
         DENSB(I)=DENSA(I)
  810 CONTINUE
C
      ENDIF
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The HF density in MO basis"
      Write(6,"(5(F12.6))")(Densa(i), i=1,N2)
#endif 
C
C Fill alpha density array for contraction with one electron integrals
      ICNT=0
      IOFF=0
      DO 710 ISYM=1,NIRREP
         DO 712 I=1,NSOIRP(ISYM)
            DO 715 J=1,I
               ICNT=ICNT+1
               DENSAA(ICNT)=DENSA((IOFF+I-1)*N+IOFF+J)
               DENSBB(ICNT)=DENSB((IOFF+I-1)*N+IOFF+J)
               IF(I.NE.J)THEN
                  DENSAA(ICNT)=2.D+00*DENSAA(ICNT)
                  DENSBB(ICNT)=2.D+00*DENSBB(ICNT)
               ENDIF
  715       CONTINUE
  712    CONTINUE
         IOFF=IOFF+NSOIRP(ISYM)
  710 CONTINUE
C
C     Multiply off diagonal elements of the alpha density matrix by two 
C     for coulomb contraction with two electron integrals
      DO 730 I=1,N
         DO 730 J=1,I-1
            CDENSA((I-1)*N+J)=2.D+00*DENSA((I-1)*N+J)
            CDENSB((I-1)*N+J)=2.D+00*DENSB((I-1)*N+J)
  730 CONTINUE
      DO 735 I=1,N
         DO 735 J=I+1,N
            CDENSA((I-1)*N+J)=2.D+00*DENSA((I-1)*N+J)
            CDENSB((I-1)*N+J)=2.D+00*DENSB((I-1)*N+J)
  735 CONTINUE
      DO 737 I=1,N
         CDENSA((I-1)*N+I)=DENSA((I-1)*N+I)
         CDENSB((I-1)*N+I)=DENSB((I-1)*N+I)
  737 CONTINUE
C
C     Read integrals and contract them with the appropriate density
C     matrix for one-electron, coulomb and exchange energies
C
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
      IF(YESNO)THEN
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',
     &        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
C
C     One electron energy
C
      E1E=0.D+00
      CALL LOCATE(LUINT,'ONEHAMIL')
    2 READ(LUINT)BUF,IBUF,NUT
      DO 5 INT=1,NUT
        E1E=E1E+DENSAA(IBUF(INT))*BUF(INT)
        E1E=E1E+DENSBB(IBUF(INT))*BUF(INT)
    5 CONTINUE
      IF(NUT.EQ.ILNBUF) GOTO 2
C
C     Two electron energy
C
      ECOULAB=0.D+00
      ECOULAA=0.D+00
      ECOULBB=0.D+00
      EEXCHAA=0.D+00
      EEXCHBB=0.D+00
C
      CALL LOCATE(LUINT,'TWOELSUP')
    1 READ(LUINT)BUF,IBUF,NUT
      DO 10 INT=1,NUT
        I=IUPKI(IBUF(INT))
        J=IUPKJ(IBUF(INT))
        K=IUPKK(IBUF(INT))
        L=IUPKL(IBUF(INT))
        X=BUF(INT)
C
C     alpha beta coulomb energy
           TMP1A=X*CDENSA((I-1)*N+J)*CDENSB((K-1)*N+L)
           IF(N*(I-1)+J.NE.N*(K-1)+L)
     &        TMP1A=TMP1A+X*CDENSA((K-1)*N+L)*CDENSB((I-1)*N+J)
C
C     alpha alpha coulomb energy
           TMP2A=X*CDENSA((I-1)*N+J)*CDENSA((K-1)*N+L)
           TMP2B=X*CDENSB((I-1)*N+J)*CDENSB((K-1)*N+L)
C
C     alpha alpha exchange energy
           TMP3A=X*DENSA((I-1)*N+L)*DENSA((J-1)*N+K)
           TMP3B=X*DENSB((I-1)*N+L)*DENSB((J-1)*N+K)
C
           IF(I.NE.J.OR.K.NE.L)THEN
              TMP3A=TMP3A+X*DENSA((I-1)*N+K)*DENSA((J-1)*N+L)
              TMP3B=TMP3B+X*DENSB((I-1)*N+K)*DENSB((J-1)*N+L)
           ENDIF
C
           IF(I.NE.J.AND.K.NE.L)THEN
              TMP3A=2.D+00*TMP3A
              TMP3B=2.D+00*TMP3B
           ENDIF
C
C     Off diagonal elements
           IF(N*(I-1)+J.NE.N*(K-1)+L)THEN
C     alpha alpha coulomb energy
              TMP2A=2.D+00*TMP2A
              TMP2B=2.D+00*TMP2B
C     alpha alpha exchange energy
              TMP3A=2.D+00*TMP3A
              TMP3B=2.D+00*TMP3B
           ENDIF
C
        ECOULAB=ECOULAB+TMP1A
        ECOULAA=ECOULAA+TMP2A
        ECOULBB=ECOULBB+TMP2B
        EEXCHAA=EEXCHAA+TMP3A
        EEXCHBB=EEXCHBB+TMP3B
C
   10 CONTINUE
      IF(NUT.EQ.ILNBUF) GOTO 1
      CLOSE(LUINT,STATUS='KEEP')
      ENDIF
C
C     If symmetry is C1 we are done
      IF(ISYMC1.EQ.1)GOTO 9999 
      CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
      IF(YESNO)THEN
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
C
C     Two electron energy
C
      CALL LOCATE(LUINT,'TWOELSUP')
 1001 READ(LUINT)BUF,IBUF,NUT
      DO 1010 INT=1,NUT
        I=IUPKI(IBUF(INT))
        J=IUPKJ(IBUF(INT))
        K=IUPKK(IBUF(INT))
        L=IUPKL(IBUF(INT))
        X=BUF(INT)
C
C     alpha alpha exchange energy
           TMP3A=X*DENSA((I-1)*N+L)*DENSA((J-1)*N+K)
           TMP3B=X*DENSB((I-1)*N+L)*DENSB((J-1)*N+K)
C
           IF(I.NE.J.OR.K.NE.L)THEN
              TMP3A=TMP3A+X*DENSA((I-1)*N+K)*DENSA((J-1)*N+L)
              TMP3B=TMP3B+X*DENSB((I-1)*N+K)*DENSB((J-1)*N+L)
           ENDIF
C
           IF(I.NE.J.AND.K.NE.L)THEN
              TMP3A=2.D+00*TMP3A
              TMP3B=2.D+00*TMP3B
           ENDIF
C
C     Off diagonal elements
           IF(N*(I-1)+J.NE.N*(K-1)+L)THEN
C     alpha alpha exchange energy
              TMP3A=2.D+00*TMP3A
              TMP3B=2.D+00*TMP3B
           ENDIF
C
        EEXCHAA=EEXCHAA+TMP3A
        EEXCHBB=EEXCHBB+TMP3B
C
 1010 CONTINUE
      IF(NUT.EQ.ILNBUF) GOTO 1001
      CLOSE(LUINT,STATUS='KEEP')
      ENDIF
C
      CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
      IF(YESNO)THEN
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
C
C     Two electron energy
C
      CALL LOCATE(LUINT,'TWOELSUP')
 2001 READ(LUINT)BUF,IBUF,NUT
      DO 2010 INT=1,NUT
        I=IUPKI(IBUF(INT))
        J=IUPKJ(IBUF(INT))
        K=IUPKK(IBUF(INT))
        L=IUPKL(IBUF(INT))
        X=BUF(INT)
C
C     alpha beta coulomb energy
           TMP1A=X*CDENSA((I-1)*N+J)*CDENSB((K-1)*N+L)
           IF(N*(I-1)+J.NE.N*(K-1)+L)
     &        TMP1A=TMP1A+X*CDENSA((K-1)*N+L)*CDENSB((I-1)*N+J)
C
C     alpha alpha coulomb energy
           TMP2A=X*CDENSA((I-1)*N+J)*CDENSA((K-1)*N+L)
           TMP2B=X*CDENSB((I-1)*N+J)*CDENSB((K-1)*N+L)
C
C     Off diagonal elements
           IF(N*(I-1)+J.NE.N*(K-1)+L)THEN
C     alpha alpha coulomb energy
              TMP2A=2.D+00*TMP2A
              TMP2B=2.D+00*TMP2B
           ENDIF
C
        ECOULAB=ECOULAB+TMP1A
        ECOULAA=ECOULAA+TMP2A
        ECOULBB=ECOULBB+TMP2B
C
 2010 CONTINUE
      IF(NUT.EQ.ILNBUF) GOTO 2001
      CLOSE(LUINT,STATUS='KEEP')
      ENDIF
C
 9999 CONTINUE
      ECOUL=ECOULAB+0.5D+00*(ECOULAA+ECOULBB)
      EEXCH=-0.5D+00*(EEXCHAA+EEXCHBB)
C
C     If finite field then nuclear dipole correction
      FNUC=0.D+00
      IF(IFLAGS(23).NE.0.OR.IFLAGS(24).NE.0.OR.IFLAGS(25).NE.0) THEN
         CALL GFNAME('VPOUT   ',FNAME,ILENGTH)
         OPEN(UNIT=30,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
         REWIND(30)
      CALL MAKMAP(N,NIRREP,NSOIRP,IMAP)
C
312      I23=IFLAGS(23)
         I24=IFLAGS(24)
         I25=IFLAGS(25)
         IF(IFLAGS(23).NE.0)THEN
            LABEL='     X  '
            IND=IFLAGS(23)
            IFLAGS(23)=0
         ELSEIF(IFLAGS(24).NE.0)THEN
            LABEL='     Y  '
            IND=IFLAGS(24)
            IFLAGS(24)=0
         ELSEIF(IFLAGS(25).NE.0)THEN
            LABEL='     Z  '
            IND=IFLAGS(25)
            IFLAGS(25)=0
         ELSE
            GOTO 23
         ENDIF
         FIELD=DFLOAT(IND)*1.E-6
C
         CALL SEEKLB(LABEL,IERR,0)
         BACKSPACE(30)
         READ(30)JUNK,DIPNUC
         FNUC=FNUC+FIELD*DIPNUC
C
123      READ(30)BUF,IBUF,NUT
         IF(NUT.EQ.-1)GOTO 312
            DO 313 INDEX=1,NUT
               IOFF=IMAP(IBUF(INDEX))
               IF(IOFF.NE.0)
     $            E1E=E1E+(DENSAA(IOFF)+DENSBB(IOFF))*BUF(INDEX)*FIELD
313         CONTINUE
         GOTO 123
C
23       CLOSE(30,STATUS='KEEP')
         IFLAGS(23)=I23
         IFLAGS(24)=I24
         IFLAGS(25)=I25
      ENDIF
C

      ESCF=REPULS+FNUC+E1E+ECOUL+EEXCH
      EHAR=ESCF-EEXCH
      WRITE(LUOUT,9020) E1E
      WRITE(LUOUT,9030) ECOUL
      WRITE(LUOUT,9040) EEXCH
      WRITE(LUOUT,9050) ESCF
 9020 FORMAT(/'The One electron energy =',F20.10)
 9030 FORMAT('The Coulomb Energy      =',F20.10)
 9040 FORMAT('The Exchange Enery      =',F20.10)
 9050 FORMAT('The SCF Energy          =',F20.10)
C
      RETURN
      END
