      SUBROUTINE BUILDA(CVECS,NBAS,VECIN,VECOUT,N,IAABEG,FOCK,LDIM1,
     &                  IUHF,LDIM3,INDR,NDET,LDIM2,SCR1,RAA,G1AA,G2AA,
     &                  SCRA,SCRB,SCRC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NOCC(16)
      COMMON /POPUL/ NOCC
#include "machsp.com"
#include "flags.com"
#include "symm2.com"

      PARAMETER(ONE=1.0D0,HALF=0.5d0,ZERO=0.0d0,ONEM=-1.0D0,
     &          TWO=2.0D0,FOUR=4.0D0)
      DATA SMALL/1.D-9/

      DIMENSION VECIN(N),VECOUT(N),SCR1(N)
      DIMENSION CVECS(NBAS,NBAS)
      DIMENSION FOCK((IUHF+1)*LDIM1)
      DIMENSION IAABEG(16),NDET(16)
      DIMENSION RAA(LDIM3*(IUHF+1)),G1AA(LDIM3*(IUHF+1)),
     $          G2AA(LDIM3*(IUHF+1))
      DIMENSION BUF(600),IBUF(600)
      DIMENSION INDR(NBAS,NBAS)
      DIMENSION SCRA(LDIM2),SCRB(LDIM2),SCRC(LDIM2)
      CHARACTER*80 FNAME
      CHARACTER*8 SCFREC(2)
      DATA SCFREC/'SCFEVCA0','SCFEVCB0'/
      LOGICAL BOPEN
C***************
C     PRELIMINARIES NECESSARY TO READ THE INTEGRAL FILES
C***************
      IUPKI(INT)=IAND(INT,IALONE)
      IUPKJ(INT)=IAND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=IAND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=IAND(ISHFT(INT,-3*IBITWD),IALONE)

      NDX=LDIM3*(IUHF+1)
      DO I=1,NDX
         RAA(I)=ZERO
         G1AA(I)=ZERO
         G2AA(I)=ZERO
      END DO

c      DO I=1,N
c      SCR1(I)=ZERO
c      END DO

c      ICOUNT=1
c      DO I=1,NIRREP
c         ISTART=IREPS(I)
c         IEND=IREPS(I)+NBFIRR(I)-1
c         DO J=ISTART,IEND
c         DO K=ISTART,IEND
c         INDR(K,J)=ICOUNT
c         ICOUNT=ICOUNT+1
c         END DO
c         END DO
c      END DO

C BUILD THE R INTERMEDIATES FROM PAPER GENERALIZED TO UHF MATRICES
C    R(ILA,ISIG)=SUM(JB) C(ILA,J)*C(ISIG,B)*VECIN(JB)
C WHERE C IS MO COEFFICIENT AND VECIN IS INPUT VECTOR
C AND J,B RUN OVER OCCUPIED/VIRTUAL SPACES RESPECTIVELY
C R IS BLOCKED BY SYMMETRY AND NEED ONE FOR ALPHA AND BETA
C IN GENERAL. I WILL WRITE THIS AS A SERIES OF MATRIX MULTIPLICATIONS (I HOPE)

      NFOCK=IUHF+1
      DO ISPIN=1,NFOCK
         CALL GETREC(1,'JOBARC',SCFREC(ISPIN),NBAS*NBAS*IINTFP,CVECS)
         DO I=1,NIRREP
         IF( NDET( I + (ISPIN-1)*8)  .NE.  0)THEN
           CALL GETBLK(CVECS,SCRA,NBFIRR(I),NBAS,IREPS(I))
           NOCCA=NOCC(I + (ISPIN-1)*8 )
           NVIRTA=NBFIRR(I)-NOCCA
           NROWA=NOCCA*NBFIRR(I)+1
           NROWB=IAABEG(I + (ISPIN-1)*8 )
           CALL MXM(SCRA(NROWA),NBFIRR(I),VECIN(NROWB),
     $              NVIRTA,SCRC,NOCCA)
           CALL TRANSP(SCRC,SCRB,NOCCA,NBFIRR(I))
           CALL MXM(SCRA(1),NBFIRR(I),SCRB,
     $              NOCCA,RAA(((ISPIN-1)*LDIM3)+ISQROF(I)),NBFIRR(I))
         END IF
         END DO
      END DO

C   COMPUTE G1 AND G2 INTERMEDIATES LISTS
C   G1(I,J)=SUM(K,L) RAA(K,L)*(I,J,K,L) =>COULOMB LIST
C   G2(I,J)=SUM(K,L) RAA(K,L)*(I,K,J,L) =>EXCHANGE LIST
C   G1 IS SYMMETRIC HOWEVER G2 IS NOT SINCE RAA(K,L)/=RAA(L,K)

C READ INTEGRALS AND PLACE THEM IN EVERY SUM
C TO WHICH THEY CONTRIBUTE.
C IIII INTEGRALS CONTRIBUTE TO G1 AND G2 INTERMEDIATES
      ISTORE=10
      INQUIRE(UNIT=ISTORE,OPENED=BOPEN) ! DO NOT OVERWRITE BOPEN
      IF (.NOT.BOPEN) THEN
         OPEN(UNIT=ISTORE,FORM='UNFORMATTED',FILE='IIII',
     &        ACCESS='SEQUENTIAL',STATUS='OLD')
      END IF
      CALL LOCATE(ISTORE,'TWOELSUP')
 4    READ(ISTORE) BUF, IBUF, NUT
C     SELECT INTEGRAL (I,J|K,L)
C***********
      DO IND=1,NUT
         IMU=IUPKI(IBUF(IND))
         INU=IUPKJ(IBUF(IND))
         ILA=IUPKK(IBUF(IND))
         ISIG=IUPKL(IBUF(IND))
         X=BUF(IND)
         SCALE=1.0D0
         IF(IMU.EQ.INU)SCALE=HALF
         IF(ILA.EQ.ISIG)SCALE=SCALE*HALF
         if(IMU.EQ.INU.and.ILA.EQ.ISIG.and.IMU.EQ.ILA)SCALE=SCALE*HALF
         if(IMU.EQ.ILA.and.INU.EQ.ISIG.and.IMU.NE.INU)SCALE=HALF
         X=X*SCALE

C PLOP THE (IMU,INU,ILA,ISIG) INTEGRAL INTO ALL MATRIX
C ELEMENTS TO WHICH IT MAKES A CONTRIBUTION
C SORRY THIS IS LOWER CASE...
C REMOVE REDUNDANT ADDITIONS AND MATRIX ELEMENT
C LOOKUPS AFTER DEBUGGING!!!!!!!!!!!!!!!!!!!!!

         i1=indr(imu,inu)
         i2=indr(ila,isig)
         i3=indr(isig,ila)
         i4=indr(inu,imu)
         t1=(raa(i2)+raa(i3))*x
         t2=(raa(i1)+raa(i4))*x
         g1aa(i1)=g1aa(i1)+t1
         g1aa(i4)=g1aa(i4)+t1
         g1aa(i2)=g1aa(i2)+t2
         g1aa(i3)=g1aa(i3)+t2

         i5=indr(imu,ila)
         i6=indr(inu,isig)
         i7=indr(imu,isig)
         i8=indr(inu,ila)
         i9=indr(ila,imu)
         i10=indr(isig,inu)
         i11=indr(ila,inu)
         i12=indr(isig,imu)

         g2aa(i5)=g2aa(i5)+raa(i6)*x
         g2aa(i7)=g2aa(i7)+raa(i8)*x
         g2aa(i8)=g2aa(i8)+raa(i7)*x
         g2aa(i6)=g2aa(i6)+raa(i5)*x
         g2aa(i9)=g2aa(i9)+raa(i10)*x
         g2aa(i11)=g2aa(i11)+raa(i12)*x
         g2aa(i12)=g2aa(i12)+raa(i11)*x
         g2aa(i10)=g2aa(i10)+raa(i9)*x

C IF UHF MAKE BETA INTERMEDIATES
         IF(IUHF.NE.0)THEN
         i1=i1+ldim3
         i2=i2+ldim3
         i3=i3+ldim3
         i4=i4+ldim3
         t1=(raa(i2)+raa(i3))*x
         t2=(raa(i1)+raa(i4))*x
         g1aa(i1)=g1aa(i1)+t1
         g1aa(i4)=g1aa(i4)+t1
         g1aa(i2)=g1aa(i2)+t2
         g1aa(i3)=g1aa(i3)+t2

         i5=i5+ldim3
         i6=i6+ldim3
         i7=i7+ldim3
         i8=i8+ldim3
         i9=i9+ldim3
         i10=i10+ldim3
         i11=i11+ldim3
         i12=i12+ldim3

         g2aa(i5)=g2aa(i5)+raa(i6)*x
         g2aa(i7)=g2aa(i7)+raa(i8)*x
         g2aa(i8)=g2aa(i8)+raa(i7)*x
         g2aa(i6)=g2aa(i6)+raa(i5)*x
         g2aa(i9)=g2aa(i9)+raa(i10)*x
         g2aa(i11)=g2aa(i11)+raa(i12)*x
         g2aa(i12)=g2aa(i12)+raa(i11)*x
         g2aa(i10)=g2aa(i10)+raa(i9)*x

         END IF

      END DO
      IF(NUT.EQ.600)GOTO 4
      CLOSE(UNIT=ISTORE,STATUS='KEEP')

c FOR NIRREP > 1, READ ALSO IIJJ AND IJIJ INTEGRAL FILES

      IF (NIRREP.GT.1) THEN
c IIJJ CONTRIBUTES ONLY TO G1 COLOUMB INTEGRALS
c       OPEN(UNIT=ISTORE,FORM='UNFORMATTED',FILE='IIJJ',
c     &ACCESS='SEQUENTIAL',STATUS='OLD')

      CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
      OPEN(UNIT=ISTORE,FILE=FNAME(1:ILENGTH),
     &     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')


      CALL LOCATE(ISTORE,'TWOELSUP')

 5    CONTINUE
      READ(ISTORE) BUF, IBUF, NUT
      DO IND=1,NUT
         IMU=IUPKI(IBUF(IND))
         INU=IUPKJ(IBUF(IND))
         ILA=IUPKK(IBUF(IND))
         ISIG=IUPKL(IBUF(IND))
         X=BUF(IND)
         SCALE=1.0D0
         IF(IMU.EQ.INU)SCALE=HALF
         IF(ILA.EQ.ISIG)SCALE=SCALE*HALF
         X=X*SCALE

C PLOP THE (IMU,INU,ILA,ISIG) INTEGRAL INTO ALL MATRIX
C ELEMENTS TO WHICH IT MAKES A CONTRIBUTION
C SORRY THIS IS LOWER CASE...
C REMOVE REDUNDANT ADDITIONS AND MATRIX ELEMENT
C LOOKUPS AFTER DEBUGGING!!!!!!!!!!!!!!!!!!!!!

         i1=indr(imu,inu)
         i2=indr(ila,isig)
         i3=indr(isig,ila)
         i4=indr(inu,imu)
         t1=(raa(i2)+raa(i3))*x
         t2=(raa(i1)+raa(i4))*x
         g1aa(i1)=g1aa(i1)+t1
         g1aa(i4)=g1aa(i4)+t1
         g1aa(i2)=g1aa(i2)+t2
         g1aa(i3)=g1aa(i3)+t2

C IF UHF MAKE BETA INTERMEDIATES
         IF(IUHF.NE.0)THEN
         i1=i1+ldim3
         i2=i2+ldim3
         i3=i3+ldim3
         i4=i4+ldim3
         t1=(raa(i2)+raa(i3))*x
         t2=(raa(i1)+raa(i4))*x
         g1aa(i1)=g1aa(i1)+t1
         g1aa(i4)=g1aa(i4)+t1
         g1aa(i2)=g1aa(i2)+t2
         g1aa(i3)=g1aa(i3)+t2
         END IF

      end do
      IF(NUT.EQ.600)GOTO 5
      CLOSE(UNIT=ISTORE,STATUS='KEEP')

c IJIJ CONTRIBUTES ONLY TO EXCHANGE
      OPEN(UNIT=ISTORE,FORM='UNFORMATTED',FILE='IJIJ',
     &ACCESS='SEQUENTIAL',STATUS='OLD')
      CALL LOCATE(ISTORE,'TWOELSUP')

 6    READ(ISTORE) BUF, IBUF, NUT
C     SELECT INTEGRAL (I,J|K,L)
C***********
      DO IND=1,NUT
         IMU=IUPKI(IBUF(IND))
         INU=IUPKJ(IBUF(IND))
         ILA=IUPKK(IBUF(IND))
         ISIG=IUPKL(IBUF(IND))
         X=BUF(IND)
         SCALE=1.0D0
         if(IMU.EQ.ILA .and. INU.EQ.ISIG)SCALE=HALF
         X=X*SCALE

C PLOP THE (IMU,INU,ILA,ISIG) INTEGRAL INTO ALL MATRIX
C ELEMENTS TO WHICH IT MAKES A CONTRIBUTION
C SORRY THIS IS LOWER CASE...
C REMOVE REDUNDANT ADDITIONS AND MATRIX ELEMENT
C LOOKUPS AFTER DEBUGGING!!!!!!!!!!!!!!!!!!!!!

         i5=indr(imu,ila)
         i6=indr(inu,isig)
         i9=indr(ila,imu)
         i10=indr(isig,inu)

         g2aa(i5)=g2aa(i5)+raa(i6)*x
         g2aa(i6)=g2aa(i6)+raa(i5)*x
         g2aa(i9)=g2aa(i9)+raa(i10)*x
         g2aa(i10)=g2aa(i10)+raa(i9)*x

C IF UHF MAKE BETA INTERMEDIATES
         IF(IUHF.NE.0)THEN
         i5=i5+ldim3
         i6=i6+ldim3
         i9=i9+ldim3
         i10=i10+ldim3

         g2aa(i5)=g2aa(i5)+raa(i6)*x
         g2aa(i6)=g2aa(i6)+raa(i5)*x
         g2aa(i9)=g2aa(i9)+raa(i10)*x
         g2aa(i10)=g2aa(i10)+raa(i9)*x
         END IF

      end do
      IF(NUT.EQ.600)GOTO 6
      CLOSE(UNIT=ISTORE,STATUS='KEEP')

c     END IF (NIRREP.GT.1)
      END IF

      DO I=1,N
         VECOUT(I)=ZERO
      END DO

C HERE WE TRANSFORM THE INTERMEDIATES TO THE MO BASIS
C AND BUILD THE VECOUT(IA) ELEMENTS

      NFOCK=IUHF+1
      DO ISPIN=1,NFOCK
         CALL GETREC(1,'JOBARC',SCFREC(ISPIN),NBAS*NBAS*IINTFP,CVECS)
         DO I=1,NIRREP
         IF( NDET( I + (ISPIN-1)*8)  .NE.  0)THEN
           CALL GETBLK(CVECS,SCRA,NBFIRR(I),NBAS,IREPS(I))
           NOCCA=NOCC(I + (ISPIN-1)*8 )
           NVIRTA=NBFIRR(I)-NOCCA
           NROWA=NOCCA*NBFIRR(I)+1
           NROWB=IAABEG(I + (ISPIN-1)*8 )
           CALL MXM(G1AA(((ISPIN-1)*LDIM3)+ISQROF(I)),NBFIRR(I),
     $              SCRA(1),NBFIRR(I),SCRC,NOCCA)
           CALL TRANSP(SCRA(NROWA),SCRB,NVIRTA,NBFIRR(I))
           CALL MXM(SCRB,NVIRTA,SCRC,NBFIRR(I),SCR1(NROWB),NOCCA)

           IF(IUHF.EQ.0)THEN
             CALL DAXPY(NDET(I),FOUR,SCR1(NROWB),1,VECOUT(NROWB),1)
           ELSE
             CALL DAXPY(NDET(I+(ISPIN-1)*8),
     &                  TWO,SCR1(NROWB),1,VECOUT(NROWB),1)
             IF(ISPIN.EQ.1)THEN
               CALL MXM(G1AA(LDIM3+ISQROF(I)),NBFIRR(I),
     $                  SCRA(1),NBFIRR(I),SCRC,NOCCA)
               CALL MXM(SCRB,NVIRTA,SCRC,NBFIRR(I),SCR1(NROWB),NOCCA)
               CALL DAXPY(NDET(I),TWO,SCR1(NROWB),1,VECOUT(NROWB),1)
             ELSE
               CALL MXM(G1AA(ISQROF(I)),NBFIRR(I),
     $                  SCRA(1),NBFIRR(I),SCRC,NOCCA)
               CALL MXM(SCRB,NVIRTA,SCRC,NBFIRR(I),SCR1(NROWB),NOCCA)
               CALL DAXPY(NDET(I+8),TWO,SCR1(NROWB),1,VECOUT(NROWB),1)
             END IF
           END IF

           CALL MXM(G2AA(((ISPIN-1)*LDIM3)+ISQROF(I)),NBFIRR(I),
     $              SCRA(1),NBFIRR(I),SCRC,NOCCA)
           CALL MXM(SCRB,NVIRTA,SCRC,NBFIRR(I),SCR1(NROWB),NOCCA)
           CALL DAXPY(NDET(I + (ISPIN-1)*8),ONEM,SCR1(NROWB),1,
     $                VECOUT(NROWB),1)

           CALL TRANSP(G2AA(((ISPIN-1)*LDIM3)+ISQROF(I)),SCRC,
     $                 NBFIRR(I),NBFIRR(I))
C USE RAA AS A SCRATCH ARRAY HERE SINCE WE DON'T NEED IT ANYMORE
           CALL MXM(SCRC,NBFIRR(I),SCRA(1),NBFIRR(I),RAA,NOCCA)
           CALL MXM(SCRB,NVIRTA,RAA,NBFIRR(I),SCR1(NROWB),NOCCA)
           CALL DAXPY(NDET(I + (ISPIN-1)*8),ONEM,SCR1(NROWB),1,
     $                VECOUT(NROWB),1)

C ADD IN ONE ELECTRON CONTRIBUTIONS

           CALL EXPND2(FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     $                 SCRA,NBFIRR(I))
           ISTART=NOCCA+1
           IEND=NBFIRR(I)
           CALL QCGETF(SCRA,NBFIRR(I),SCRB,NVIRTA,ISTART,IEND)
           CALL MXM(SCRB,NVIRTA,VECIN(NROWB),NVIRTA,SCR1(NROWB),NOCCA)
           CALL DAXPY(NDET(I + (ISPIN-1)*8),ONE,SCR1(NROWB),1,
     $                VECOUT(NROWB),1)

           ISTART=1
           IEND=NOCCA
           CALL QCGETF(SCRA,NBFIRR(I),SCRB,NOCCA,ISTART,IEND)
           CALL MXM(VECIN(NROWB),NVIRTA,SCRB,NOCCA,SCR1(NROWB),NOCCA)
           CALL DAXPY(NDET(I + (ISPIN-1)*8),ONEM,SCR1(NROWB),1,
     $                VECOUT(NROWB),1)

         END IF
         END DO
       END DO

C IF RHF COPY THE UPPER HALF (ALPHA) INTO THE LOWER HALF (BETA)
C SINCE WE ARE NOT USING SPIN ADAPTED CONFIGURATIONS

       IF (IUHF.EQ.0) THEN
          NDX=0
          DO I=1,NIRREP
             NDX=NDX+NDET(I)
          END DO
          CALL DCOPY(NDX,VECOUT(1),1,VECOUT(NDX+1),1)
       END IF

c       DO I=1,N
c       IF(DABS(VECOUT(I)).LT.SMALL)VECOUT(I)=ZERO
c       END DO

c Re-open IIII if it was open on entry (as with FOCK=PK).
      IF (BOPEN) THEN
         OPEN(UNIT=ISTORE,FORM='UNFORMATTED',FILE='IIII',
     &        ACCESS='SEQUENTIAL',STATUS='OLD')
      END IF

       RETURN
       END

