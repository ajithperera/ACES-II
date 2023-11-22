#include<aces.h>
#define _IIII
c THIS ROUTINE LOADS THE AO INTEGRALS FROM THE CORRESPONDING
c INTEGRAL FILE AND CONSTRUCTS THE FOCK MATRIX (ARRAY F) USING
c THE DENSITY MATRIX SUPPLIED IN D (UHF VERSION).
c
c COLOUMB TERMS:
c
c    FA(I,J) = FA(I,J) + DT(K,L) (IJ|KL)
c    FA(K,L) = FA(K,L) + DT(I,J) (IJ|KL)
c
c    FB(I,J) = FB(I,J) + DT(K,L) (IJ|KL)
c    FB(K,L) = FB(K,L) + DT(I,J) (IJ|KL)
c
c EXCHANGE TERMS:
c
c    FA(I,K) = FA(I,K) - 1/2 DA(J,L) (IJ|KL)
c    FA(J,L) = FA(J,L) - 1/2 DA(I,K) (IJ|KL)
c    FA(I,L) = FA(I,L) - 1/2 DA(J,K) (IJ|KL)
c    FA(J,K) = FA(J,K) - 1/2 DA(I,L) (IJ|KL)
c
c    FB(I,K) = FB(I,K) - 1/2 DB(J,L) (IJ|KL)
c    FB(J,L) = FB(J,L) - 1/2 DB(I,K) (IJ|KL)
c    FB(I,L) = FB(I,L) - 1/2 DB(J,K) (IJ|KL)
c    FB(J,K) = FB(J,K) - 1/2 DB(I,L) (IJ|KL)
c
      SUBROUTINE BUILT_FRAC_FOCK_AO(FA,FB,DA,DB,DT,BUF,IBUF,
     &                              IMAP,NBAST,NBAS,NTOTAL,
     &                              ILNBUF,LUINT) 

c   FA - SYMMETRY BLOCKED FOCK MATRIX (ALPHA)
c   FB - SYMMETRY BLOCKED FOCK MATRIX (BETA)
c   DA - SYMMETRY BLOCKED DENSITY MATRIX (ALPHA)
c   DB - SYMMETRY BLOCKED DENSITY MATRIX (BETA)
c
C  ARGUMENT LIST
      INTEGER NTOTAL, NBAST, NBAS(8), ILNBUF, LUINT
      DOUBLE PRECISION FA(NTOTAL), FB(NTOTAL)
      DOUBLE PRECISION DA(NTOTAL), DB(NTOTAL), DT(NTOTAL)
      DOUBLE PRECISION BUF(ILNBUF)
      INTEGER IBUF(ILNBUF), IMAP(NBAST,NBAST)
C
c---------------------------------------------------------------
c INTERNAL VARIABLES FOR/FROM VMOL

      DOUBLE PRECISION X
      INTEGER I,J,K,L,IOFF,IOFFT
      INTEGER IJ,IK,IL,JK,JL,KL,INDI,INDJ,INDK, INDL
      INTEGER IRREP,INDS,ILENGTH
      INTEGER NAOBUF,NUMINT,NUT
      CHARACTER*80 FNAME
      DOUBLE PRECISION TWOM,ONEM,HALFM,ZERO,FOURTH,HALF,ONE,TWO
      PARAMETER (TWOM=-2.0D0, ONEM=-1.0D0)
      PARAMETER (HALFM=-0.5D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (FOURTH=0.25D0, HALF=0.5D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)

c COMMON BLOCKS
#include "machsp.com"
#include "flags.com"

c FUNCTION DEFINITIONS
      INTEGER IAND,IOR,ISHFT
      INTEGER INT, IUPKI, IUPKJ, IUPKK, IUPKL, IPACK, INDX
      IUPKI(INT)=IAND(INT,IALONE)
      IUPKJ(INT)=IAND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=IAND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=IAND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=IOR(IOR(IOR(I,ISHFT(J,IBITWD)),
     &                             ISHFT(K,2*IBITWD)),
     &                             ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2

#ifdef _IIII
#ifdef _DEBUG_LVLM
        write(*,*) "The total density matrix"
        write(*,"(6(1x,f10.5))") (Dt(i),i=1,ntotal)
#endif 
      CALL GETREC(20,"JOBARC","NIRREP  ",1,NIRREP)
      IOFF=0
      IOFFT=0
      DO IRREP=1,NIRREP
         DO I=1,NBAS(IRREP)
            DO J=1,NBAS(IRREP)
               INDS=INDX(MAX(I,J),MIN(I,J))+IOFFT
                IMAP(I+IOFF,J+IOFF)=INDS
            END DO
          END DO
          IOFF=IOFF+NBAS(IRREP)
          IOFFT=IOFFT+NBAS(IRREP)*(NBAS(IRREP)+1)/2
      END DO
      CALL DZERO(FA,NTOTAL)
      CALL DZERO(FB,NTOTAL)

c IIII CONTRIBUTES TO BOTH COLOUMB AND EXCHANGE INTEGRALS
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      NAOBUF=0
      NUMINT=0
      CALL LOCATE(LUINT,'TWOELSUP')
 1    CONTINUE
         READ(LUINT) BUF, IBUF, NUT
         NAOBUF=NAOBUF+1
         DO INT=1,NUT
            X=BUF(INT)
            INDI=IUPKI(IBUF(INT))
            INDJ=IUPKJ(IBUF(INT))
            INDK=IUPKK(IBUF(INT))
            INDL=IUPKL(IBUF(INT))
            IJ=IMAP(INDJ,INDI)
            KL=IMAP(INDL,INDK)
            IK=IMAP(INDK,INDI)
            JL=IMAP(INDL,INDJ)
            IL=IMAP(INDL,INDI)
            JK=IMAP(INDJ,INDK)
            IF (INDI.EQ.INDJ) X=X*HALF
            IF (INDK.EQ.INDL) X=X*HALF
            IF (IJ.EQ.KL)     X=X*HALF
            FA(IJ)=FA(IJ)+DT(KL)*X
            FA(KL)=FA(KL)+DT(IJ)*X
            FB(IJ)=FB(IJ)+DT(KL)*X
            FB(KL)=FB(KL)+DT(IJ)*X
 
            FA(IK)=FA(IK)+HALFM*DA(JL)*X
            FA(JL)=FA(JL)+HALFM*DA(IK)*X
            FA(IL)=FA(IL)+HALfM*DA(JK)*X
            FA(JK)=FA(JK)+HALFM*DA(IL)*X
            FB(IK)=FB(IK)+HALFM*DB(JL)*X
            FB(JL)=FB(JL)+HALFM*DB(IK)*X
            FB(IL)=FB(IL)+HALfM*DB(JK)*X
            FB(JK)=FB(JK)+HALFM*DB(IL)*X

  50     continue

         END DO
         NUMINT=NUMINT+NUT
      IF (NUT.NE.-1) GOTO 1
      CLOSE(UNIT=LUINT,STATUS='KEEP')

#endif /* _IIII */

c SCALE DIAGONAL PARTS OF F BY A FACTOR OF TWO

      CALL SCALEF(FA,NTOTAL,NBAS)
      CALL SCALEF(FB,NTOTAL,NBAS)
      CALL XSCAL(NTOTAL,TWO,FA,1)
      CALL XSCAL(NTOTAL,TWO,FB,1)
#ifdef _DEBUG_LVLM
      write(*,*) "The Fock matrix in AO basis"
      write(*,"(6(1x,f10.5))") (Fa(i),i=1,ntotal)
#endif 
c
      RETURN
      END

