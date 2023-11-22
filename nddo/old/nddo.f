      subroutine NDDO(NBAS,NATOMS,IRATIO,NBFIRR,NIRREP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


C         **************************************
C         *        NDDO initial guess          *
C         **************************************
C
C         Written by Carlos Taylor
C
C 
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD


      DIMENSION OVLAP(NBAS*NBAS),S2(nbas,nbas)
      DIMENSION XYZ(3,NATOMS),NUCCHG(NATOMS),ISHL(20)
     $,BCOORD(3)

C     SET A MAX ARRAY SIZE ON THE NUMBER OF CONTRACTION COEFFICIENTS
C     AND NUMBER OF ATOMIC ORBITAL BASIS FUNCTIONS FOR THAT SHELL
C     I WOULD HOPE THIS WILL HANDLE MOST CASES
      PARAMETER(NCMAX=20,NBMAX=30)
      DIMENSION CTEMP(NCMAX,NBMAX),ETEMP(NCMAX)
C     THE FOLLOWING ARRAY WILL HOLD EACH BASIS FN'S CONTRACTION COEF
C     AND EXPONENTS IN THE CORRESPONDING ROW
      DIMENSION CONTR(NCMAX,NBAS),EXPON(NCMAX,NBAS),COORB(3,NBAS)
C     NGAUSS IS THE NUMBER OR PRIMITIVES FOR EACH BASIS FUNCTION WHICH
C     IS REQUIRED FOR THE OVERLAP INTEGRALS BETWEEN STO AND GTO 
      DIMENSION NGAUSS(NBAS)
C     LONE WILL HOLD THE TYPE OF FUNCTION (S=1,X=2,Y=3,Z=4
C     XX=5,YY=6,ZZ=7,xy=8,xz=9,yz=10) 
      DIMENSION LONE(NBAS),IANGX(6),ISALC(8),NBFIRR(8)
      CHARACTER*80 XLINE
      CHARACTER*8 ABC,LABEL(10)
      DATA LABEL/'S','X','Y','Z','XX','YY','ZZ','XY','XZ','YZ'/
      DATA IANGX/5,8,9,6,10,7/
      PARAMETER(ZERO=0.0D0)

      CHARACTER*80 FNAME
      LOGICAL YESNO2
      CHARACTER*8 TITLE(24)

C
C     GET THE GTO BASIS FUNCTION INFORMATION SO THAT I CAN
C     BUILD THE GTO-NDDO OVERLAP MATRIX 
C     THE BASIS INFO IS STORED ON THE 'MOL' FILE
C     READ THE INFO AS INSTRUCTED BY AJITH'S SUBROUTINES
C

      OPEN(UNIT=10, FILE='MOL_NDDO')
      REWIND(10)
C READ FIRST FIVE LINES OF MOL FILE WHICH WE DON'T NEED
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
      READ(10,'(A)') XLINE
C     ICOUNT IS A COUNTER FOR THE TOTAL NUMBER OF BASIS FUNCTIONS
      ICOUNT=0
C LOOP OVER ATOMS AND READ BASIS FUNCTIONS FOR EACH CENTER
C (SEE MOL FILE FOR SYNTAX)
      CALL GETREC(-1,'JOBARC','ATOMCHRG',NATOMS,NUCCHG)
      DO IATOM=1,NATOMS
         if(nucchg(iatom).eq.0)goto 321
  

         READ(10, 1110) ZNUC, IJUNK, NSHL,(ISHL(I),I=1,NSHL)
 1110    FORMAT(F20.1,8I5)
C
         READ(10,1115) ABC,(BCOORD(I), I = 1, 3)
 1115    FORMAT(A4,3F20.12)
C            LOOP OVER SHELLS ON EACH ATOM
              DO ISHELL=1,NSHL
C            LOOP OVER SECTIONS FOR EACH SHELL (CASES WHERE EACH
C            SHELL IS SPLIT OVER MORE THAN ONE SECTION)
              DO JSHELL=1,ISHL(ISHELL)
                   READ(10,1120) NP1, NAO
 1120              FORMAT(2I5)
C                       LOOP OVER PRIMITIVES IN EACH SHELL
                        DO J=1,NP1
                        READ(10,*)ETEMP(J),(CTEMP(J,K),K=1,NAO)
C                        WRITE(*,*)ETEMP(J),(CTEMP(J,K),K=1,NAO)
                        END DO

C     PUT THE CONTRACTION COEFS AND EXPONENTS IN BIG ARRAY 
                    IF(ISHELL.EQ.1)THEN
C              THESE ARE S ORBITALS SO ONLY ONE OF EACH
                          DO IAO=1,NAO
                            ICOUNT=ICOUNT+1
                            IROW=0
                             DO JPR=1,NP1
                              IF(CTEMP(JPR,IAO).NE.ZERO)THEN
                              IROW=IROW+1
                              CONTR(IROW,ICOUNT)=CTEMP(JPR,IAO)
                              EXPON(IROW,ICOUNT)=ETEMP(JPR)
                              END IF
                             END DO
                            NGAUSS(ICOUNT)=IROW
                            LONE(ICOUNT)=1
                            COORB(1,ICOUNT)=BCOORD(1)
                            COORB(2,ICOUNT)=BCOORD(2)
                            COORB(3,ICOUNT)=BCOORD(3)
                          END DO
                    ELSEIF(ISHELL.EQ.2)THEN
C           ISHELL=2 MEANS P ORBITALS SO NEED 3 OF EACH P ORBITAL
C  APPARENTLY ACES STORES THE ORBITALS AS X,X',Y,Y',Z,Z'
C  AS OPPOSED TO X,Y,Z,X',Y',Z' (WHICH IS THE ORDERING THAT IS 
C  INDICATED IN THE EIGENVECTOR PRINTOUT).
C  ANYWAY, I WILL REWRITE THIS SECTION.
                  DO L=1,3 
C     FOR GOD'S SAKE, DON'T EVER CHANGE THE ORDER OF THE FOLLOWING 
C     TWO LINES OR THE COUNTER WILL BE ALL MESSED UP
                     ICOLUMN=ICOUNT+1
                     ICOUNT=ICOUNT+NAO
                          DO IAO=1,NAO
                               ISPOT=ICOLUMN+IAO-1
                               IROW=0
                                  DO JPR=1,NP1
                                     IF(CTEMP(JPR,IAO).NE.ZERO)THEN
                                     IROW=IROW+1
                                     CONTR(IROW,ISPOT)=CTEMP(JPR,IAO)
                                     EXPON(IROW,ISPOT)=ETEMP(JPR)
                                     END IF
                                   END DO
                             NGAUSS(ISPOT)=IROW
                             LONE(ISPOT)=L+1
                             COORB(1,ISPOT)=BCOORD(1)
                             COORB(2,ISPOT)=BCOORD(2)
                             COORB(3,ISPOT)=BCOORD(3)
                          END DO
                    END DO

                ELSEIF(ISHELL.EQ.3)THEN
C           ISHELL=3 MEANS d ORBITALS SO NEED 6 OF EACH D ORBITAL
C           THE CARTESIAN TO SPHERICAL TRANSFORM WILL HAVE TO BE
C           DONE BY ACES

                    DO L=1,6
C     FOR GOD'S SAKE, DON'T EVER CHANGE THE ORDER OF THE FOLLOWING 
C     TWO LINES OR THE COUNTER WILL BE ALL MESSED UP
                     ICOLUMN=ICOUNT+1
                     ICOUNT=ICOUNT+NAO
                          DO IAO=1,NAO
                               ISPOT=ICOLUMN+IAO-1
                               IROW=0
                                  DO JPR=1,NP1
                                     IF(CTEMP(JPR,IAO).NE.ZERO)THEN
                                     IROW=IROW+1
                                     CONTR(IROW,ISPOT)=CTEMP(JPR,IAO)
                                     EXPON(IROW,ISPOT)=ETEMP(JPR)
                                     END IF
                                   END DO
                             NGAUSS(ISPOT)=IROW
                             LONE(ISPOT)=IANGX(L)
                             COORB(1,ISPOT)=BCOORD(1)
                             COORB(2,ISPOT)=BCOORD(2)
                             COORB(3,ISPOT)=BCOORD(3)
                          END DO
                    END DO
                END IF
C END LOOP OVER SECTIONS OF EACH SHELL
             END DO
C END LOOP OVER SHELLS
             END DO




C END LOOP OVER ATOMS
 321         continue
      END DO


C     PRINT THE FOLLOWING TO DEBUG:
C     LIST OF COORDINATES, TYPE(S,X,Y,Z,ETC),CONTRACTION COEFS
C     AND EXPONENTS FOR EACH BASIS FUNCTION IN THE COMPUTATIONAL ORDER


c       WRITE(*,*)'TOTAL NUMBER OF CONTRACTED GAUSSIANS',ICOUNT
c      do i=1,ICOUNT
c      write(*,*)'basis function number',i
c      write(*,*)'coordinates',coorb(1,i),coorb(2,i),coorb(3,i)
c      write(*,*)'type  ',label(lone(i))
c      write(*,*)'------------------------'
c      do j=1,ngauss(i)
c      write(*,*)expon(j,i),contr(j,i)
c      end do
c      end do

C     PUT THE BASIS INFORMATION ON DISK 
C
      OPEN(UNIT=50,FILE='GTO')
      WRITE(50,*)NBAS
      WRITE(50,*)NCMAX
      WRITE(50,*)LONE
      WRITE(50,*)NGAUSS
      WRITE(50,*)CONTR
      WRITE(50,*)EXPON
      WRITE(50,*)COORB
      CLOSE(50)
      
      IRATIO=IINTFP


         LUINT=10
         CALL GFNAME('IIII    ',FNAME,ILENGTH)
         INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO2)
         IF (YESNO2) THEN
            OPEN(LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &           ACCESS='SEQUENTIAL')

         READ(LUINT)TITLE,NIRREP,(NBFIRR(I),I=1,NIRREP),REPULS,
     &                 IIJUNK
            CLOSE(LUINT,STATUS='KEEP')
         END IF

         



      return
      end 





