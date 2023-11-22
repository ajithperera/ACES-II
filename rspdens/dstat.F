      SUBROUTINE DSTAT(String,V,LENGTH1,LENGTH2,SYMMETRIZE,TYPE,ISPIN)
C
C RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
C     ZQ(1)  Largest absolute magnitude
C     ZQ(2)  Smallest absolute magnitude
C     ZQ(3)  Sum 
C   

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LENGTH1*LENGTH2),ZQ(4)
      DIMENSION U(LENGTH1*LENGTH2)
      LOGICAL SYMMETRIZE 
      INTEGER POP,VRT,DIRPRD 
      INTEGER FIN1,FIN2
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1(2),NF2(2)

      CHARACTER*7 STRING
      CHARACTER*2 TYPE

      D_OLD_MAX = 0.0D0
      D_OLD_MIN = 1.0D0
      ZQ(3)     = 0.0D0
      ZQ(4)     = 0.0D0
      Ioff      = 1
      Joff      = 1
      IF (SYMMETRIZE) THEN
         IF (TYPE .EQ. "OO") CALL SYMMET5(V,U,POP(1,ISPIN))
         IF (TYPE .EQ. "VV") CALL SYMMET5(V,U,VRT(1,ISPIN))
      ENDIF 

      DO IRREP=1,NIRREP
         IF (TYPE .EQ. "OO") THEN
            FIN1 = POP(IRREP,ISPIN)
            FIN2 = POP(IRREP,ISPIN)
         ELSEIF (TYPE .EQ. "VV") THEN
            FIN1 = VRT(IRREP,ISPIN)
            FIN2 = VRT(IRREP,ISPIN)
         ELSEIF (TYPE .EQ. "VO") THEN
            FIN1 = POP(IRREP,ISPIN)
            FIN2 = VRT(IRREP,ISPIN)
         ENDIF 

         Call output(V(ioff),1,Fin1,1,Fin2,Fin1,Fin2,1)

         Call GETSTAT(v(Ioff),FIN1,FIN2,D_NEW_MAX,
     &                D_NEW_MIN,D_OLD_MAX,D_OLD_MIN,
     &                IMAX,JMAX,IMIN,JMIN,IRREP,IRREP_MIN,
     &                IRREP_MAX,ZQ)
         Ioff = Ioff + Fin1*Fin2
      ENDDO

      ZQ(1) = D_OLD_MAX 
      IF (D_OLD_MIN .EQ. 1.0D0) THEN
         ZQ(2) = 0.0D0
      ELSE
         ZQ(2) = D_OLD_MIN
      ENDIF 

      WRITE(6,"(a)") String 
      WRITE(6,"(a,E24.12,3I4)") "Maximum,indices and irrep of maximun:",
     &                           ZQ(1),IMAX,JMAX,IRREP_MAX
      WRITE(6,"(a,E24.12,3I4)") "Minumum,indices and irrep of minimum:",
     &                           ZQ(2),IMIN,JMIN,IRREP_MIN
      WRITE(6,"(a,E24.12)")     "Sum D(p,q)                          :",
     &                           ZQ(3)
      WRITE(6,"(a,E24.12)")     "Sum D(p,q)^2                        :",
     &                           ZQ(4)
      WRITE(6,"(a)")            "p,q"
      WRITE(6,*)

      RETURN
      END
