 
      SUBROUTINE PUTSTF(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      CALL PUTREC(20,'JOBARC','NOCCORB ',2,NOCCO)
      CALL PUTREC(20,'JOBARC','NVRTORB ',2,NVRTO)
      CALL PUTREC(20,'JOBARC','UHFRHF  ',1,IUHF)
      NORBA=NOCCO(1)+NVRTO(1)
      NORBB=NOCCO(2)+NVRTO(2)
      NSIZA=(NORBA*(NORBA+1))/2
      NSIZB=(NORBB*(NORBB+1))/2
      IR1=1
      IR2=IR1+NORBA*NORBA*IINTFP
      IR3=IR2+NORBA*NORBA*IINTFP
C Noticed the FOCKLIST arguments have changed. I am not doing
C anything here since focklist call is commented here.
C Ajith Perera, 05/2014.

C      CALL FOCKLIST(ICORE(IR1),ICORE(IR2),NORBA,IUHF)
      RETURN
      END
