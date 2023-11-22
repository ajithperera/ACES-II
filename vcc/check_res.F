      Subroutine Check_Res(Work,Length,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null 
      Data onem,One /-1.0,1.0/
      Integer Dissyt 
#include "sympop.com"
#include "sym.com"
#include "syminf.com"
#include "parcc.com"


      INTEGER FSDPDAN,FSDPDNA,FSDPDAA,FSDPDIN,FSDPDNI,FSDPDII,FSDPDAI,
     $   FSDPDIA
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      Irrepx = 1

#ifdef DEBUG_LVLM
C T2-ABIJ
      Length_61=IDSYMSZ(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
      If (Iuhf .NE. 0) Length_62=IDSYMSZ(IRREPX,ISYTYP(1,62),
     +                                   ISYTYP(2,62))
      Length_63=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))

      Write(6,*) 
      Call Getall(Work, Length_61, Irrepx, 61)
      Call checksum("T2-AAAA :",Work,Length_61)
      If (Iuhf .ne. 0) Then
      Call Getall(Work, Length_62, Irrepx, 62)
      Call checksum("T2-BBBB :",Work,Length_62)
      Endif 
      Call Getall(Work, Length_63, Irrepx, 63)
      Call checksum("T2-ABAB :",Work,Length_63)

C T1-AI
      Length_aa =Irpdpd(Irrepx,9)
      call getlst(Work,1,1,1,1,90)
      Call checksum("T1-AA   :",Work,Length_aa)
      If (Iuhf .Ne. 0) Then 
      Length_bb =Irpdpd(Irrepx,10)
      call getlst(Work,1,1,1,2,90)
      Call checksum("T1-BB   :",Work,Length_bb)
      Endif 

      Write(6,*) 

      Write(6,"(a)") " The T1aa residual in list 90-3"
      call getlst(Work,1,1,1,3,90)

      Ioff = 1
      Do irrepr = 1, Nirrep
         irrepl = Dirprd(Irrepr, Irrepx)

         Nrow = vrt(Irrepl,1)
         Ncol = pop(Irrepr,1)

	 call output(work(ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
	 Ioff = Ioff + Nrow*Ncol
      Enddo

      Write(6,"(a)") " The T1bb residual in list 90-4"
      call getlst(Work,1,1,1,4,90)

      Ioff = 1
      Do irrepr = 1, Nirrep
         irrepl = Dirprd(Irrepr, Irrepx)

         Nrow = vrt(Irrepl,2)
         Ncol = pop(Irrepr,2)

	 call output(work(ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
	 Ioff = Ioff + Nrow*Ncol
      Enddo

      Do Ispin = 1, Iuhf+1
         Ioff = 1
         Do irrepr = 1, Nirrep
            irrepl = Dirprd(Irrepr, Irrepx)
     
             Nrow = Irpdpd(Irrepl,Ispin)
             Ncol = Irpdpd(Irrepr,Ispin+2)
             Call Getlst(Work(ioff),1,Ncol,1,Irrepr,60+Ispin)
             call output(work(ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
             Ioff = Ioff + Nrow*Ncol
         Enddo 
      Enddo 

      Ioff = 1
      Write(6,"(a)") " The T2 residual in list 63"
      Do irrepr = 1, Nirrep
         irrepl = Dirprd(Irrepr, Irrepx)

         Nrow = Irpdpd(Irrepl,13)
         Ncol = Irpdpd(Irrepr,14)

         Call Getlst(Work(ioff),1,Ncol,1,Irrepr,63)
         call output(work(ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
         Ioff = Ioff + Nrow*Ncol
      Enddo

#endif 
      
      DO 100 IRREP=1,NIRREP
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,63))
         IF(NUMSYT.EQ.0) GOTO 100
         DISSYT=FSDPDII(IRREP,ISYTYP(1,63))
         IF(DISSYT.EQ.0) GOTO 100
           write(*,'(A,I1,3X,A,I2,3X,A,I2,3X,A,I2,3X,A,I2)')
     $          "SPIN=",ISPIN,"IRREP=",IRREP,
     $          "NUMSYT=",NUMSYT,"DISSYT=",DISSYT
            CALL FSGET(Work,1,NUMSYT,1,IRREP,63,'IIAA')
            CALL OUTPUT(Work,1,DISSYT,1,NUMSYT,DISSYT,NUMSYT,1)
 100  CONTINUE

      Return
      End

