      Program libc4
 
      Implicit double precision(A-H,O-Z)
      Logical T2_fine

#include "icore.com"
#include "istart.com"
#include "machsp.com"

      Data Ione, Dzero /1, 0.0D0/
     
      Call Crapsi_c4(Icore,I0,Icrsiz,Iuhf,0)

      Write(6,"(10x,a)") " -----------Interface to CFOUR------"
      Write(6,*) 

      Call Getrec_c4(20,"JOBARC","NATOMS ", Ione, Natoms)
      Call Getrec_c4(20,"JOBACR","COORD  ", 3*Natoms,Icore(I0))

      Write(6,"(a,1x,I2)") " Hello CFOUR world!"
 
      Write(6,"(a,1x,I3)") " The number of  atoms: ", Natoms 
      Write(6,"(a)")       " Molecular geomery   :" 
      Call output(Icore(I0),1,Natoms,1,3,Natoms,3,1)
      
      Distance=Ddot(3*Natoms,Icore(I0),1,Icore(I0),1)
      If (Natoms .Gt. 0 .And. Distance .Gt. Dzero) Then 
         Write(6,*)
         Write(6,"(a)") " Communication with JOBARC is established!"
      Endif 
     
      Call Getrec_C4(20,"JOBARC",'SCFENEG ',Ione*IINTFP,ESCF) 
      Call Getrec_C4(20,"JOBARC",'TOTENERG',Ione*IINTFP,ETOT) 

      Ediff = Etot - Escf 
   
      If (Ediff .Ne. Dzero) Then
         Call Check_t2(Icore(I0),Icrsiz,Iuhf,T2_Fine)
      Endif 

      If (T2_fine) Then
         Write(6,*)
         Write(6,"(a)") " Communication with MOIO is established!"
      Endif 

      Write(6,*)
      Call Crapso_c4()

CSSS      Call Aces_io_fin_c4
CSSS      Call aces_ja_fin_c4

      Stop
      End
