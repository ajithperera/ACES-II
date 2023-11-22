      Subroutine Reord_2zmat(Atmmass,Atmlabel,Iatmchrg,Imap,
     +                       Atmmass_zmtord,Atmlabel_zmtord,
     +                       Iatmchrg_zmtord,Natoms,Nreals)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Imap(Natoms)
      Dimension Atmmass(Natoms)
      Dimension Atmmass_zmtord(Nreals)
      Dimension Iatmchrg(Natoms)
      Dimension Iatmchrg_zmtord(Nreals)
      Character*5 Atmlabel(Natoms)
      Character*5 Atmlabel_zmtord(Nreals)

      Data Dnull,Izero /0.0D0,0/

      

      Do Iatmvml = 1, Natoms 
         Iatmzmt = Imap(Iatmvml)
         If (Iatmzmt .Ne. Izero) Then 
             Atmlabel_zmtord(Iatmzmt) = Atmlabel(Iatmvml)
             Iatmchrg_zmtord(Iatmzmt) = Iatmchrg(Iatmvml)
             Atmmass_zmtord(Iatmzmt)  = Atmmass(Iatmvml)
         Endif 
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(2a)") " The input ordered atom labels, nuclear charges",
     +                " and masses"
      write(6,"(6(4x,5a))") (Atmlabel_zmtord(i),I=1,Nreals)
      write(6,"(6(2x,I3))") (Iatmchrg_zmtord(i),I=1,Nreals)
      write(6,"(6(1x,F10.5))") (Atmmass_zmtord(i),I=1,Nreals)
      write(6,*) 
#endif

      Return
      End 


