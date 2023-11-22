#include "program.h"
      subroutine a2rd_vib(natoms,iatchrg,freq,freqco,dnormmd,Nimag,
     &                    Nvib,iunit,Write_molden_file)
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nelement=103)
      parameter (angtob=1.88972598772D0)
      character *80 wrk
      character *10 Label
      Character  *4 PTGRP
      Character *11 Type(3*Natoms)
      character*2 celemol(nelement)
      logical yesno,Write_molden_file,Not_done
      inTEger iatchrg(natoms)
      DOuble precision freq(nvib),freqco(3,natoms),
     &  dnormmd(nvib,3,natoms)

#include "machsp.com"
#include "flags.com"
      data (celemol(i),i=1,nelement)
     & /' H','He',
     & 'Li','Be',' B',' C',' N',' O',' F','Ne',
     & 'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     & ' K','Ca',
     & 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     & 'Ga','Ge','As','Se','Br','Kr',
     & 'Rb','Sr',
     & ' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     & 'In','Sn','Sb','Te',' I','Xe',
     & 'Cs','Ba','La',
     & 'Ce','Pr','Nd','Pm','Sm','Eu','Gd',
     & 'Tb','Dy','Ho','Er','Tm','Yb','Lu',
     & 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     & 'Tl','Pb','Bi','Po','At','Rn',
     & 'Fr','Ra','Ac',
     & 'Th','Pa',' U','Np','Pu','Am','Cm',
     & 'Bk','Cf','Es','Fm','Md','No','Lr'/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C If Writing the MOLDEN files then read the NORMCO file otherwise
C read the FRQARC file.
C
      If (Write_molden_file) Then
C
C Reset the number of frequncies to include rotational and
C vibrational modes since the NORMCO file contains all of
C them.
C
        Nvib  = 3*Natoms
        Nimag = 0
        INQUIRE(FILE='NORMCO',EXIST=YESNO)
        IF(YESNO)THEN
          OPEN(UNIT=4,FILE='NORMCO',FORM='FORMATTED',STATUS='OLD')
          rewind(4)
        else
          write(6,10)
10        format(T3,'@RDVIB-F, NORMCO file not present.')
          call errex
        ENDIF
C
C This is an error, molden visualization does not need mass weigthed 
C coordnates. Ajth Perera, 05/2010
        read(4,*)wrk
        do iatom=1,natoms
CSSS           read(4,*)(freqco(i,iatom),i=1,3)
           read(4,*) tmp
        end do 
      
        do 30 ivib=1,nvib
           read(4,*)wrk      
           read(4,'(a)') wrk(1:21)
           read(unit=wrk(1:20),fmt=*) freqy
           if (wrk(21:21) .eq. 'i') Then
              freqy = -freqy 
              Nimag = Nimag + 1
           endif
           freq(ivib) = freqy      
           read(4,*)wrk
           do 40 iatom=1,natoms
              read(4,*)(dnormmd(ivib,i,iatom),i=1,3)
#ifdef _DEBUG_LVLM1
              Write(6, "(F10.5)") (dnormmd(ivib,i,iatom),i=1,3)
#endif
         
   40      continue
   30   continue
C 
C Use the VIB_TYPE record to write frq. that correspond only to vibrations.
C Others are written as zeros. 06/2018.

#ifdef _USE_CFOUR
        Call getcrec(0,'JOBARC','VIB_TYPE',Length,Ijunk)
        If (Length .Le. 0) Then
        Write(6,"(2a)") " Please modify CFOUR to add a new JOBARC",
     &                 " record VIB_TYPE to joda (vib1.f) and rerun."
        Write(6,"(2a)") " The vibration type is already determined",
     &                  " in vib1.f and stored in the variable TYPE."
        Call Errex
        Else
        Call GETREC_C4(20,'JOBARC','VIB_TYPE',Nvib*IINTFP,TYPE)
        ENdif 
#else
        Call getcrec(20,'JOBARC','VIB_TYPE',Nvib*11*IINTFP,TYPE)
#endif 
C
        close(unit=4,status='KEEP')    

        do 50 ivib=1,nvib
           If (type(ivib) .eq." VIBRATION ") then
           if (freq(ivib) .lt. 0.0D0) freq(ivib)= -freq(ivib)
           write(iunit,60)' ',freq(ivib)
           else
           write(iunit,60)' ',0.0D0
           endif 
60         format(A,f9.4)
50      continue
        write(iunit,70)'[FR-COORD]'
70      format(A)
        do 80 iatom=1,natoms
           write(iunit,90)' ',celemol(iatchrg(iatom)),freqco(1,iatom),
     &     freqco(2,iatom),freqco(3,iatom) 
90         format(A,A2,3x,f9.6,3x,f9.6,3x,f9.6)
80      continue
        write(iunit,100)'[FR-NORM-COORD]'
100     format(A)
        do 110 ivib=1,nvib
           write(iunit,115)'vibration    ',ivib
115        format(A,I2)
           do 120 iatom=1,natoms
              write(iunit,130)'   ',dnormmd(ivib,1,iatom),
     &        dnormmd(ivib,2,iatom),dnormmd(ivib,3,iatom)
130           format(A,f9.6,4x,f9.6,4x,f9.6)
120        continue 
110     continue
C
      Else
C
        INQUIRE(FILE='FRQARC',EXIST=YESNO)
        IF(YESNO)THEN
          OPEN(UNIT=5,FILE='FRQARC',FORM='UNFORMATTED',STATUS='OLD')
          rewind(5)
        else
          write(6,20)
20        format(T3,'@RDVIB-F, The FRQARC file not present.')
          call errex
        ENDIF
        READ(5) LABEL, PTGRP, NRX, NIMAG, ZMAS, ZIX, ZIY, ZIZ,
     &          (Freq(IVIB), IVIB=1, Nvib)
        Nimag = 0
        Nvib_read  = 3*Natoms 

C where the molecule assume principle axis orientations, the sym=none
C Cartesian frequency calculations (ACES III) this is not assured. 
C The following logic is suppose to take care of that issue. 
C Ajith Perera 12/2019.
 
        INQUIRE(FILE='NORMCO',EXIST=YESNO)
        IF(YESNO)THEN
          OPEN(UNIT=4,FILE='NORMCO',FORM='FORMATTED',STATUS='OLD')
          rewind(4)
        else
          write(6,10)
          call errex
        ENDIF

        Not_done = .True.
        Do While (Not_Done) 
           Read(4,*) Wrk
           If (Wrk(1:9) .EQ. "Projected") Not_done = .False.
        Enddo 

        do ivib=1,Nvib_read
           read(4,'(a)') wrk(1:21)
           read(unit=wrk(1:20),fmt=*) freqy
           if (wrk(21:21) .eq. 'i') Then
              freqy = -freqy
           endif
           freq(ivib) = freqy
         Enddo 

         Call getcrec(20,'JOBARC','VIB_TYPE',Nvib_read*11*IINTFP,TYPE)

         Ikeep = 0
         do ivib=1,Nvib_read
           If (type(ivib) .eq." VIBRATION ") then
              Ikeep = Ikeep + 1
              If (freq(ivib) .lt. 0.0D0) Then 
                  Nimag = Nimag + 1
                  freq(Ikeep)= freq(Ivib)
              Endif 
              Freq(Ikeep) = Freq(Ivib)
           Endif   
         Enddo 

#ifdef _DEBUG_LVLM
        Write(6,"(a)") "TGRP, NRX, NIMAG, ZMAS, ZIX, ZIY, ZIZ"
        Print*,LABEL, PTGRP, NRX, NIMAG, ZMAS, ZIX, ZIY, ZIZ
#endif
        Write(6,*) 
        Write(6,"(2a)") " Vibrational frequencies read from",
     &                  " NORMCO file."
        Write(6,"(5(1x,F15.5))") (Freq(IVIB), IVIB=1,Ikeep )

       If (Ikeep .NE. Nvib) Then
            Write(6,"(a)") " There is an internal inconsistency in",
     &                     " in the number of vibrations" 
            Write(6,"(a,I3)") "The number of vibrations read    :", 
     &                         Ikeep
            Write(6,"(a,I3)") "The number of vibrations computed:", 
     &                         Ivib
            Call Errex 
       Endif 
C The frequencies in FRQARC are before vibrational projections.
C While this should not have created problems for frequency calculations
C Scale coordinates to bohr from angstrom, do not know what
C this is at the moment, Ajith Perera, 04/2011
C
CSSSS        call dscal(3*natoms,angtob,coord,1)
        close(unit=4,status='KEEP')    
        close(unit=5,status='KEEP')
C
      Endif
C
      Return
      End
