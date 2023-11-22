      Subroutine RdOrd(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
      Implicit Integer (A-Z)
      Write(*,*) '@RDORD: This version of ACES was built without',
     &           ' a MOLCAS license;'
      Write(*,*) '        therefore, the keyword INTEGRALS cannot',
     &           ' take the value\n',
     &           '        "SEWARD".'
      Stop 1
      Return
      End
