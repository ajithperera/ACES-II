
      SUBROUTINE READ_AO_INTS(W,dBuf,iBuf,iSymAO,iOffAO,
     &                        iOff_Irrep_4AO,nBas, nSize,
     &                        iLnBuf,iUnit,bLast,DoKL,DoIJ)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL bLast, DoKL, DoIJ, Full_diagonal, Half_diagonal
c
      DIMENSION iBuf(iLnBuf),dBuf(iLnBuf),W(1)
      DIMENSION nBas(8),iOffAO(8),iOff_Irrep_4AO(8)
      DIMENSION iSymAO(100)
c
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,iAlone,iBitWd
      COMMON /FLAGS/  iFlags(100)
      COMMON /FLAGS2/ iFlags2(500)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)

c ----------------------------------------------------------------------
      Indx(i, j) = j+(i*(i-1))/2
      NnP1o2(I)  = (i*(i+1))/2

c   o initialize W array
      If (DoKL) call zero(W,nSize)

c   o increment Load1 counter
C
      nAOInt = 0
      NUT    = iLnBuf
      do while (NUT.eq.iLnBuf)

c   o read in integrals
      read(iUnit) dBuf, iBuf, NUT
      do int = 1, NUT

c      o extract X value and IX, JX, KX, LX indices
         X    = dBuf(int)
         iTmp = iBuf(int)
         IX   = iand(      iTmp,           iAlone)
         JX   = iand(ishft(iTmp,-  iBitWd),iAlone)
         KX   = iand(ishft(iTmp,-2*iBitWd),iAlone)
         LX   = iand(ishft(iTmp,-3*iBitWd),iAlone)

c      o get irrep of integral
         IrrepX = iSymAO(IX)
C
c      o get offsets within basis functions
         iTmp   = iOffAO(IrrepX)
         IX     = IX - iTmp
         JX     = JX - iTmp
         KX     = KX - iTmp
         LX     = LX - iTmp

c      o condition the packed indices
         I = max(IX,JX)
         J = min(IX,JX)
         K = max(KX,LX)
         L = min(KX,LX)
         Full_diagonal = (I.eq.j) .and. (I.eq.k) .and. (i.eq.l) 
     &                    .and. (j.eq.k) .and. (j.eq.l) .and. 
     &                    (k.eq.l)
         Half_diagonal = (i.eq.k) .and. (j.eq.l)
C
#ifdef _DEBUG_LVLM1
         write(*,"(t3,a,1x,I3,1x,i3,1x,i3,1x,i3)"),
     &         "IX,JX,KX,LX ", IX,JX, KX, LX
         write(*,"(t3,a,1x,I3,1x,i3,1x,i3,1x,i3)"),
     &         "I,J,K,L ", I,J, K, L
#endif
C
          IF (DoKL) THEN
          IOff_AO = iOff_Irrep_4AO(Irrepx)
          INDKL = Indx(K, L)
          nBasX  = nBas(IrrepX)
          iAdrIJ = I*(I-1)/2 + J
          iAdr   = IAdrIJ +  NnP1o2(nBasX)*(INDKL - 1)
          W(iAdr+Ioff_AO) = X
C
#ifdef _DEBUG_LVLM1
          write(*,"(t3,a,1x,I5,1X,F12.7)"),"AdrI-AdrJ and X",
     &             iAdr, X
          Write(6,*)
#endif
          Endif

          IF (DoIJ) Then
          IOff_AO = iOff_Irrep_4AO(Irrepx)
          INDIJ = Indx(I, J)
          nBasX  = nBas(IrrepX)
          iAdrKL = K*(K-1)/2 + L 
          iAdr   = IAdrKL +  NnP1o2(nBasX)*(INDIJ - 1)
          If (.NOT. (Full_diagonal .or.  Half_diagonal))
     &       W(iAdr+Ioff_AO) = W(iAdr+Ioff_AO) + X
#ifdef _DEBUG_LVLM1
          write(*,"(t3,a,1x,I5,1X,F12.7)"),"AdrK-AdrL and X",
     &             iAdr, X
          Write(6,*)
#endif
           Endif
C
c     end do int = 1, NUT
      end do

      nAOInt = nAOInt + NUT
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Print*, "The AO/MO(if old way) integrals"
      Ioff = 0
      Do Irrep = 1, Nirrep
         Nao = NBAS(IrreP)
         KLPAIRS = NAO*(NAO+1)/2
         IJPAIRS = NAO*(NAO+1)/2
C
         Do KL = 1, KLPAIRS
            Print*, "The KL Pair = ", KL
            Write(*,'(4(1X,F12.7)))'), (W(Ioff + I), I=1, IJPAIRS)
CSSS            Call checksum2("AOLAD1", W(ioff), IJPAIRS)
            Ioff = Ioff + IJPAIRS
         End Do
      End Do 
#endif

c     end do while (NUT.eq.iLnBuf)
      end do

      if (NUT.eq.-1) then
         nAOInt = nAOInt + 1
         NUT = 0
      end if
      nLoad1 = nAOInt

      if (bLast) then
         if (iFlags(18).ge.3.or.
     &       iFlags(26).eq.1.or.
     &       iFlags(93).eq.2.or.
     &       (iSpin.eq.1.and.iFlags(11).gt.0).or.
     &       iFlags2(103).eq.1
     &      ) then
            close(unit=iUnit,status='KEEP')
         else
            close(unit=iUnit,status='KEEP')
         end if
      else
         rewind(iUnit)
         call locate(iUnit,'TWOELSUP')
      end if

      write(*,'(t3,a,i8,a)')
     &           '@READ_AO_INTS: ',nAOInt,' integrals read in from 
     &           file IIII'

      return
c     end subroutine load1
      end

