
c THIS ROUTINE LOADS AO INTEGRALS FROM THE INTEGRAL FILE IIII AND
c TRANSFORMS THE FIRST INDEX.

      SUBROUTINE LOAD1(CMO,W,W2,dBuf,iBuf,iSymAO,nBas,
     &                 nFirst,nStart,nEnd,iSize,iOffAO,iOffI,
     &                 nSize,iLnBuf,iSpin,iUnit,bLast)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL bLast
c
      DIMENSION iBuf(iLnBuf),dBuf(iLnBuf),W(1),W2(1),CMO(1)
      DIMENSION nBas(8),nFirst(8),iSize(8),iOffAO(8),iOffI(8)
      DIMENSION nStart(8),nEnd(8)
      DIMENSION iSymAO(100)
c
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,iAlone,iBitWd
      COMMON /FLAGS/  iFlags(100)
      COMMON /FLAGS2/ iFlags2(500)
      COMMON /AOOFST/ INDOCC(8,2)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /VTINFO/ nPass1,nPass2,nPass3,nPass4,
     &                nLoad1,nLoad2,nLoad3,nLoad4,
     &                NWRIT1, NWRIT2, NWRIT3, NWRIT4,
     &                NWRIT1A,NWRIT2A,NWRIT3A,NWRIT4A,
     &                NWRIT1B,NWRIT2B,NWRIT3B,NWRIT4B

c ----------------------------------------------------------------------
      Indx(i, j) = j+(i*(i-1))/2
      NnP1o2(I)  = (i*(i+1))/2

#ifdef _DEBUG_LVL0
      Print*, "The nstartx and nendx arrays"
      Print*, (nStart(I), i=1, 8)
      Print*, (nEnd(I), i=1, 8)
      Write(6,*)
      Print*, "The Offset for AO and MO"
      Print*, (iOffAO(I), i=1, 8)
      Print*, (iOffI(I), i=1, 8)
      Write(6,*) 
      Print*, "The Nbas Nfirst and Isize arrays"
      Print*, (Nbas(I), i=1, 8)
      Print*, (NFirst(I), i=1, 8)
      Print*, (isize(I), i=1, 8)
      Write(6,*)
      Print*, "The Nsize", nSize
#endif
   

c   o initialize W array
      call zero(W,nSize)

c   o increment Load1 counter
      nPass1 = nPass1 + 1

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
         Print*, Irrepx

c      o process the occupied orbitals within IrrepX (if any)

#ifdef _NONOSKIP
         nStartX = nStart(IrrepX)
         nEndX   = nEnd(IrrepX)
#endif
C
         if (nStartX.le.nEndX) then

c         o get offsets within basis functions
            iTmp   = iOffAO(IrrepX)
            IX     = IX - iTmp
            JX     = JX - iTmp
            KX     = KX - iTmp
            LX     = LX - iTmp

c         o condition the packed indices
            I = max(IX,JX)
            J = min(IX,JX)
            K = max(KX,LX)
            L = min(KX,LX)
#ifdef _DEBUG_LVL0
         write(*,"(t3,a,1x,I3,1x,i3,1x,i3,1x,i3)"),
     &         "IX,JX,KX,LX ", IX,JX, KX, LX
         write(*,"(t3,a,1x,I3,1x,i3,1x,i3,1x,i3)"),
     &         "I,J,K,L ", I,J, K, L
#endif
C
c         o scale the integral
#ifdef _NONOSKIP
            if (I.eq.K.and.J.eq.L) X = X * 0.5d0
            if (K.eq.L) then
               X1 = X * 0.5d0
            else
               X1 = X
            end if
            if (I.eq.J) then
               X2 = X * 0.5d0
            else
               X2 = X
            end if
#endif

c         o determine number of basis functions and offsets within IrrepX
C
#ifdef _NONOSKIP
            iSizeX = iSize(IrrepX)
            nBasX  = nBas(IrrepX)
            iOffC  = INDOCC(IrrepX,iSpin) + nBasX*(nStartX-1) - 1
            iOffW  = iOffI(IrrepX)
            iOffW  = iOffAO(IrrepX)
#endif

c DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
c THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
c
c   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
c   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
c
            INDKL = Indx(K, L)

#ifdef _NONOSKIP
#ifdef _DEBUG_LVL0
            write(*,"(t3,a,1x,I5,1x,i5)"),
     &           "Nstart and Nend", NstartX, NendX
#endif
#endif
C
            nBasX  = nBas(IrrepX)
            iAdrIJ = I*(I-1)/2 + J
            iAdr   = IAdrIJ +  NnP1o2(nBasX)*(INDKL - 1)
C
#ifdef _DEBUG_LVL0
             write(*,"(t3,a,1x,I5,1X,F12.7)"),"AdrI-AdrJ and X",
     &             iAdr, X
             Write(6,*)
#endif
C
#ifdef _READ_ONLY
                 W(iAdr) = X
#endif
c
#ifdef _NONOSKIP
c   NU, SIGMA >= RHO,  I

            IND = nBasX * ( I*(I-1)/2 + (J-1) )
C
C Given pair of I and J (mu,nu), go over K and L indices
C and compute Int(mu,nu;lam,sig)*MO(lam,i) and Int(mu,nu;lam,sgi)*
C MO(sig,i) and accumulate them in iAdr which is KX +IND or LX + IND
C advanced by  nmo*(nmo+1)/2 for each mo (nmo: #of mos in a irrep).
C So, for a given Int, we proccess 2*nmo Ints (addresses) (one for each K,L)
C and the transformed elements have INt*CMO(lam,i) and INt*CMO(sg,i) and
C written in to nmo address locations for each K and K.. As the file get
C read, for given I, J all lam and sig get transformed and added to the same 
C nmo locations. So, at the end we have, 2*nmo addreses locations for
C given I, J.

c         o TRANSFORM INDEX L TO MO BASIS AND INCREMENT L(KX,IX,JX;I)
            iOffCL = iOffC + LX
            iAdr   = iOffW + KX + IND
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + X1*CMO(iOffCL)
#ifdef _DEBUG_LVL0
               Write(*,'(t3,a,I100,3(1X,F12.7))') "Transforming L", 
     &              iAdr, X1, CMO(iOffCL), W(iadr)
               Write(6,*)
#endif 
               iOffCL = iOffCL + nBasX
               iAdr   = iAdr   + iSizeX
            end do
            
c         o TRANSFORM INDEX K TO MO BASIS AND INCREMENT L(LX,JX,IX;I)
            iOffCK = iOffC + KX
            iAdr   = iOffW + LX + IND
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + X1*CMO(iOffCK)
#ifdef _DEBUG_LVL0
               Write(*,'(t3,a,I100,3(1X,F12.7))') "Transforming K", 
     &              iAdr, X1, CMO(iOffCK), W(iadr)
               Write(6,*)
#endif 
               iOffCK = iOffCK + nBasX
               iAdr   = iAdr   + iSizeX
            end do
            
            IND = nBasX * ( K*(K-1)/2 + (L-1) )
C
C Given pair of K and L (lam,sig), go over I and J indices
C and compute Int(mu,nu;lam,sig)*MO(mu,i) and Int(mu,nu;lam,sgi)*
C MO(nu,i) and accumulate them in iAdr which is IX +IND or JX + IND
C advanced by  nmo*(nmo+1)/2 for each mo (nmo: #of mos in a irrep).
C So, for a given Int, we proccess 2*nmo Ints (addresses) (one for each I, J)
C and the transformed elements have INt*CMO(mu,i) and INt*CMO(nu,i) and
C written in to nmo address locations for each K and L. As the file get
C read, for given K, L all mu and nu get transformed and added to the same
C nmo locations. So, at the end we have, 2*nmo addreses locations for
C given K, L.
            
c         o TRANSFORM INDEX J TO MO BASIS AND INCREMENT L(IX,KX,LX;I)
            iOffCJ = iOffC + JX
            iAdr   = iOffW + IX + IND
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + X2*CMO(iOffCJ)
#ifdef _DEBUG_LVL0
               Write(*,'(t3,a,I100,3(1X,F12.7))') "Transforming J", 
     &              iAdr, X2, CMO(iOffCJ), W(iadr)
               Write(6,*)
#endif 
               iOffCJ = iOffCJ + nBasX
               iAdr   = iAdr   + iSizeX
            end do
            
c         o TRANSFORM INDEX I TO MO BASIS AND INCREMENT L(JX,LX,KX;I)
            iOffCI = iOffC + IX
            iAdr   = iOffW + JX + IND
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + X2*CMO(iOffCI)
#ifdef _DEBUG_LVL0
               Write(*,'(t3,a,I100,3(1X,F12.7))') "Transforming I", 
     &              iAdr, X2, CMO(iOffCI), W(iadr)
               Write(6,*)
#endif 
               iOffCI = iOffCI + nBasX
               iAdr   = iAdr   + iSizeX
            end do
         Write(6,*)
C
C Since IJ are faster running than KL, the final first indexed transformed
C integrals are orderd as for given K>L all I>J. 
C
#endif
C
c        end if (nStartX.le.nEndX)
         end if

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
         Nmo = NFirst(IrreP)
         KLPAIRS = NAO*(NAO+1)/2
         IJPAIRS = NAO*(NAO+1)/2
#ifdef _NONOSKIP
         IJPAIRS = NAO*NMO
#endif
         print*, IJPAIRS
         Do KL = 1, KLPAIRS
            Print*, "The KL Pair = ", KL
            Write(*,'(4(1X,F12.7)))'), (W(Ioff + I), I=1, IJPAIRS)
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
     &           '@LOAD1: ',nAOInt,' integrals read in from 
     &           file IIII'

      return
c     end subroutine load1
      end

