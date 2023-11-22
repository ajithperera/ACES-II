
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
      COMMON /VTINFO/ nPass1,nPass2,nPass3,nPass4,
     &                nLoad1,nLoad2,nLoad3,nLoad4,
     &                NWRIT1, NWRIT2, NWRIT3, NWRIT4,
     &                NWRIT1A,NWRIT2A,NWRIT3A,NWRIT4A,
     &                NWRIT1B,NWRIT2B,NWRIT3B,NWRIT4B

c ----------------------------------------------------------------------

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

c      o process the occupied orbitals within IrrepX (if any)
         nStartX = nStart(IrrepX)
         nEndX   = nEnd(IrrepX)
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

c         o scale the integral
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

c         o determine number of basis functions and offsets within IrrepX
            iSizeX = iSize(IrrepX)
            nBasX  = nBas(IrrepX)
            iOffC  = INDOCC(IrrepX,iSpin) + nBasX*(nStartX-1) - 1
            iOffW  = iOffI(IrrepX)

c DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
c THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
c
c   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
c   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
c
c HOWEVER, WE NEED ONLY FOUR SINCE WE STORE THE AS
c   NU, SIGMA >= RHO,  I

            IND = nBasX * ( I*(I-1)/2 + (J-1) )

c         o TRANSFORM INDEX L TO MO BASIS AND INCREMENT L(KX,IX,JX;I)
            iOffCL = iOffC + LX
            iAdr   = iOffW + KX + IND
            do iCount = nStartX, nEndX
               W2(iCount) = X1*CMO(iOffCL)
               iOffCL = iOffCL + nBasX
C               iOffCL = iOffCL + 1
C               iAdr   = iAdr   + iSizeX
            end do
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + W2(ICount) 
C               iOffCL = iOffCL + nBasX
C               iOffCL = iOffCL + 1
               iAdr   = iAdr   + iSizeX
            end do

c         o TRANSFORM INDEX K TO MO BASIS AND INCREMENT L(LX,JX,IX;I)
            iOffCK = iOffC + KX
            iAdr   = iOffW + LX + IND
            do iCount = nStartX, nEndX
               W2(iCount) = X1*CMO(iOffCK)
               iOffCK = iOffCK + nBasX
C               iOffCK = iOffCK +  1
C               iAdr   = iAdr   + iSizeX
            end do
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + W2(iCount)
               iOffCK = iOffCK + nBasX
C               iOffCK = iOffCK +  1
               iAdr   = iAdr   + iSizeX
            end do

            IND = nBasX * ( K*(K-1)/2 + (L-1) )

c         o TRANSFORM INDEX J TO MO BASIS AND INCREMENT L(IX,KX,LX;I)
            iOffCJ = iOffC + JX
            iAdr   = iOffW + IX + IND
            do iCount = nStartX, nEndX
               W2(iCount) =  X2*CMO(iOffCJ)
               iOffCJ = iOffCJ + nBasX
C               iOffCJ = iOffCJ + 1
C               iAdr   = iAdr   + iSizeX
            end do
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + W2(iCount)
C               iOffCJ = iOffCJ + nBasX
C               iOffCJ = iOffCJ + 1
               iAdr   = iAdr   + iSizeX
            end do

c         o TRANSFORM INDEX I TO MO BASIS AND INCREMENT L(JX,LX,KX;I)
            iOffCI = iOffC + IX
            iAdr   = iOffW + JX + IND
            do iCount = nStartX, nEndX
               W2(iCount) = X2*CMO(iOffCI)
               iOffCI = iOffCI + nBasX
C               iOffCI = iOffCI + 1
C               iAdr   = iAdr   + iSizeX
            end do
            do iCount = nStartX, nEndX
               W(iAdr) = W(iAdr) + W2(ICount)
C               iOffCI = iOffCI + nBasX
C               iOffCI = iOffCI + 1
               iAdr   = iAdr   + iSizeX
            end do

c        end if (nStartX.le.nEndX)
         end if

c     end do int = 1, NUT
      end do

      nAOInt = nAOInt + NUT

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

c      write(*,'(t3,a,i8,a)')
c                   '@LOAD1: ',nAOInt,' integrals read in from file IIII'

      return
c     end subroutine load1
      end

