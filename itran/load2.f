
c THIS ROUTINE LOADS AO INTEGRALS FROM THE INTEGRAL FILE IIJJ AND
c TRANSFORMS THE FIRST INDEX (PARTIAL TRANSFORMATION ONLY)

      SUBROUTINE LOAD2(CMO,W,W2,dBuf,iBuf,iSymAO,nBas,
     &                 nOcc,nStart,nEnd,iSize,iSizT,iOffAO,iOffI,
     &                 nSize,iLnBuf,iSpin,iUnit,bLast)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL bLast
c
      DIMENSION iBuf(iLnBuf),dBuf(iLnBuf),W(1),W2(1),CMO(1)
      DIMENSION nBas(8),nOcc(8),iSize(8),iOffAO(8),iOffI(8)
      DIMENSION iSizT(8,8),nStart(8),nEnd(8)
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

c   o increment Load2 counter
      nPass2 = nPass2 + 1

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

c      o get irrep of integral and offset within basis functions
         IrrepX1 = iSymAO(IX)
         IrrepX2 = iSymAO(KX)
         iOff1   = iOffAO(IrrepX1)
         iOff2   = iOffAO(IrrepX2)

c      o scale the indices
         IX = IX - iOff1
         JX = JX - iOff1
         KX = KX - iOff2
         LX = LX - iOff2

c      o condition the packed indices
         I = MAX(IX,JX)
         J = MIN(IX,JX)
         K = MAX(KX,LX)
         L = MIN(KX,LX)

c      o process the occupied orbitals within IrrepX
         nStartX1 = nStart(IrrepX1)
         nStartX2 = nStart(IrrepX2)
         nEndX1   = nEnd(IrrepX1)
         nEndX2   = nEnd(IrrepX2)

c      o determine number of basis functions and offsets within IrrepX
         iSizeX1 = iSize(IrrepX1)
         iSizeX2 = iSize(IrrepX2)
         nBasX1  = nBas(IrrepX1)
         nBasX2  = nBas(IrrepX2)
         iOffC1  = INDOCC(IrrepX1,iSpin) + nBasX1*(nStartX1-1) - 1
         iOffC2  = INDOCC(IrrepX2,iSpin) + nBasX2*(nStartX2-1) - 1
         iOffW1  = iOffI(IrrepX1)
         iOffW2  = iOffI(IrrepX2)

c DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
c THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
c
c   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
c   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
c
c HOWEVER, WE NEED ONLY FOUR SINCE WE STORE THE AS
c   NU, SIGMA >= RHO,  I

         IND = nBasX2 * ( iSizT(IrrepX2,IrrepX1) + I*(I-1)/2 + (J-1) )
         if (K.eq.L) then
            X1 = X * 0.5d0
         else
            X1 = X
         end if

c      o TRANSFORM INDEX L TO MO BASIS AND INCREMENT L(KX,IX,JX;I)
         iOffC2L = iOffC2 + LX
         iAdr    = iOffW2 + KX + IND
         do iOcc = nStartX2, nEndX2
            W(iAdr) = W(iAdr) + X1*CMO(iOffC2L)
            iOffC2L = iOffC2L + nBasX2
            iAdr    = iAdr    + iSizeX2
         end do

c      o TRANSFORM INDEX K TO MO BASIS AND INCREMENT L(LX,JX,IX;I)
         iOffC2K = iOffC2 + KX
         iAdr    = iOffW2 + LX + IND
         do iOcc = nStartX2, nEndX2
            W(iAdr) = W(iAdr) + X1*CMO(iOffC2K)
            iOffC2K = iOffC2K + nBasX2
            iAdr    = iAdr    + iSizeX2
         end do

         IND = nBasX1 * ( iSizT(IrrepX1,IrrepX2) + K*(K-1)/2 + (L-1) )
         if (I.eq.J) then
            X2 = X * 0.5d0
         else
            X2 = X
         end if

c      o TRANSFORM INDEX J TO MO BASIS AND INCREMENT L(IX,KX,LX;I)
         iOffC1J = iOffC1 + JX
         iAdr    = iOffW1 + IX + IND
         do iOcc = nStartX1, nEndX1
            W(iAdr) = W(iAdr) + X2*CMO(iOffC1J)
            iOffC1J = iOffC1J + nBasX1
            iAdr    = iAdr    + iSizeX1
         end do

c      o TRANSFORM INDEX I TO MO BASIS AND INCREMENT L(JX,LX,KX;I)
         iOffC1I = iOffC1 + IX
         iAdr    = iOffW1 + JX + IND
         do iOcc = nStartX1, nEndX1
            W(iAdr) = W(iAdr) + X2*CMO(iOffC1I)
            iOffC1I = iOffC1I + nBasX1
            iAdr    = iAdr    + iSizeX1
         end do

c     end do int = 1, NUT
      end do

      nAOInt = nAOInt + NUT

c     end do while (NUT.eq.iLnBuf)
      end do

      if (NUT.eq.-1) then
         nAOInt = nAOInt + 1
         NUT = 0
      end if
      nLoad2 = nAOInt

      if (bLast) then
         if (iFlags(18).ge.3.or.
     &       iFlags(26).eq.1.or.
     &       iFlags(93).eq.2.or.
     &       (iSpin.eq.1.and.iFlags(11).gt.0).or.
     &       iFlags2(103).eq.1
     &      ) then
            close(unit=iUnit,status='KEEP')
         else
            close(unit=iUnit,status='DELETE')
         end if
      else
         rewind(iUnit)
         call locate(iUnit,'TWOELSUP')
      end if

      return
c     end subroutine load2
      end

