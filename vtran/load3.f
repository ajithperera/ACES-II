
c THIS ROUTINE LOADS AO INTEGRALS FROM THE INTEGRAL FILE IJIJ AND
c TRANSFORMS THE FIRST INDEX (PARTIAL TRANSFORMATION ONLY)

      SUBROUTINE LOAD3(CMO,W,W2,dBuf,iBuf,iSymAO,nBas,
     &                 nOcc,nStart,nEnd,iOffAO,iOff4,iSize3,
     &                 nSize,iLnBuf,iSpin,iUnit,bLast,iOffset)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL bLast
c
      DIMENSION iBuf(iLnBuf),dBuf(iLnBuf),W(1),W2(1),CMO(1)
      DIMENSION nBas(8),nOcc(8),iOffAO(8),iOff4(8,8)
      DIMENSION iSize3(8,8),nStart(8),nEnd(8)
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

c   o increment Load3 counter
      nPass3 = nPass3 + 1

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
         IrrepX2 = iSymAO(JX)
         IrrepX3 = iSymAO(KX)
         IrrepX4 = iSymAO(LX)
         if (IrrepX3.ne.IrrepX1) then
            iTmp = LX
            LX   = KX
            KX   = iTmp
         end if
         if (IrrepX1.gt.IrrepX2) then
            iTmp = IX
            IX   = JX
            JX   = iTmp
            iTmp = KX
            KX   = LX
            LX   = iTmp
            iTmp    = IrrepX1
            IrrepX1 = IrrepX2
            IrrepX2 = iTmp
         end if
         iOff1 = iOffAO(IrrepX1)
         iOff2 = iOffAO(IrrepX2)

c      o scale the indices
         I = IX - iOff1
         J = JX - iOff2
         K = KX - iOff1
         L = LX - iOff2

c      o process the occupied orbitals within IrrepX
         nStartX1 = nStart(IrrepX1)
         nStartX2 = nStart(IrrepX2)
         nEndX1   = nEnd(IrrepX1)
         nEndX2   = nEnd(IrrepX2)

c      o determine number of basis functions and offsets within IrrepX
         iSize32 = iSize3(IrrepX2,IrrepX1)
         iSize31 = iSize3(IrrepX1,IrrepX2)
         nBasX1  = nBas(IrrepX1)
         nBasX2  = nBas(IrrepX2)
         iOffC1  = INDOCC(IrrepX1,iSpin) + nBasX1*(nStartX1-1) - 1
         iOffC2  = INDOCC(IrrepX2,iSpin) + nBasX2*(nStartX2-1) - 1
         iOffW42 = iOff4(IrrepX2,IrrepX1)
         iOffW41 = iOff4(IrrepX1,IrrepX2)

c      o scale the integral
         if (I.eq.K.and.J.eq.L) X = X * 0.5d0

c DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
c THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
c
c   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
c   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
c
c HOWEVER, WE NEED ONLY FOUR SINCE WE STORE THE AS
c   NU, SIGMA >= RHO,  I

         IND = nBasX1*(J-1) + (I-1)

c      o TRANSFORM INDEX L TO MO BASIS AND INCREMENT L(KX,IX,JX;I)
         iOffC2L = iOffC2  + L
         iAdr    = iOffW42 + K + nBasX1*IND
         do iOcc = nStartX2, nEndX2
            W(iAdr) = W(iAdr) + X*CMO(iOffC2L)
            iOffC2L = iOffC2L + nBasX2
            iAdr    = iAdr    + iSize32
         end do

c      o TRANSFORM INDEX K TO MO BASIS AND INCREMENT L(LX,JX,IX;I)
         iOffC1K = iOffC1  + K
         iAdr    = iOffW41 + L + nBasX2*IND
         do iOcc = nStartX1, nEndX1
            W(iAdr) = W(iAdr) + X*CMO(iOffC1K)
            iOffC1K = iOffC1K + nBasX1
            iAdr    = iAdr    + iSize31
         end do

         IND = nBasX1*(L-1) + (K-1)

c      o TRANSFORM INDEX J TO MO BASIS AND INCREMENT L(IX,KX,LX;I)
         iOffC2J = iOffC2  + J
         iAdr    = iOffW42 + I + nBasX1*IND
         do iOcc = nStartX2, nEndX2
            W(iAdr) = W(iAdr) + X*CMO(iOffC2J)
            iOffC2J = iOffC2J + nBasX2
            iAdr    = iAdr    + iSize32
         end do

c      o TRANSFORM INDEX I TO MO BASIS AND INCREMENT L(JX,LX,KX;I)
         iOffC1I = iOffC1  + I
         iAdr    = iOffW41 + J + nBasX2*IND
         do iOcc = nStartX1, nEndX1
            W(iAdr) = W(iAdr) + X*CMO(iOffC1I)
            iOffC1I = iOffC1I + nBasX1
            iAdr    = iAdr    + iSize31
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
      nLoad3 = nAOInt

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
c     end subroutine load3
      end

