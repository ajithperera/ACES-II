
c THIS ROUTINE LOADS AO INTEGRALS FROM THE INTEGRAL FILE IJKL AND
c TRANSFORMS THE FIRST INDEX (PARTIAL TRANSFORMATION ONLY)

      SUBROUTINE LOAD4(CMO,W,W2,dBuf,iBuf,iSymAO,nBas,
     &                 nOcc,nStart,nEnd,iOffAO,iOf4,iSize3,iSize2,
     &                 nSize,iLnBuf,iSpin,iUnit,bLast,iOffset)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL bLast
c
      DIMENSION iBuf(iLnBuf),dBuf(iLnBuf),W(1),W2(1),CMO(1)
      DIMENSION nBas(8),nOcc(8),iOffAO(8),iOf4(8,8),iSize2(8,8,8)
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

c   o increment Load4 counter
      nPass4 = nPass4 + 1

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
         iOff1   = iOffAO(IrrepX1)
         iOff2   = iOffAO(IrrepX2)
         iOff3   = iOffAO(IrrepX3)
         iOff4   = iOffAO(IrrepX4)

c      o scale the indices
         I = IX - iOff1
         J = JX - iOff2
         K = KX - iOff3
         L = LX - iOff4

c      o process the occupied orbitals within IrrepX
         nStartX1 = nStart(IrrepX1)
         nStartX2 = nStart(IrrepX2)
         nStartX3 = nStart(IrrepX3)
         nStartX4 = nStart(IrrepX4)
         nEndX1   = nEnd(IrrepX1)
         nEndX2   = nEnd(IrrepX2)
         nEndX3   = nEnd(IrrepX3)
         nEndX4   = nEnd(IrrepX4)

c      o determine number of basis functions and offsets within IrrepX
         iSize32 = iSize3(IrrepX2,IrrepX1)
         iSize31 = iSize3(IrrepX1,IrrepX2)
         iSize34 = iSize3(IrrepX4,IrrepX3)
         iSize33 = iSize3(IrrepX3,IrrepX4)
         nBasX1  = nBas(IrrepX1)
         nBasX2  = nBas(IrrepX2)
         nBasX3  = nBas(IrrepX3)
         nBasX4  = nBas(IrrepX4)
         iOffC1  = INDOCC(IrrepX1,iSpin) + nBasX1*(nStartX1-1) - 1
         iOffC2  = INDOCC(IrrepX2,iSpin) + nBasX2*(nStartX2-1) - 1
         iOffC3  = INDOCC(IrrepX3,iSpin) + nBasX3*(nStartX3-1) - 1
         iOffC4  = INDOCC(IrrepX4,iSpin) + nBasX4*(nStartX4-1) - 1
         iOffW42 = iOf4(IrrepX2,IrrepX1)
         iOffW41 = iOf4(IrrepX1,IrrepX2)
         iOffW44 = iOf4(IrrepX4,IrrepX3)
         iOffW43 = iOf4(IrrepX3,IrrepX4)

c DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
c THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
c
c   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
c   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
c
c HOWEVER, WE NEED ONLY FOUR SINCE WE STORE THE AS
c   NU, SIGMA >= RHO,  I

         if (IrrepX1.lt.IrrepX2) then
            IND = nBasX1*(J-1) + (I-1) + iSize2(IrrepX3,IrrepX1,IrrepX2)
         else
            IND = nBasX2*(I-1) + (J-1) + iSize2(IrrepX3,IrrepX1,IrrepX2)
         end if

c      o TRANSFORM INDEX L TO MO BASIS AND INCREMENT L(KX,IX,JX;I)
         iOffC4L = iOffC4  + L
         iAdr    = iOffW44 + K + nBasX3*IND
         do iOcc = nStartX4, nEndX4
            W(iAdr) = W(iAdr) + X*CMO(iOffC4L)
            iOffC4L = iOffC4L + nBasX4
            iAdr    = iAdr    + iSize34
         end do

c      o TRANSFORM INDEX K TO MO BASIS AND INCREMENT L(LX,JX,IX;I)
         iOffC3K = iOffC3  + K
         iAdr    = iOffW43 + L + nBasX4*IND
         do iOcc = nStartX3, nEndX3
            W(iAdr) = W(iAdr) + X*CMO(iOffC3K)
            iOffC3K = iOffC3K + nBasX3
            iAdr    = iAdr    + iSize33
         end do

         if (IrrepX3.lt.IrrepX4) then
            IND = nBasX3*(L-1) + (K-1)
         else
            IND = nBasX4*(K-1) + (L-1)
         end if

c      o TRANSFORM INDEX J TO MO BASIS AND INCREMENT L(IX,KX,LX;I)
         iOffC2J = iOffC2  + J
         iAdr    = iOffW42 + I + nBasX1*(
     &                           IND + iSize2(IrrepX1,IrrepX3,IrrepX4) )
         do iOcc = nStartX2, nEndX2
            W(iAdr) = W(iAdr) + X*CMO(iOffC2J)
            iOffC2J = iOffC2J + nBasX2
            iAdr    = iAdr    + iSize32
         end do

c      o TRANSFORM INDEX I TO MO BASIS AND INCREMENT L(JX,LX,KX;I)
         iOffC1I = iOffC1  + I
         iAdr    = iOffW41 + J + nBasX2*(
     &                           IND + iSize2(IrrepX2,IrrepX3,IrrepX4) )
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
      nLoad4 = nAOInt

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
c     end subroutine load4
      end

