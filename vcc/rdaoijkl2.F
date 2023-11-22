
c This routine loads the AO integrals from the IJKL file (all four indices have
c different symmetries) and contracts them with the T2 amplitudes.

      subroutine rdaoijkl2(T2,dZ,
     &                     dBuf,iBuf,
     &                     IK0,IL0,
     &                     JK0,JL0,
     &                     KI0,KJ0,
     &                     dVal,
     &                     iSym,iType,iWhere,iAOSym,iMap,
     &                     iLnBuf,iUnit,iUHF,nAO,iMax,irp_x,
     &                     iBchOff,iBchLen)
      implicit none

c ARGUMENTS
      integer iMax, iLnBuf
      double precision T2(iMax), dZ(iMax), dBuf(iLnBuf)
      integer iBuf(iLnBuf)
      integer IK0(8*iLnBuf), IL0(8*iLnBuf)
      integer JK0(8*iLnBuf), JL0(8*iLnBuf)
      integer KI0(8*iLnBuf), KJ0(8*iLnBuf)
      double precision dVal(8*iLnBuf)
      integer iSym(8*iLnBuf), iType(8*iLnBuf), iWhere(8*iLnBuf)
      integer iAOSym(*), iMap(*)
      integer iUnit, iUHF, nAO, irp_x
      integer iBchOff(8), iBchLen(8)

c INTERNAL VARIABLES
      double precision dX
      integer ndx, iOff,  iOff1, iOff2
      integer IKOff, LenIK, lenik2, iRest
      integer iSpin
      integer irp_ij, irp_ik, irp_il, irp_ab
      integer i, ndx_i, iTmp,     ij, ik, il
      integer j, ndx_j, jTmp,         jk, jl
      integer k, ndx_k, kTmp, ki, kj,     kl
      integer l, ndx_l, lTmp, li, lj
      integer n
      integer int, nMatch, NUT, iCount

      INTEGER        IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &               IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON /AOSYM/ IAOPOP,   IOFFAO,   IOFFV,     IOFFO,
     &               IRPDPDAO,   IRPDPDAOS,   ISTART,     ISTARTMO
#include "machsp.com"
#include "sympop.com"
#include "syminf.com"
#include "flags.com"

      integer iupki, iupkj, iupkk, iupkl, indx
      IUPKI(iTmp)=iand(      iTmp,           IALONE)
      IUPKJ(iTmp)=iand(ishft(iTmp,-  IBITWD),IALONE)
      IUPKK(iTmp)=iand(ishft(iTmp,-2*IBITWD),IALONE)
      IUPKL(iTmp)=iand(ishft(iTmp,-3*IBITWD),IALONE)
      INDX(I,J)=J+(I*(I-1))/2

c ----------------------------------------------------------------------

      call locate(iUnit,'TWOELSUP')

c ----------------------------------------------------------------------

      if (iUHF.ne.0) then

         NUT = 0
         do while (NUT.ne.-1)

            read(iUnit) dBuf, iBuf, NUT

            iCount = 0
            do int = 1, NUT

               ndx_i = IUPKI(iBuf(int))
               ndx_j = IUPKJ(iBuf(int))
               ndx_k = IUPKK(iBuf(int))
               ndx_l = IUPKL(iBuf(int))

               ij = INDX(max(ndx_j,ndx_i),min(ndx_j,ndx_i))
               kl = INDX(max(ndx_l,ndx_k),min(ndx_l,ndx_k))

               if (ndx_i.ne.ndx_j) then

                  if (ndx_k.ne.ndx_l) then

                     if (ij.ne.kl) then

                        irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))
                        irp_il = dirprd(iAOSym(ndx_i),iAOSym(ndx_l))

                                     iTmp = nAO*(ndx_i-1)
                        ki = ndx_k + iTmp
                        li = ndx_l + iTmp
                                     jTmp = nAO*(ndx_j-1)
                        kj = ndx_k + jTmp
                        lj = ndx_l + jTmp
                                     kTmp = nAO*(ndx_k-1)
                        ik = ndx_i + kTmp
                        jk = ndx_j + kTmp
                                     lTmp = nAO*(ndx_l-1)
                        il = ndx_i + lTmp
                        jl = ndx_j + lTmp

                        if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ik).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ik)
                           JL0(iCount) = iMap(jl)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(jl).gt.iBchOff(irp_ik)).and.
     &                      (iMap(jl).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(jl)
                           JL0(iCount) = iMap(ik)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(ki).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ki).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ki)
                           JL0(iCount) = iMap(lj)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(lj).gt.iBchOff(irp_ik)).and.
     &                      (iMap(lj).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(lj)
                           JL0(iCount) = iMap(ki)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(il).gt.iBchOff(irp_il)).and.
     &                      (iMap(il).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(il)
                           JL0(iCount) = iMap(jk)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(jk).gt.iBchOff(irp_il)).and.
     &                      (iMap(jk).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(jk)
                           JL0(iCount) = iMap(il)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(li).gt.iBchOff(irp_il)).and.
     &                      (iMap(li).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(li)
                           JL0(iCount) = iMap(kj)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(kj).gt.iBchOff(irp_il)).and.
     &                      (iMap(kj).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(kj)
                           JL0(iCount) = iMap(li)
                           dVal(iCount) = dBuf(int)
                        end if

c                    else if (ij.eq.kl) then
                     else

                        irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))
                        irp_il = dirprd(iAOSym(ndx_i),iAOSym(ndx_l))

                                     kTmp = nAO*(ndx_k-1)
                        ik = ndx_i + kTmp
                        jk = ndx_j + kTmp
                                     lTmp = nAO*(ndx_l-1)
                        il = ndx_i + lTmp
                        jl = ndx_j + lTmp

                        if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ik).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ik)
                           JL0(iCount) = iMap(jl)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(jl).gt.iBchOff(irp_ik)).and.
     &                      (iMap(jl).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(jl)
                           JL0(iCount) = iMap(ik)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(il).gt.iBchOff(irp_il)).and.
     &                      (iMap(il).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(il)
                           JL0(iCount) = iMap(jk)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(jk).gt.iBchOff(irp_il)).and.
     &                      (iMap(jk).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(jk)
                           JL0(iCount) = iMap(il)
                           dVal(iCount) = dBuf(int)
                        end if

c                    end if (ij.ne.kl)
                     end if

c                 else if (ndx_k.eq.ndx_l) then
                  else

                     irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))

                     ki = ndx_k + nAO*(ndx_i-1)
                     lj = ndx_l + nAO*(ndx_j-1)
                     ik = ndx_i + nAO*(ndx_k-1)
                     jl = ndx_j + nAO*(ndx_l-1)

                     if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ik).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ik)
                        JL0(iCount) = iMap(jl)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(jl).gt.iBchOff(irp_ik)).and.
     &                   (iMap(jl).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(jl)
                        JL0(iCount) = iMap(ik)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(ki).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ki).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ki)
                        JL0(iCount) = iMap(lj)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(lj).gt.iBchOff(irp_ik)).and.
     &                   (iMap(lj).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(lj)
                        JL0(iCount) = iMap(ki)
                        dVal(iCount) = dBuf(int)
                     end if

c                 end if (ndx_k.ne.ndx_l)
                  end if

c              else if (ndx_i.eq.ndx_j) then
               else

                  if (ndx_k.ne.ndx_l) then

                     irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))
                     irp_il = dirprd(iAOSym(ndx_i),iAOSym(ndx_l))

                                  iTmp = nAO*(ndx_i-1)
                     ki = ndx_k + iTmp
                     li = ndx_l + iTmp
                                  jTmp = nAO*(ndx_j-1)
                     kj = ndx_k + jTmp
                     lj = ndx_l + jTmp
                                  kTmp = nAO*(ndx_k-1)
                     ik = ndx_i + kTmp
                     jk = ndx_j + kTmp
                                  lTmp = nAO*(ndx_l-1)
                     il = ndx_i + lTmp
                     jl = ndx_j + lTmp

                     if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ik).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ik)
                        JL0(iCount) = iMap(jl)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(ki).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ki).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ki)
                        JL0(iCount) = iMap(lj)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(il).gt.iBchOff(irp_il)).and.
     &                   (iMap(il).le.iBchOff(irp_il)+
     &                                iBchLen(irp_il) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_il
                        IK0(iCount) = iMap(il)
                        JL0(iCount) = iMap(jk)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(li).gt.iBchOff(irp_il)).and.
     &                   (iMap(li).le.iBchOff(irp_il)+
     &                                iBchLen(irp_il) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_il
                        IK0(iCount) = iMap(li)
                        JL0(iCount) = iMap(kj)
                        dVal(iCount) = dBuf(int)
                     end if

c                 else if (ndx_k.eq.ndx_l) then
                  else

                     if (ij.eq.kl) then

                        irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))

                        ik = ndx_i + nAO*(ndx_k-1)
                        jl = ndx_j + nAO*(ndx_l-1)

                        if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ik).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ik)
                           JL0(iCount) = iMap(jl)
                           dVal(iCount) = dBuf(int)
                        end if

c                    else if (ij.ne.kl) then
                     else

                        irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))

                        ik = ndx_i + nAO*(ndx_k-1)
                        ki = ndx_k + nAO*(ndx_i-1)
                        jl = ndx_j + nAO*(ndx_l-1)
                        lj = ndx_l + nAO*(ndx_j-1)

                        if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ik).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ik)
                           JL0(iCount) = iMap(jl)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(ki).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ki).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ki)
                           JL0(iCount) = iMap(lj)
                           dVal(iCount) = dBuf(int)
                        end if

c                    end if (ij.eq.kl)
                     end if

c                 end if (ndx_k.ne.ndx_l)
                  end if

c              end if (ndx_i.ne.ndx_j)
               end if

c           end do int = 1, NUT
            end do

c            do irp_ij = 1, nirrep
c               irp_ab = dirprd(irp_ij,irp_x)
c               call wheneq(iCount,iSym,1,irp_ab,iWhere,nMatch)
c               do iSpin = 1, 3
c                  IKOff = ISTARTMO(irp_ij,iSpin)
c                  LenIK = irpdpd(irp_ij,isytyp(2,43+iSpin))
c                  do int = 1, nMatch
c                     ndx = iWhere(int)
c                     dX = dVal(ndx)
c                     ik = IK0(ndx)-1
c                     jl = JL0(ndx)-1
c                     iOff1 = LenIK*jl + IKOff
c                     iOff2 = LenIK*ik + IKOff
c                     do i = 1, LenIK
c                        dZ(i+iOff2) = dZ(i+iOff2) + dX*T2(i+iOff1)
c                     end do
cc                 end do int = 1, nMatch
c                  end do
cc              end do iSpin = 1, 3
c               end do
cc           end do irp_ij = 1, nirrep
c            end do

            do iSpin = 1, 3
               do ndx = 1, iCount
                  irp_ab = iSym(ndx)
                  irp_ij = dirprd(irp_ab,irp_x)
                  IKOff = ISTARTMO(irp_ij,iSpin)
                  LenIK = irpdpd(irp_ij,isytyp(2,43+iSpin))
                  iOff1 = LenIK*(JL0(ndx)-1) + IKOff
                  iOff2 = LenIK*(IK0(ndx)-1) + IKOff
                  dX = dVal(ndx)
                  do i = 1, LenIK
                     dZ(i+iOff2) = dZ(i+iOff2) + dX*T2(i+iOff1)
                  end do
               end do
            end do

c        end do while (NUT.ne.-1)
         end do

c ----------------------------------------------------------------------

c     else if (iUHF.eq.0) then
      else

         NUT = 0
         do while (NUT.ne.-1)

            read(iUnit) dBuf, iBuf, NUT

            iCount = 0
            do int = 1, NUT

               ndx_i = IUPKI(iBuf(int))
               ndx_j = IUPKJ(iBuf(int))
               ndx_k = IUPKK(iBuf(int))
               ndx_l = IUPKL(iBuf(int))

               ij = INDX(max(ndx_j,ndx_i),min(ndx_j,ndx_i))
               kl = INDX(max(ndx_l,ndx_k),min(ndx_l,ndx_k))

               if (ndx_i.ne.ndx_j) then

                  if (ndx_k.ne.ndx_l) then

                     if (ij.ne.kl) then

                        irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))
                        irp_il = dirprd(iAOSym(ndx_i),iAOSym(ndx_l))

                                     kTmp = nAO*(ndx_k-1)
                        ik = ndx_i + kTmp
                        jk = ndx_j + kTmp
                                     lTmp = nAO*(ndx_l-1)
                        il = ndx_i + lTmp
                        jl = ndx_j + lTmp

                        if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ik).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ik)
                           JL0(iCount) = iMap(jl)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(jl).gt.iBchOff(irp_ik)).and.
     &                      (iMap(jl).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(jl)
                           JL0(iCount) = iMap(ik)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(il).gt.iBchOff(irp_il)).and.
     &                      (iMap(il).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(il)
                           JL0(iCount) = iMap(jk)
                           dVal(iCount) = dBuf(int)
                        end if

                        if ((iMap(jk).gt.iBchOff(irp_il)).and.
     &                      (iMap(jk).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(jk)
                           JL0(iCount) = iMap(il)
                           dVal(iCount) = dBuf(int)
                        end if

c                    else if (ij.eq.kl) then
                     else

                        irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))
                        irp_il = dirprd(iAOSym(ndx_i),iAOSym(ndx_l))

                                     kTmp = nAO*(ndx_k-1)
                        ik = ndx_i + kTmp
                        jk = ndx_j + kTmp
                                     lTmp = nAO*(ndx_l-1)
                        il = ndx_i + lTmp
                        jl = ndx_j + lTmp

                        if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                      (iMap(ik).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(ik)
                           JL0(iCount) = iMap(jl)
                           dVal(iCount) = dBuf(int)*0.5d0
                        end if

                        if ((iMap(jl).gt.iBchOff(irp_ik)).and.
     &                      (iMap(jl).le.iBchOff(irp_ik)+
     &                                   iBchLen(irp_ik) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_ik
                           IK0(iCount) = iMap(jl)
                           JL0(iCount) = iMap(ik)
                           dVal(iCount) = dBuf(int)*0.5d0
                        end if

                        if ((iMap(il).gt.iBchOff(irp_il)).and.
     &                      (iMap(il).le.iBchOff(irp_il)+
     &                                   iBchLen(irp_il) )    ) then
                           iCount = iCount + 1
                           iSym(iCount) = irp_il
                           IK0(iCount) = iMap(il)
                           JL0(iCount) = iMap(jk)
                           dVal(iCount) = dBuf(int)
                        end if

c                    end if (ij.ne.kl)
                     end if

c                 else if (ndx_k.eq.ndx_l) then
                  else

                     irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))

                     ik = ndx_i + nAO*(ndx_k-1)
                     jl = ndx_j + nAO*(ndx_l-1)

                     if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ik).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ik)
                        JL0(iCount) = iMap(jl)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(jl).gt.iBchOff(irp_ik)).and.
     &                   (iMap(jl).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(jl)
                        JL0(iCount) = iMap(ik)
                        dVal(iCount) = dBuf(int)
                     end if

c                 end if (ndx_k.ne.ndx_l)
                  end if

c              else if (ndx_i.eq.ndx_j) then
               else

                  if (ndx_k.ne.ndx_l) then

                     irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))
                     irp_il = dirprd(iAOSym(ndx_i),iAOSym(ndx_l))

                                  kTmp = nAO*(ndx_k-1)
                     ik = ndx_i + kTmp
                     jk = ndx_j + kTmp
                                  lTmp = nAO*(ndx_l-1)
                     jl = ndx_j + lTmp
                     il = ndx_i + lTmp

                     if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ik).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ik)
                        JL0(iCount) = iMap(jl)
                        dVal(iCount) = dBuf(int)
                     end if

                     if ((iMap(il).gt.iBchOff(irp_il)).and.
     &                   (iMap(il).le.iBchOff(irp_il)+
     &                                iBchLen(irp_il) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_il
                        IK0(iCount) = iMap(il)
                        JL0(iCount) = iMap(jk)
                        dVal(iCount) = dBuf(int)
                     end if

c                 else if (ndx_k.eq.ndx_l) then
                  else

                     irp_ik = dirprd(iAOSym(ndx_i),iAOSym(ndx_k))

                     ik = ndx_i + nAO*(ndx_k-1)
                     jl = ndx_j + nAO*(ndx_l-1)

                     if ((iMap(ik).gt.iBchOff(irp_ik)).and.
     &                   (iMap(ik).le.iBchOff(irp_ik)+
     &                                iBchLen(irp_ik) )    ) then
                        iCount = iCount + 1
                        iSym(iCount) = irp_ik
                        IK0(iCount) = iMap(ik)
                        JL0(iCount) = iMap(jl)
                        if (ij.eq.kl) then
                           dVal(iCount) = dBuf(int) * 0.5d0
                        else
                           dVal(iCount) = dBuf(int)
                        end if
                     end if

c                 end if (ndx_k.ne.ndx_l)
                  end if

c              end if (ndx_i.ne.ndx_j)
               end if

c           end do int = 1, NUT
            end do

c            do irp_ij = 1, nirrep
c               irp_ab = dirprd(irp_ij,irp_x)
c               call wheneq(iCount,iSym,1,irp_ab,iWhere,nMatch)
c               IKOff = ISTARTMO(irp_ij,3)
c               LenIK = irpdpd(irp_ij,isytyp(2,46))
c               do int = 1, nMatch
c                  ndx = iWhere(int)
c                  iOff1 = IKOff + LenIK*(JL0(ndx)-1)
c                  iOff2 = IKOff + LenIK*(IK0(ndx)-1)
c                  dX = dVal(ndx)
c                  do i = 1, LenIK
c                     dZ(i+iOff2) = dZ(i+iOff2) + dX*T2(i+iOff1)
c                  end do
c               end do
c            end do

            do ndx = 1, iCount
               irp_ab = iSym(ndx)
               irp_ij = dirprd(irp_ab,irp_x)
               IKOff = ISTARTMO(irp_ij,3)
               LenIK = irpdpd(irp_ij,isytyp(2,46))
               iOff1 = IKOff + LenIK*(JL0(ndx)-1)
               iOff2 = IKOff + LenIK*(IK0(ndx)-1)
               dX = dVal(ndx)
               do i = 1, LenIK
                  dZ(i+iOff2) = dZ(i+iOff2) + dX*T2(i+iOff1)
               end do
            end do

c        end do while (NUT.ne.-1)
         end do

c ----------------------------------------------------------------------

c     end if (iUHF.ne.0)
      end if

      return
c     end subroutine rdaoijkl2
      end

