      subroutine oct(atomchrg,atmvc,rij,
     &    radgrid,relsiz,rtmp,rint,radpt,rwt,
     &    bslrd,gridxyz,grdangpts)

c This routine sets up the numerical integration grid

#include <aces.h>
#include <io.par>
#include <const.par>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>
#include <int.com>
      integer
     &    atomchrg(natoms),radgrid(int_numradpts),
     &    grdangpts(numgrid)
      M_REAL
     &    atmvc(natoms,natoms,3),
     &    rij(natoms,natoms),
     &    relsiz(natoms,natoms),rtmp(natoms),
     &    rint(natoms,maxangpts),radpt(int_numradpts),
     &    rwt(int_numradpts,maxangpts),bslrd(natoms),
     &    gridxyz(3,maxangpts,numgrid)

      integer
     &    iradpt,iangpt,i,midpt,j,jcntr
      M_REAL
     &    xnradp1,rad,xiradpt,eps,
     &    xn,xinlngth,xinl,xmid,xi,zz,p1,p2,xj,p3,pp,zz1,
     &    rwttmp,xoulngth

      call callstack_push('OCT')

      eps=three*ten**(-14)
      xnradp1=dble(int_numradpts+1)

c Use Slater radii for radial integration scaling
      if(int_radscal.eq.int_radscal_slater)then
        rad=atmrad(atomchrg(icntr))

c all molecular centers use the same integration radius
      else
        rad=one
      endif

c If using rigid partitioning
      if(int_parttyp.eq.int_parttyp_rigid)then

        write(stdout,9000)
 9000   format(t3,'@OCT - rigid partitioning not yet supported.')
        call errex
        return

c$$$c Determine the radial limit of integration
c$$$
c$$$        rintmin=ten**10
c$$$        do 40 iangpt=1,maxangpts
c$$$          do 20 iatom=1,natoms
c$$$c           if(iatom.eq.icntr)goto 20
c$$$            rtmp(iatom)=zero
c$$$            do 10 i=1,3
c$$$              rtmp(iatom)=rtmp(iatom)+(atmvc(icntr,iatom,i)+
c$$$     &            rstep*gridxyz(i,iangpt,maxanggrd))**2
c$$$   10       continue
c$$$            rtmp(iatom)=sqrt(rtmp(iatom))
c$$$   20     continue
c$$$
c$$$          do 30 iatom=1,natoms
c$$$            if(iatom.eq.icntr)goto 30
c$$$            cosangl=(rtmp(icntr)**2 + rij(icntr,iatom)**2 -
c$$$     &          rtmp(iatom)**2)/
c$$$     &          (two*rtmp(icntr)*rij(icntr,iatom))
c$$$            if(cosangl.gt.1.d-10)then
c$$$              rint(iatom,iangpt)=
c$$$     &            relsiz(icntr,iatom)*rij(icntr,iatom)/cosangl
c$$$            else
c$$$              rint(iatom,iangpt)=ten**10
c$$$            endif
c$$$            if(rint(iatom,iangpt).lt.rintmin)
c$$$     &          rintmin=rint(iatom,iangpt)
c$$$   30     continue
c$$$   40   continue
c$$$
c$$$c Determine radial point and weight
c$$$
c$$$        iset=0
c$$$        radpt(1)=rad / ((xnradp1-one)**2)
c$$$        bndry=(two*rad*xnradp1) / ((xnradp1-one)**3)
c$$$        rwt(1,1)=bndry*radpt(1)**2
c$$$        do 140 iradpt=2,int_numradpts
c$$$          xiradpt=dble(iradpt)
c$$$          radpt(iradpt)=(rad*xiradpt**2) / ((xnradp1-xiradpt)**2)
c$$$          dis=(two*rad*xnradp1*xiradpt) / ((xnradp1-xiradpt)**3)
c$$$          bndry2=bndry+dis
c$$$
c$$$c The radial boundary is less than the shortest radial line,
c$$$c determine the radial weight in the normal fashion
c$$$          if(bndry2.lt.rintmin)then
c$$$            rwt(iradpt,1)=dis*radpt(iradpt)**2
c$$$            bndry=bndry2
c$$$
c$$$          else
c$$$
c$$$c The radial boundary is greater than the shortest radial line, set
c$$$c the angular grid for a smooth transition to the top angular grid
c$$$c which is used throughout the interstitial region
c$$$
c$$$            if(iset.eq.1)goto 110
c$$$            radgrid(iradpt-1)=itopgrd
c$$$            do 50 iangpt=2,maxangpts
c$$$              rwt(iradpt-1,iangpt)=rwt(iradpt-1,1)
c$$$   50       continue
c$$$            icnt=2
c$$$            do 80 ii=1,itopgrd-1
c$$$              if(radgrid(iradpt-icnt).lt.itopgrd-ii)
c$$$     &            radgrid(iradpt-icnt)=itopgrd-ii
c$$$              do 60 iangpt=2,grdangpts(radgrid(iradpt-icnt))
c$$$                rwt(iradpt-icnt,iangpt)=rwt(iradpt-icnt,1)
c$$$   60         continue
c$$$              icnt=icnt+1
c$$$              if(radgrid(iradpt-icnt).lt.itopgrd-ii)
c$$$     &            radgrid(iradpt-icnt)=itopgrd-ii
c$$$              do 70 iangpt=2,grdangpts(radgrid(iradpt-icnt))
c$$$                rwt(iradpt-icnt,iangpt)=rwt(iradpt-icnt,1)
c$$$   70         continue
c$$$              icnt=icnt+1
c$$$   80       continue
c$$$
c$$$c Fill out the angular weights for the rest of the radial points
c$$$            do 100 ii=1,iradpt-icnt
c$$$              do 90 iangpt=2,grdangpts(radgrid(ii))
c$$$                rwt(ii,iangpt)=rwt(ii,1)
c$$$   90         continue
c$$$  100       continue
c$$$
c$$$            iset=1
c$$$
c$$$  110       continue
c$$$
c$$$            radgrid(iradpt)=itopgrd
c$$$
c$$$c Radial lines will now end at different lengths for different angles
c$$$            do 130 iangpt=1,grdangpts(itopgrd)
c$$$
c$$$c Determine if the radial boundary is less than the current radial line
c$$$
c$$$              rinttmp=ten**10
c$$$              iset2=0
c$$$              do 120 iatom=1,natoms
c$$$                if(iatom.eq.icntr)goto 120
c$$$                if(rint(iatom,iangpt).lt.bndry2)then
c$$$                  if(rint(iatom,iangpt).lt.rinttmp)then
c$$$                    rinttmp=rint(iatom,iangpt)
c$$$                    iset2=1
c$$$                  endif
c$$$                endif
c$$$  120         continue
c$$$
c$$$c The radial boundary is less than the current radial line,
c$$$c determine the radial weight in the normal fashion
c$$$
c$$$              if(iset2.eq.0)then
c$$$                rwt(iradpt,iangpt)=dis*radpt(iradpt)**2
c$$$              else
c$$$
c$$$c The current radial line is greater than the radial boundary, adjust
c$$$c the radial weight accordingly
c$$$
c$$$                dis2=rinttmp-bndry
c$$$                if(dis2.lt.zero)then
c$$$                  rwt(iradpt,iangpt)=zero
c$$$                else
c$$$                  rwt(iradpt,iangpt)=dis2*radpt(iradpt)**2
c$$$                endif
c$$$              endif
c$$$  130       continue
c$$$
c$$$            bndry=bndry2
c$$$
c$$$          endif
c$$$
c$$$  140   continue

      else

c For fuzzy partitioning (or atoms)

c For Handy radial integration

        if (int_radtyp.eq.int_radtyp_handy) then
          do 160 iradpt=1,int_numradpts
            xiradpt=dble(iradpt)
            radpt(iradpt)=(rad*xiradpt**2) / ((xnradp1-xiradpt)**2)
            rwt(iradpt,1)=(two*rad**3*xnradp1*xiradpt**5) /
     &          ((xnradp1-xiradpt)**7)
            do 150 iangpt=2,grdangpts(radgrid(iradpt))
              rwt(iradpt,iangpt)=rwt(iradpt,1)
  150       continue
  160     continue

        else

c Determine inner radial integration (Gauss-Legendre)

          if(natoms.eq.1)then
            xinlngth=int_radlimit
            numradpts=int_numradpts
            goto 510
          endif

          numradpts=int_numradpts/2
          xinlngth=ten**10

          do 310 i=1,natoms
            if(icntr.eq.i) goto 310
            if(rij(icntr,i).lt.xinlngth) then
              jcntr=i
              xinlngth=rij(icntr,i)
            endif
  310     continue

          if(xinlngth.gt.three)xinlngth=two

c          xinlngth=half*xinlngth
          xinlngth=(bslrd(icntr)/(bslrd(icntr)+bslrd(jcntr)))
     &      *xinlngth

  510     continue

          xn=dble(numradpts)
          midpt=(numradpts+1)/2
          xmid=half*xinlngth
          xinl=half*xinlngth

          do 320 i=1,midpt
            xi=dble(i)
            zz=cos(pi*(xi-one/four)/(xn+half))

  330       continue

            p1=one
            p2=zero

            do 340 j=1,numradpts
              xj=dble(j)
              p3=p2
              p2=p1
              p1=((two*xj-one)*zz*p2-(xj-one)*p3)/xj
  340       continue

            pp=xn*(zz*p1-p2)/(zz*zz-one)
            zz1=zz
            zz=zz1-p1/pp
            if(abs(zz-zz1).gt.eps)go to 330
            radpt(i)=xmid-xinl*zz
            radpt(numradpts+1-i)=xmid+xinl*zz
            rwttmp=two*xinl/((one-zz*zz)*pp*pp)
            rwt(i,1)=(radpt(i)**2)*rwttmp
            rwt(numradpts+1-i,1)=(radpt(numradpts+1-i)**2)*rwttmp
            do 350 iangpt=2,grdangpts(radgrid(i))
              rwt(i,iangpt)=rwt(i,1)
  350       continue
            do 360 iangpt=2,grdangpts(radgrid(numradpts+1-i))
              rwt(numradpts+1-i,iangpt)=rwt(numradpts+1-i,1)
  360       continue
  320     continue

          if(natoms.eq.1) goto 520

c Determine outer radial integration (Gauss-Legendre)

          xoulngth=11*one

          xmid=half*(xinlngth+xoulngth)
          xinl=half*(xoulngth-xinlngth)

          do 420 i=1,midpt
            xi=dble(i)
            zz=cos(pi*(xi-one/four)/(xn+half))

  430       continue

            p1=one
            p2=zero

            do 440 j=1,numradpts
              xj=dble(j)
              p3=p2
              p2=p1
              p1=((two*xj-one)*zz*p2-(xj-one)*p3)/xj
  440       continue

            pp=xn*(zz*p1-p2)/(zz*zz-one)
            zz1=zz
            zz=zz1-p1/pp
            if(abs(zz-zz1).gt.eps)go to 430
            radpt(numradpts+i)=xmid-xinl*zz
            radpt(2*numradpts+1-i)=xmid+xinl*zz
            rwttmp=two*xinl/((one-zz*zz)*pp*pp)
            rwt(numradpts+i,1)=(radpt(numradpts+i)**2)*rwttmp
            rwt(2*numradpts+1-i,1)=(radpt(2*numradpts+1-i)**2)*
     &          rwttmp
            do 450 iangpt=2,grdangpts(radgrid(numradpts+i))
              rwt(numradpts+i,iangpt)=rwt(numradpts+i,1)
  450       continue
            do 460 iangpt=2,grdangpts(radgrid(2*numradpts+1-i))
              rwt(2*numradpts+1-i,iangpt)=rwt(2*numradpts+1-i,1)
  460       continue
  420     continue

  520     continue

c Determine outer radial integration (Handy method)
c          xinradp1=dble(numradpts+1)
c          do 370 iradpt=1,numradpts
c            xiradpt=dble(iradpt)
c            radpt(iradpt+numradpts)=(xiradpt**2)/((xinradp1-xiradpt)**2)
c     &          + xinlngth
c            rwt(iradpt+numradpts,1)=((two*xinradp1*xiradpt) /
c     &          ((xinradp1-xiradpt)**3))*radpt(iradpt+numradpts)**2
c            do 380 iangpt=2,grdangpts(radgrid(iradpt+numradpts))
c              rwt(iradpt+numradpts,iangpt)=rwt(iradpt+numradpts,1)
c 380        continue
c 370      continue

        endif

      endif

      call callstack_pop
      return
      end
