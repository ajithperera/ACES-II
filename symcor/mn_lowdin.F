      subroutine mn_lowdin(pvec, qvec, n, np, nq, scr, mxcor)
c     
c     input are two sets of non-orthogonal orbitals
c     pvec(n,np) and qvec(n,nq) with np+nq = n
c     
c     the pvec are lowdin orthogonalized.
c     the qvec are orthogonalized to pvec
c     the qvec are lowdin orthogonalized
c     
c     the input and output vectors are contained in pvec and qvec.
c     
      implicit none
c     
      integer n, np, nq, mxcor
      double precision
     $     scr(mxcor), pvec(n,np), qvec(n,np), sqrt, one, zero
      integer is, iscr, idum, ip, iq, 
     $     ieval, ievec, iexcp, ievalsel, ioff, i,
     $     length
      logical print
c     
      print = .false.
      one = 1.0d0
      zero = 0.0d0
c     
c     we call the eig routine from aces2 ..... See if this leads to troubles....
c     
      if (np .gt. 0) then
c     
c     calculate overlap matrix from pvec
c     
         is = 1
         call dgemm('T', 'N', np, np, n, one, pvec, n,
     $        pvec, n, zero, scr(is), np)
c     
c     form s^-1/2 in np*np basis
c     
         idum = is + np*np
         iscr = idum + np*np
         call eig(scr(is), scr(idum), np, np, 0)
c     
         ioff=is
         do  i=1,np
            scr(ioff)=one/sqrt(scr(ioff))
            ioff=ioff+np+1
         enddo
c     
c     now form S^(-1/2) matrix from eigenvalues
c     
         call dgemm('N','N',np,np,np,one,scr(idum),np,
     &        scr(is),np,zero,scr(iscr),np)
         call dgemm('N','T',np,np,np,one,scr(iscr),np,
     &        scr(idum),np,zero,scr(is),np)
c     
c     transform pvec to orthogonal representation
c     
         call dgemm('N','N',n,np,np,one,pvec,n,
     &        scr(is),np,zero,scr(iscr),n)
c     
c     copy to pvec. That is done
c     
         call DCOPY(n*np, scr(iscr), 1, pvec, 1)
c     
c     mn         write(6,*) ' @Lowdin: pvec:'
c     mn         call output(pvec, 1, n, 1, np, n, np, 1)
c     
c     check orthogonality
c     
         call dgemm('T', 'N', np, np, n, one, pvec, n,
     $        pvec, n, zero, scr(is), np)
c     
         if (print) then
            write(6,*) ' @Lowdin: orthogonality of pvec'
            call output(scr(is), 1, np, 1, np, np, np, 1)
         endif
         if (nq .gt. 0) then
c     
c     construct projector on pvec
c     
            ip = 1
            iq = ip + n*n
c     
            call dgemm('N', 'T', n, n, np, one, pvec, n,
     $           pvec, n, zero, scr(ip), n)
c     
c     project: qvec = qvec - P * qvec
c     
            call DCOPY(n*nq, qvec, 1, scr(iq), 1)
            call dgemm('N', 'N', n, nq, n, -one, scr(ip), n,
     $           scr(iq), n, one, qvec, n)
c     
            if (print) then
               write(6,*) ' @Lowdin: qvec after projection'
               call output(qvec, 1, n, 1, nq, n, nq, 1)
            endif
c     
            call dgemm('T', 'N', np, nq, n, one, pvec, n,
     $           qvec, n, zero, scr(ip), np)
c     
            if (print) then
               write(6,*) ' @Lowdin: orthogonality of pvec and qvec'
               call output(scr(ip), 1, np, 1, nq, np, nq, 1)
            endif

c     
c     now Lowdin orthogonalize qvec as before
c     
         endif
      endif
c     
      if (nq .gt. 0) then
c     
c     calculate overlap matrix from qvec
c     
         is = 1
         call dgemm('T', 'N', nq, nq, n, one, qvec, n,
     $        qvec, n, zero, scr(is), nq)
c     
c     mn         write(6,*) ' overlap matrix qvec'
c     mn         call output(scr(is), 1, nq, 1, nq, nq, nq, 1)
c     
c     form s^-1/2 in nq*nq basis
c     
         idum = is + nq*nq
         iscr = idum + nq*nq
         call eig(scr(is), scr(idum), nq, nq, 0)

c     
c     mn         write(6,*) ' eigenvectors S'
c     mn         call output(scr(idum), 1, nq, 1, nq, nq, nq, 1)
c     
c     check orthonormality eigenvectors ...
c     
         call dgemm('T', 'N', nq, nq, nq, one, scr(idum), nq,
     $        scr(idum), nq, zero, scr(iscr), nq)
c     
c     mn         write(6,*) ' orthogonality eigenvectors of S'
c     mn         call output(scr(iscr), 1, nq, 1, nq, nq, nq, 1)
c     
         ioff=is
         do  i=1,nq
            scr(ioff)=one/sqrt(scr(ioff))
            ioff=ioff+nq+1
         enddo
c     
c     write(6,*) '   We reform S, not S^(-1/2) !!'
c     
c     mn         write(6,*) ' eigenvalues of S^(-1/2) for qvec'
c     mn         call output(scr(is), 1, nq, 1, nq, nq, nq, 1)
c     mn         write(6,*) ' eigenvectors of S '
c     mn         call output(scr(idum), 1, nq, 1, nq, nq, nq, 1)
c     
c     now form S^(-1/2) matrix from eigenvalues
c     
         call dgemm('N','N',nq,nq,nq,one,scr(idum),nq,
     &        scr(is),nq,zero,scr(iscr),nq)
         call dgemm('N','T',nq,nq,nq,one,scr(iscr),nq,
     &        scr(idum),nq,zero,scr(is),nq)
c     
c     mn         write(6,*) ' reconstructed S '
c     mn         call output(scr(is), 1, nq, 1, nq, nq, nq, 1)
c     
c     transform qvec to orthogonal representation
c     
         call dgemm('N','N',n,nq,nq,one,qvec,n,
     &        scr(is),nq,zero,scr(iscr),n)
c     
c     copy to qvec.
c     
         call DCOPY(n*nq, scr(iscr), 1, qvec, 1)      
c     
         if (print) then
            write(6,*) ' @Lowdin: qvec after orthogonalization'
            call output(qvec, 1, n, 1, nq, n, nq, 1)
         endif
c     
c     
c     check orthogonality
c     
         call dgemm('T', 'N', nq, nq, n, one, qvec, n,
     $        qvec, n, zero, scr(iscr), nq)
c     
c     mn         write(6,*) ' @Lowdin: orthogonality of qvec'
c     mn         call output(scr(iscr), 1, nq, 1, nq, nq, nq, 1)
c     
         call dgemm('T', 'N', np, nq, n, one, pvec, n,
     $        qvec, n, zero, scr(iscr), np)
c     
         if (print) then
            write(6,*) ' @Lowdin: orthogonality of pvec and qvec'
            call output(scr(iscr), 1, np, 1, nq, np, nq, 1)
         endif

      endif
c     
c     We now have a complete set of orthonormal vectors contained in
c     pvec, qvec
c     
      return
      end
