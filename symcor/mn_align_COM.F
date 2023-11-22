      subroutine mn_align_COM(natoms, Q, atmass)
c     
c     translate cartesian coordianates to COM system.
c
      implicit none
      double precision CMX, CMY, CMZ, MOLWT, CM(3)
      integer i, j, natoms, luout
      double precision q(3*natoms), atmass(natoms)
c     
      logical print
      print = .false.
      luout = 6
C     
      if (print) then
         Write(6,"(a, 5I8)"), "The number of atoms", Natoms
         Write(6,*)
         Print*, "The Cartesians before translations to CM"
         call output(Q, 1, 3, 1, natoms, 3, natoms, 1)
c     Write(6, "(3F10.5)"), (Q(I), I=1, 3*NATOMS)
         Write(6,*)
      endif
C     
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      MOLWT=0.D0
      DO I = 1,NATOMS
         CMX = ATMASS(I)*Q(3*I-2)+CMX
         CMY = ATMASS(I)*Q(3*I-1)+CMY
         CMZ = ATMASS(I)*Q(3*I)+CMZ
         MOLWT = MOLWT+ATMASS(I)
      ENDDO
      IF (MOLWT.LT.1.0D-8) THEN
         WRITE(LUOUT,*) '@SYMMETRY: No real atoms in Z-matrix.'
         CALL ERREX
      END IF
      CM(1) = CMX/MOLWT
      CM(2) = CMY/MOLWT
      CM(3) = CMZ/MOLWT
      DO  30 I = 1,NATOMS
         DO 301 J = 0,2
            Q(3*I-J) = Q(3*I-J)-CM(3-J)
 301     continue
 30   continue
c     
c     jdw 5/26/95
c     
      IF (print) THEN
         write(6,*)
         write(6,*) ' @symmetry-i, Coordinates after  COM shift '
         call output(Q, 1, 3, 1, natoms, 3, natoms, 1)
c     do 1010 i=1,natoms
c     write(6,'(3F20.12)') q(3*i-2),q(3*i-1),q(3*i)
c     1010       continue
         write(6,*)
      ENDIF
c     
      return
      end
