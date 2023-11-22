      subroutine chklist(iScr,iFam)
      implicit integer (a-z)
      double precision snrm2, x
      integer aces_list_rows, aces_list_cols
 100  format(T4,'SUBLIST',T14,'LIST',T30,' NORM ')
      write(*,100)
      x = 0.d0
      do iGrp = 1, 10
         nrows = aces_list_rows(iGrp,iFam)
         ncols = aces_list_cols(iGrp,iFam)
         if ((ncols.ne.0).and.(nrows.ne.0)) then
            call getlst(iScr,1,ncols,1,iGrp,iFam)
            x = x + snrm2(nrows*ncols,iScr,1)**2
         end if
      end do
      write(*,'(T7,I3,T26,F20.10)') iFam, x
      return
      end
