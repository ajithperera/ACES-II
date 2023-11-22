      subroutine oepenrscl (
     &     nbas, nocc, type,
     &     e, x)
c
c This routine scales the given intermediate with the energy
c denominators.
c
c Argument TYPE specifies what type of intermediate is to be scaled.
c It takes the following values: 
c   'PP  ' : X (a,b) = X (a, b) / (e_a - e_b) for a \neq b;
c   'PH  ' : X (a,i) = X (a, i) / (e_i - e_a);
c   'HH  ' : X (i,j) = X (i, j) / (e_i - e_j) for i \neq j;
c   'PPPH' : X (ij, ab) = X (ij, ab) / e_ij^ab, et cetera
c
c   'dePH' : X (a,i) = X (a,i) * (e_i - e_a) - "descale"
c
c In: nbas, nocc, e, x, type
c Out: x is modified
c
c Igor Schweigert, Sep 2003
c $Id: oepenrscl.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c 
      implicit none
c
c     Arguments
c 
      integer
     &     nbas, nocc
c 
      double precision
     &     e (*), x (*)
c 
      character*4
     &     type
c
c     Local variables
c 
      integer
     &     n
c
c     Define the one- and two-electron offset operators
c
#define DEFINE_OPERATORS
#include <oepintind.com>
#define N_BAS nbas
#define N_OCC nocc
#define DEFINE_VARIABLES
#include <oepintind.com>
c     
c     Messages to the user
c 
 1000 format(
     &     t3,'@OEPENRSCL-F: Argument TYPE has an illegal value ',
     &     a4,'.')
c
c     Check what intermediate to scale and scale it. So far only IA intermediates.
c     
      if (type  .eq. 'PP  ') then
         do n = 1, n_pp
            j = i2_pp (n)
            i = i1_pp (n, j)
            x (n) = x (n) / (e (i+n_h) - e (j+n_h))
         enddo
      elseif (type .eq. 'PH' ) then
         do n = 1, n_ph
            j = i2_ph (n)
            i = i1_ph (n, j)
            x (n) = x (n) / (e (j) - e (i+n_h))
         enddo
      elseif (type .eq. 'HH  ') then
         do n = 1, n_hh
            j = i2_hh (n)
            i = i1_hh (n, j)
            x (n)= x (n) / (e (i) - e (j))
         enddo
      elseif (type .eq. 'dePH') then
         do n = 1, n_ph
            j = i2_ph (n)
            i = i1_ph (n, j)
            x (n) = x (n) * (e (j) - e (i+n_h))
         enddo
      else
         write (6, 1000) type
         call errex
      endif
c 
      return
      end

         

