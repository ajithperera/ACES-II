
      SUBROUTINE INCENT_CH(VLIST,NOCMX,I,IATOM)
      IMPLICIT NONE
      CHARACTER*8 VLIST,ASYMB
      INTEGER NOCMX,I,IATOM,J
      DIMENSION VLIST(NOCMX,10),ASYMB(103)
      DATA (ASYMB(J),J = 1,103)
     1/'      H ', '      He', '      Li', '      Be', '      B ',
     2 '      C ', '      N ', '      O ', '      F ', '      Ne',
     3 '      Na', '      Mg', '      Al', '      Si', '      P ',
     4 '      S ', '      Cl', '      Ar', '      K ', '      Ca',
     5 '      Sc', '      Ti', '      V ', '      Cr', '      Mn',
     6 '      Fe', '      Co', '      Ni', '      Cu', '      Zn',
     7 '      Ga', '      Ge', '      As', '      Se', '      Br',
     8 '      Kr', '      Rb', '      Sr', '      Y ', '      Zr',
     9 '      Nb', '      Mo', '      Tc', '      Ru', '      Rh',
     O '      Pd', '      Ag', '      Cd', '      In', '      Sn',
     1 '      Sb', '      Te', '      I ', '      Xe', '      Cs',
     2 '      Ba', '      La', '      Ce', '      Pr', '      Nd',
     3 '      Pm', '      Sm', '      Eu', '      Gd', '      Tb',
     4 '      Dy', '      Ho', '      Er', '      Tm', '      Yb',
     5 '      Lu', '      Hf', '      Ta', '      W ', '      Re',
     6 '      Os', '      Ir', '      Pt', '      Au', '      Hg',
     7 '      Tl', '      Pb', '      Bi', '      Po', '      At',
     8 '      Rn', '      Fr', '      Ra', '      Ac', '      Th',
     9 '      Pa', '      U ', '      Np', '      Pu', '      Am',
     O '      Cm', '      Bk', '      Cf', '      Es', '      Fm',
     1 '      Md', '      No', '      Lr' /
C
      IF(IATOM.GT.0)THEN
       VLIST(I, 5) = ASYMB(IATOM)
       VLIST(I, 9) = '        '
       VLIST(I,10) = '        '
      ELSE
       VLIST(I, 5) = '        '
       VLIST(I, 9) = '        '
       VLIST(I,10) = '        '
      ENDIF
C
      RETURN
      END
