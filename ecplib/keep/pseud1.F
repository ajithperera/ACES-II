     Subroutine Pseud1(xa, ya, za, xb, yb, zb, xc, yc, zc, expa,
    &                  expb, 

     Implicit Double Precision (A-H, O-Z)

     
     


C
C Coordinate CA and CB (Vec(C-A) and Vec(C-B))
C
     Xca = Xc - Xa
     Yca = Yc - Ya
     Zca = Zc - Za

     Xcb = Xc - Xb
     Ycb = Yc - Yb
     Zcb = Zc - Zb
C 
     Ca = Dsqrt(Xca**2 + Yca**2 + Zca**2)
     Cb = Dsqrt(Xcb**2 + Ycb**2 + Zcb**2)


