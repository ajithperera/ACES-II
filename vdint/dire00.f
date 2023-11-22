      SUBROUTINE DIRE00
      IMPLICIT DOUBLE PRECISION (F,S), LOGICAL (D)
      DOUBLE PRECISION
     *        SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     *        SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z
      LOGICAL TWOCEN, THRCEN, FOUCEN, DERONE, DERTWO
      COMMON /EXPCOM/ SIGN1X, SIGN1Y, SIGN1Z, SIGN2X, SIGN2Y, SIGN2Z,
     *                SIGN3X, SIGN3Y, SIGN3Z, SIGN4X, SIGN4Y, SIGN4Z,
     *                NCENT1, NCENT2, NCENT3, NCENT4,
     *                ISO1,   ISO2,   ISO3,   ISO4,
     *                DERONE, DERTWO, TWOCEN, THRCEN, FOUCEN,
     *                NINTYP, NCCINT
      COMMON /DEODIR/ DEXA00, DEXB00, DEXC00, DEXD00,
     *                DEYA00, DEYB00, DEYC00, DEYD00,
     *                DEZA00, DEZB00, DEZC00, DEZD00,
     *                DEXXAA, DEXXAB, DEXXBB,
     *                DEXXCC, DEXXCD, DEXXDD,
     *                DEXYAA, DEXYAB, DEXYBA, DEXYBB,
     *                DEXYCC, DEXYCD, DEXYDC, DEXYDD,
     *                DEXZAA, DEXZAB, DEXZBA, DEXZBB,
     *                DEXZCC, DEXZCD, DEXZDC, DEXZDD,
     *                DEYYAA, DEYYAB, DEYYBB,
     *                DEYYCC, DEYYCD, DEYYDD,
     *                DEYZAA, DEYZAB, DEYZBA, DEYZBB,
     *                DEYZCC, DEYZCD, DEYZDC, DEYZDD,
     *                DEZZAA, DEZZAB, DEZZBB,
     *                DEZZCC, DEZZCD, DEZZDD
      COMMON /DETDIR/ DEXXAC, DEXYAC, DEXZAC,
     *                DEXYCA, DEYYAC, DEYZAC,
     *                DEXZCA, DEYZCA, DEZZAC,
     *                DEXXAD, DEXYAD, DEXZAD,
     *                DEXYDA, DEYYAD, DEYZAD,
     *                DEXZDA, DEYZDA, DEZZAD,
     *                DEXXBC, DEXYBC, DEXZBC,
     *                DEXYCB, DEYYBC, DEYZBC,
     *                DEXZCB, DEYZCB, DEZZBC,
     *                DEXXBD, DEXYBD, DEXZBD,
     *                DEXYDB, DEYYBD, DEYZBD,
     *                DEXZDB, DEYZDB, DEZZBD
      COMMON /DHODIR/ DHXP00, DHYP00, DHZP00,
     *                DHXQ00, DHYQ00, DHZQ00,
     *                DHXXPP, DHXYPP, DHXZPP, DHYYPP, DHYZPP, DHZZPP,
     *                DHXXQQ, DHXYQQ, DHXZQQ, DHYYQQ, DHYZQQ, DHZZQQ
      COMMON /DEOADR/ IEXA00, IEXB00, IEXC00, IEXD00,
     *                IEYA00, IEYB00, IEYC00, IEYD00,
     *                IEZA00, IEZB00, IEZC00, IEZD00,
     *                IEXXAA, IEXXAB, IEXXBB,
     *                IEXXCC, IEXXCD, IEXXDD,
     *                IEXYAA, IEXYAB, IEXYBA, IEXYBB,
     *                IEXYCC, IEXYCD, IEXYDC, IEXYDD,
     *                IEXZAA, IEXZAB, IEXZBA, IEXZBB,
     *                IEXZCC, IEXZCD, IEXZDC, IEXZDD,
     *                IEYYAA, IEYYAB, IEYYBB,
     *                IEYYCC, IEYYCD, IEYYDD,
     *                IEYZAA, IEYZAB, IEYZBA, IEYZBB,
     *                IEYZCC, IEYZCD, IEYZDC, IEYZDD,
     *                IEZZAA, IEZZAB, IEZZBB,
     *                IEZZCC, IEZZCD, IEZZDD
      COMMON /DETADR/ IEXXAC, IEXYAC, IEXZAC,
     *                IEXYCA, IEYYAC, IEYZAC,
     *                IEXZCA, IEYZCA, IEZZAC,
     *                IEXXAD, IEXYAD, IEXZAD,
     *                IEXYDA, IEYYAD, IEYZAD,
     *                IEXZDA, IEYZDA, IEZZAD,
     *                IEXXBC, IEXYBC, IEXZBC,
     *                IEXYCB, IEYYBC, IEYZBC,
     *                IEXZCB, IEYZCB, IEZZBC,
     *                IEXXBD, IEXYBD, IEXZBD,
     *                IEXYDB, IEYYBD, IEYZBD,
     *                IEXZDB, IEYZDB, IEZZBD
      COMMON /DHOADR/ IHXP00, IHYP00, IHZP00,
     *                IHXQ00, IHYQ00, IHZQ00,
     *                IHXXPP, IHXYPP, IHXZPP, IHYYPP, IHYZPP, IHZZPP,
     *                IHXXQQ, IHXYQQ, IHXZQQ, IHYYQQ, IHYZQQ, IHZZQQ
      COMMON /DEOFAC/ FEXA00, FEXB00, FEXC00, FEXD00,
     *                FEYA00, FEYB00, FEYC00, FEYD00,
     *                FEZA00, FEZB00, FEZC00, FEZD00,
     *                FEXXAA, FEXXAB, FEXXBB,
     *                FEXXCC, FEXXCD, FEXXDD,
     *                FEXYAA, FEXYAB, FEXYBA, FEXYBB,
     *                FEXYCC, FEXYCD, FEXYDC, FEXYDD,
     *                FEXZAA, FEXZAB, FEXZBA, FEXZBB,
     *                FEXZCC, FEXZCD, FEXZDC, FEXZDD,
     *                FEYYAA, FEYYAB, FEYYBB,
     *                FEYYCC, FEYYCD, FEYYDD,
     *                FEYZAA, FEYZAB, FEYZBA, FEYZBB,
     *                FEYZCC, FEYZCD, FEYZDC, FEYZDD,
     *                FEZZAA, FEZZAB, FEZZBB,
     *                FEZZCC, FEZZCD, FEZZDD
      COMMON /DETFAC/ FEXXAC, FEXYAC, FEXZAC,
     *                FEXYCA, FEYYAC, FEYZAC,
     *                FEXZCA, FEYZCA, FEZZAC,
     *                FEXXAD, FEXYAD, FEXZAD,
     *                FEXYDA, FEYYAD, FEYZAD,
     *                FEXZDA, FEYZDA, FEZZAD,
     *                FEXXBC, FEXYBC, FEXZBC,
     *                FEXYCB, FEYYBC, FEYZBC,
     *                FEXZCB, FEYZCB, FEZZBC,
     *                FEXXBD, FEXYBD, FEXZBD,
     *                FEXYDB, FEYYBD, FEYZBD,
     *                FEXZDB, FEYZDB, FEZZBD
      COMMON /DHOFAC/ FHXP00, FHYP00, FHZP00,
     *                FHXQ00, FHYQ00, FHZQ00,
     *                FHXXPP, FHXYPP, FHXZPP, FHYYPP, FHYZPP, FHZZPP,
     *                FHXXQQ, FHXYQQ, FHXZQQ, FHYYQQ, FHYZQQ, FHZZQQ
      ENTRY DIREA0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXA00 = .TRUE.
         DEYA00 = .TRUE.
         DEZA00 = .TRUE.
         IEXA00 = (INT1 - 1)*NCCINT
         IEYA00 = (INT2 - 1)*NCCINT
         IEZA00 = (INT3 - 1)*NCCINT
         FEXA00 = FAC1
         FEYA00 = FAC2
         FEZA00 = FAC3
      RETURN
      ENTRY DIREB0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXB00 = .TRUE.
         DEYB00 = .TRUE.
         DEZB00 = .TRUE.
         IEXB00 = (INT1 - 1)*NCCINT
         IEYB00 = (INT2 - 1)*NCCINT
         IEZB00 = (INT3 - 1)*NCCINT
         FEXB00 = FAC1
         FEYB00 = FAC2
         FEZB00 = FAC3
      RETURN
      ENTRY DIREC0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXC00 = .TRUE.
         DEYC00 = .TRUE.
         DEZC00 = .TRUE.
         IEXC00 = (INT1 - 1)*NCCINT
         IEYC00 = (INT2 - 1)*NCCINT
         IEZC00 = (INT3 - 1)*NCCINT
         FEXC00 = FAC1
         FEYC00 = FAC2
         FEZC00 = FAC3
      RETURN
      ENTRY DIRED0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DEXD00 = .TRUE.
         DEYD00 = .TRUE.
         DEZD00 = .TRUE.
         IEXD00 = (INT1 - 1)*NCCINT
         IEYD00 = (INT2 - 1)*NCCINT
         IEZD00 = (INT3 - 1)*NCCINT
         FEXD00 = FAC1
         FEYD00 = FAC2
         FEZD00 = FAC3
      RETURN
      ENTRY DIREP0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DHXP00 = .TRUE.
         DHYP00 = .TRUE.
         DHZP00 = .TRUE.
         IHXP00 = (INT1 - 1)*NCCINT
         IHYP00 = (INT2 - 1)*NCCINT
         IHZP00 = (INT3 - 1)*NCCINT
         FHXP00 = FAC1
         FHYP00 = FAC2
         FHZP00 = FAC3
      RETURN
      ENTRY DIREQ0(INT1,INT2,INT3,FAC1,FAC2,FAC3)
         DHXQ00 = .TRUE.
         DHYQ00 = .TRUE.
         DHZQ00 = .TRUE.
         IHXQ00 = (INT1 - 1)*NCCINT
         IHYQ00 = (INT2 - 1)*NCCINT
         IHZQ00 = (INT3 - 1)*NCCINT
         FHXQ00 = FAC1
         FHYQ00 = FAC2
         FHZQ00 = FAC3
      RETURN
      ENTRY DIREAA(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXAA = .TRUE.
         DEXYAA = .TRUE.
         DEXZAA = .TRUE.
         DEYYAA = .TRUE.
         DEYZAA = .TRUE.
         DEZZAA = .TRUE.
         IEXXAA = (INT1 - 1)*NCCINT
         IEXYAA = (INT2 - 1)*NCCINT
         IEXZAA = (INT3 - 1)*NCCINT
         IEYYAA = (INT4 - 1)*NCCINT
         IEYZAA = (INT5 - 1)*NCCINT
         IEZZAA = (INT6 - 1)*NCCINT
         FEXXAA = FAC1
         FEXYAA = FAC2
         FEXZAA = FAC3
         FEYYAA = FAC4
         FEYZAA = FAC5
         FEZZAA = FAC6
      RETURN
      ENTRY DIREBB(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXBB = .TRUE.
         DEXYBB = .TRUE.
         DEXZBB = .TRUE.
         DEYYBB = .TRUE.
         DEYZBB = .TRUE.
         DEZZBB = .TRUE.
         IEXXBB = (INT1 - 1)*NCCINT
         IEXYBB = (INT2 - 1)*NCCINT
         IEXZBB = (INT3 - 1)*NCCINT
         IEYYBB = (INT4 - 1)*NCCINT
         IEYZBB = (INT5 - 1)*NCCINT
         IEZZBB = (INT6 - 1)*NCCINT
         FEXXBB = FAC1
         FEXYBB = FAC2
         FEXZBB = FAC3
         FEYYBB = FAC4
         FEYZBB = FAC5
         FEZZBB = FAC6
      RETURN
      ENTRY DIRECC(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXCC = .TRUE.
         DEXYCC = .TRUE.
         DEXZCC = .TRUE.
         DEYYCC = .TRUE.
         DEYZCC = .TRUE.
         DEZZCC = .TRUE.
         IEXXCC = (INT1 - 1)*NCCINT
         IEXYCC = (INT2 - 1)*NCCINT
         IEXZCC = (INT3 - 1)*NCCINT
         IEYYCC = (INT4 - 1)*NCCINT
         IEYZCC = (INT5 - 1)*NCCINT
         IEZZCC = (INT6 - 1)*NCCINT
         FEXXCC = FAC1
         FEXYCC = FAC2
         FEXZCC = FAC3
         FEYYCC = FAC4
         FEYZCC = FAC5
         FEZZCC = FAC6
      RETURN
      ENTRY DIREDD(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DEXXDD = .TRUE.
         DEXYDD = .TRUE.
         DEXZDD = .TRUE.
         DEYYDD = .TRUE.
         DEYZDD = .TRUE.
         DEZZDD = .TRUE.
         IEXXDD = (INT1 - 1)*NCCINT
         IEXYDD = (INT2 - 1)*NCCINT
         IEXZDD = (INT3 - 1)*NCCINT
         IEYYDD = (INT4 - 1)*NCCINT
         IEYZDD = (INT5 - 1)*NCCINT
         IEZZDD = (INT6 - 1)*NCCINT
         FEXXDD = FAC1
         FEXYDD = FAC2
         FEXZDD = FAC3
         FEYYDD = FAC4
         FEYZDD = FAC5
         FEZZDD = FAC6
      RETURN
      ENTRY DIREPP(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DHXXPP = .TRUE.
         DHXYPP = .TRUE.
         DHXZPP = .TRUE.
         DHYYPP = .TRUE.
         DHYZPP = .TRUE.
         DHZZPP = .TRUE.
         IHXXPP = (INT1 - 1)*NCCINT
         IHXYPP = (INT2 - 1)*NCCINT
         IHXZPP = (INT3 - 1)*NCCINT
         IHYYPP = (INT4 - 1)*NCCINT
         IHYZPP = (INT5 - 1)*NCCINT
         IHZZPP = (INT6 - 1)*NCCINT
         FHXXPP = FAC1
         FHXYPP = FAC2
         FHXZPP = FAC3
         FHYYPP = FAC4
         FHYZPP = FAC5
         FHZZPP = FAC6
      RETURN
      ENTRY DIREQQ(INT1,INT2,INT3,INT4,INT5,INT6,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6)
         DHXXQQ = .TRUE.
         DHXYQQ = .TRUE.
         DHXZQQ = .TRUE.
         DHYYQQ = .TRUE.
         DHYZQQ = .TRUE.
         DHZZQQ = .TRUE.
         IHXXQQ = (INT1 - 1)*NCCINT
         IHXYQQ = (INT2 - 1)*NCCINT
         IHXZQQ = (INT3 - 1)*NCCINT
         IHYYQQ = (INT4 - 1)*NCCINT
         IHYZQQ = (INT5 - 1)*NCCINT
         IHZZQQ = (INT6 - 1)*NCCINT
         FHXXQQ = FAC1
         FHXYQQ = FAC2
         FHXZQQ = FAC3
         FHYYQQ = FAC4
         FHYZQQ = FAC5
         FHZZQQ = FAC6
      RETURN
      ENTRY DIREAB(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXAB = .TRUE.
         DEXYAB = .TRUE.
         DEXZAB = .TRUE.
         DEXYBA = .TRUE.
         DEYYAB = .TRUE.
         DEYZAB = .TRUE.
         DEXZBA = .TRUE.
         DEYZBA = .TRUE.
         DEZZAB = .TRUE.
         IEXXAB = (INT1 - 1)*NCCINT
         IEXYAB = (INT2 - 1)*NCCINT
         IEXZAB = (INT3 - 1)*NCCINT
         IEXYBA = (INT4 - 1)*NCCINT
         IEYYAB = (INT5 - 1)*NCCINT
         IEYZAB = (INT6 - 1)*NCCINT
         IEXZBA = (INT7 - 1)*NCCINT
         IEYZBA = (INT8 - 1)*NCCINT
         IEZZAB = (INT9 - 1)*NCCINT
         FEXXAB = FAC1
         FEXYAB = FAC2
         FEXZAB = FAC3
         FEXYBA = FAC4
         FEYYAB = FAC5
         FEYZAB = FAC6
         FEXZBA = FAC7
         FEYZBA = FAC8
         FEZZAB = FAC9
      RETURN
      ENTRY DIREAC(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXAC = .TRUE.
         DEXYAC = .TRUE.
         DEXZAC = .TRUE.
         DEXYCA = .TRUE.
         DEYYAC = .TRUE.
         DEYZAC = .TRUE.
         DEXZCA = .TRUE.
         DEYZCA = .TRUE.
         DEZZAC = .TRUE.
         IEXXAC = (INT1 - 1)*NCCINT
         IEXYAC = (INT2 - 1)*NCCINT
         IEXZAC = (INT3 - 1)*NCCINT
         IEXYCA = (INT4 - 1)*NCCINT
         IEYYAC = (INT5 - 1)*NCCINT
         IEYZAC = (INT6 - 1)*NCCINT
         IEXZCA = (INT7 - 1)*NCCINT
         IEYZCA = (INT8 - 1)*NCCINT
         IEZZAC = (INT9 - 1)*NCCINT
         FEXXAC = FAC1
         FEXYAC = FAC2
         FEXZAC = FAC3
         FEXYCA = FAC4
         FEYYAC = FAC5
         FEYZAC = FAC6
         FEXZCA = FAC7
         FEYZCA = FAC8
         FEZZAC = FAC9
      RETURN
      ENTRY DIREAD(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXAD = .TRUE.
         DEXYAD = .TRUE.
         DEXZAD = .TRUE.
         DEXYDA = .TRUE.
         DEYYAD = .TRUE.
         DEYZAD = .TRUE.
         DEXZDA = .TRUE.
         DEYZDA = .TRUE.
         DEZZAD = .TRUE.
         IEXXAD = (INT1 - 1)*NCCINT
         IEXYAD = (INT2 - 1)*NCCINT
         IEXZAD = (INT3 - 1)*NCCINT
         IEXYDA = (INT4 - 1)*NCCINT
         IEYYAD = (INT5 - 1)*NCCINT
         IEYZAD = (INT6 - 1)*NCCINT
         IEXZDA = (INT7 - 1)*NCCINT
         IEYZDA = (INT8 - 1)*NCCINT
         IEZZAD = (INT9 - 1)*NCCINT
         FEXXAD = FAC1
         FEXYAD = FAC2
         FEXZAD = FAC3
         FEXYDA = FAC4
         FEYYAD = FAC5
         FEYZAD = FAC6
         FEXZDA = FAC7
         FEYZDA = FAC8
         FEZZAD = FAC9
      RETURN
      ENTRY DIREBC(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXBC = .TRUE.
         DEXYBC = .TRUE.
         DEXZBC = .TRUE.
         DEXYCB = .TRUE.
         DEYYBC = .TRUE.
         DEYZBC = .TRUE.
         DEXZCB = .TRUE.
         DEYZCB = .TRUE.
         DEZZBC = .TRUE.
         IEXXBC = (INT1 - 1)*NCCINT
         IEXYBC = (INT2 - 1)*NCCINT
         IEXZBC = (INT3 - 1)*NCCINT
         IEXYCB = (INT4 - 1)*NCCINT
         IEYYBC = (INT5 - 1)*NCCINT
         IEYZBC = (INT6 - 1)*NCCINT
         IEXZCB = (INT7 - 1)*NCCINT
         IEYZCB = (INT8 - 1)*NCCINT
         IEZZBC = (INT9 - 1)*NCCINT
         FEXXBC = FAC1
         FEXYBC = FAC2
         FEXZBC = FAC3
         FEXYCB = FAC4
         FEYYBC = FAC5
         FEYZBC = FAC6
         FEXZCB = FAC7
         FEYZCB = FAC8
         FEZZBC = FAC9
      RETURN
      ENTRY DIREBD(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXBD = .TRUE.
         DEXYBD = .TRUE.
         DEXZBD = .TRUE.
         DEXYDB = .TRUE.
         DEYYBD = .TRUE.
         DEYZBD = .TRUE.
         DEXZDB = .TRUE.
         DEYZDB = .TRUE.
         DEZZBD = .TRUE.
         IEXXBD = (INT1 - 1)*NCCINT
         IEXYBD = (INT2 - 1)*NCCINT
         IEXZBD = (INT3 - 1)*NCCINT
         IEXYDB = (INT4 - 1)*NCCINT
         IEYYBD = (INT5 - 1)*NCCINT
         IEYZBD = (INT6 - 1)*NCCINT
         IEXZDB = (INT7 - 1)*NCCINT
         IEYZDB = (INT8 - 1)*NCCINT
         IEZZBD = (INT9 - 1)*NCCINT
         FEXXBD = FAC1
         FEXYBD = FAC2
         FEXZBD = FAC3
         FEXYDB = FAC4
         FEYYBD = FAC5
         FEYZBD = FAC6
         FEXZDB = FAC7
         FEYZDB = FAC8
         FEZZBD = FAC9
      RETURN
      ENTRY DIRECD(INT1,INT2,INT3,INT4,INT5,INT6,INT7,INT8,INT9,
     *             FAC1,FAC2,FAC3,FAC4,FAC5,FAC6,FAC7,FAC8,FAC9)
         DEXXCD = .TRUE.
         DEXYCD = .TRUE.
         DEXZCD = .TRUE.
         DEXYDC = .TRUE.
         DEYYCD = .TRUE.
         DEYZCD = .TRUE.
         DEXZDC = .TRUE.
         DEYZDC = .TRUE.
         DEZZCD = .TRUE.
         IEXXCD = (INT1 - 1)*NCCINT
         IEXYCD = (INT2 - 1)*NCCINT
         IEXZCD = (INT3 - 1)*NCCINT
         IEXYDC = (INT4 - 1)*NCCINT
         IEYYCD = (INT5 - 1)*NCCINT
         IEYZCD = (INT6 - 1)*NCCINT
         IEXZDC = (INT7 - 1)*NCCINT
         IEYZDC = (INT8 - 1)*NCCINT
         IEZZCD = (INT9 - 1)*NCCINT
         FEXXCD = FAC1
         FEXYCD = FAC2
         FEXZCD = FAC3
         FEXYDC = FAC4
         FEYYCD = FAC5
         FEYZCD = FAC6
         FEXZDC = FAC7
         FEYZDC = FAC8
         FEZZCD = FAC9
      RETURN
      END
