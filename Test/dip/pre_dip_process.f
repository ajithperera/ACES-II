





































































































































































































































































































































































































































































































































      Program pre_dip_process
       
      ImpliciT Double Precision (A-H, O-Z)

      Logical act_spc_cc
      Character*4 act
      Dimension Nocc(16)

c COMMON BLOCKS







c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end























c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end





c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end





c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end





c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end





c files.com : begin
      integer        luout, moints
      common /files/ luout, moints
c files.com : end



      Luout = 6
C
      call aces_init(icore,i0,icrsiz,iuhf,.True.)
      act_spc_cc = .true.
      act = "READ"
      if (act_spc_cc .and. act .eq."READ")
     &    call modf_amps(ICORE(I0),Icrsiz,IUHF,0,.TRUE.,'T',act,"AMPS")

C
c Read in the new occ numbers or other orbital manipulations from the
C GUESS file. 
c
      Write(6,*) "I am Here"
      Call Readgs(Nocc, Nirrep, Iuhf)
      Call Putrec(20,'JOBARC','OCCUPYA0',Nirrep,NOCC(1))
      If (Iuhf .gt. 0) Call Putrec(20,'JOBARC','OCCUPYB0',Nirrep,
     &                             NOCC(9))
C
C Set the scf_maxcyc to zero and rerun the SCF to built n+2 Fock matrix
C
      iflags(16) = 0
      Call Putrec(20, "JOBARC", "IFLAGS  ", 100, Iflags)
      Call Runit("xvscf") 
C

      Call aces_fin
       
      End