c***********************************************************************
c	subroutine that allows the use of switches.                      *
c     it is called at the beginning of jacobian.f, residual.f & rates.f*
c     Any switches need to be programmed here by hand in fortran       *
c     In maple you may use a dummy variable (e.g. H1), with which you  *
c     associate a condition here (e.g. H1=1 if oversaturated, 0 if not)*
c                                                                      *
c     if no switches are used, then provide this subroutine as a dummy *
c     (empty) procudure.                                               *
c                                                                      *
c     CM, June 2002                                                    *
c***********************************************************************

	subroutine switches(j)

      include 'common_geo.inc'
      include 'common.inc'

	integer j
	real*8 omega

c	goto 99
	omega = sp(8,j)*sp(4,j)/keqcal
	
!	if (omega.gt.1.d0.or.j.gt.(nx/2)) then
	if (omega.gt.1.d0) then
		h1 = 0.d0
	else
		h1 = 1.d0
	end if	

c	omega_mn = sp(12,j)*sp(4,j)/K_mnco3   ! version shallow
c      omega_fe = sp(14,j)*sp(4,j)/K_feco3  ! version shallow	
c	omega_fs = sp(14,j)*sp(17,j)/(sp(6,j)*K_fes) ! version shallow
c	omega_caco3 = sp(5,j)*sp(11,j)/(10**(-1.0*pk23))
	omega_mn = sp(10,j)*sp(4,j)/K_mnco3
      omega_fe = sp(12,j)*sp(4,j)/K_feco3
	omega_fs = sp(12,j)*sp(15,j)/(sp(6,j)*K_fes)
	

	if (omega_mn.gt.1.d0) then
c		sw22 = 1.0 ! version shallow
		sw17 = 1.0
	else
c		sw22= 0.0 ! version shallow
		sw17= 0.0
	endif
	if (omega_fe.gt.1.d0) then
c		sw23 = 1.d0 ! version shallow
		sw18 = 1.0
	else
c		sw23 = 0.d0 ! version shallow
		sw18 = 0.0
	endif
	if (omega_fs.gt.1.d0) then
c		sw24 = 1.d0 ! version shallow
		sw19 = 1.0
	else
c		sw24 = 0.d0 ! version shallow
		sw19 = 0.0
	endif

c	if (omega_caco3.gt.1.0) then
c		sw23 = 1.0
c	else
c		sw23 = 0.0
c	endif
	
c	if (j.eq.1) write(*,*) omega_fs, sw22, sw23, sw24
c	pause	

	if (sp(1,j) .gt. kmo2) then
		ho2 = 1.d0
		hno3 = 0.0
		hmn4 = 0.0
		hfe3 = 0.0
		hso4 = 0.0
	else if (sp(2,j) .gt. kmno3) then
		ho2 = 0.d0
		hno3 = 1.d0
		hno3f2 = 1.d0
		hmn4 = 0.0
		hfe3 = 0.0
		hso4 = 0.0
c	else if (sp(13,j) .gt. kmno2) then ! version shallow
	else if (sp(11,j) .gt. kmno2) then
		ho2 = 0.d0
		hno3 = 1.d0
		hno3f2 = 0.d0
		hmn4 = 1.d0
		hmn4f2 = 1.d0
		hfe3 = 0.0
		hso4 = 0.0
c	else if (sp(15,j) .gt. kmfeoh3) then ! version shallow
	else if (sp(13,j) .gt. kmfeoh3) then
		ho2 = 0.d0
		hno3 = 1.d0
		hno3f2 = 0.d0
		hmn4 = 1.d0
		hmn4f2 = 0.d0
		hfe3 = 1.0 
		hfe3f2 = 1.d0
		hso4 = 0.0
c	else if (sp(11,j) .gt. kmso4) then ! version shallow
	else if (sp(9,j) .gt. kmso4) then
		ho2 = 0.d0
		hno3 = 1.d0
		hno3f2 = 0.d0
		hmn4 = 1.d0
		hmn4f2 = 0.d0
		hfe3 = 1.0 
		hfe3f2 = 0.d0
		hso4 = 1.0 
		hso4f2 = 1.d0
	else
c		pause 'ran out of TEAs'
c		write(*,*) 'ran out of TEAs', j, sp(1,j), sp(7,j)
		ho2 = 0.d0
		hno3 = 1.d0
		hno3f2 = 0.d0
		hmn4 = 1.d0
		hmn4f2 = 0.d0
		hfe3 = 1.0 
		hfe3f2 = 0.d0
		hso4 = 1.0
		hso4f2 = 0.d0
	endif


99	continue	
      return
	end



