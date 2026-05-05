c    $Id: main.f 16 2007-10-18 12:36:47Z centler $
      program main
c***********************************************************************
c  this code is the fortran side of the biogeochemical reaction        *
c  network simulator, which consist of a maple interface and a fortran *
c  backbone for problem setup and solving, respectively                *
c  please do not alter the source code without talking to the authors  *
c  (see header below)                                                  *
c                                                                      *
c  credits for program parts taken or adapted from the literature are  *
c  given in the individual subroutines                                 *
c                                                                      *
c  July 2002, CM                                                       *
c***********************************************************************
      include 'common_geo.inc'
      include 'common.inc'
      include 'common_drive.inc'
c***********************************************************************
c  some explanations for variables defined in common blocks            *
c  - common_geo.inc                                                    *
c  nsolid:     number of solid species                                 *
c  ndiss:          number of dissolved species                         *
c  ncomp:          total number of species (nsolid+ndiss)              *
c  nreac:         total number of reactions                            *
c  nx:        gridcell number (note that concentrations are calculated *
c             at every second point, mixing parameters at the nodes    *
c             in between. node 1 and nx are the boundary points        *
c  kinetics:     kinetic rate constants and equilibrium constants      *
c  physics:     constants describing the physical environment          *
c             plus discretization                                      *
c  timestuff:     end time of simulation. value is defined in basic.f  *
c  - common.inc                                                        *
c  tspt:          sp, co and spold contain concentration profiles      *
c             disp, dsol_0, f_T deal with diffusion/dispersion coeff.  *
c               r are the reaction rates                               *
c  bound:     defines type of boundary conditions (value given in spb) *
c  - common_opt.inc                                                    *
c  nopt:          number of optimized parameters                       *
c  ntopt:          number of times measurements were made              *
c  ntotparam:     total number of parameters (see transferfw/back)     *
c  par:          entire list of parameters                             *
c  timemeas:     times measurements were taken                         *
c  idpar:          identifcation number of the optimized parameters    *
c  - common_meas.inc                                                   *
c  maxxmeas:     maximum numbers of measurements in a single profile   *
c  maxspmeas:     maximum number of species measured at the same time  *
c  idspmeas:     identification number of the species measured         *
c  nrxmeas:     number of measurements in a profile                    *
c  nrspmeas:     number of species measured at a given time            *
c  spmeas:     measured concentration profile                          *
c  xmeas:          depths of measurements                              *
c  sigmeas:     standard deviations of the measurements                *
c  - common_drive.inc                                                  *
c  see drivervalues.f                                                  *
c***********************************************************************

!!c      elapsed_time = TIMEF()

       write(*,*)
      write(*,*) '_________________________BRNS________________________'
      write(*,*) '               reactive transport model              '
      write(*,*) '                   Version 2.0 - UFZ                 '
      write(*,*) '_____________________________________________________'
      write(*,*) '  contributers:                                      '
      write(*,*) '                                                     '
      write(*,*) '  P.Regnier and co-workers: basic RTM & concept      '
      write(*,*) '                                                     '
      write(*,*) '  Florian Centler, Martin Thullner:                  '
      write(*,*) '                coupling and generalization concept  '
      write(*,*) '                                                     '
      write(*,*) '  Florian Centler: Maple 10+ and DLL-versions        '
      write(*,*) '                                                     '
      write(*,*) '  All Rights Reserved                                '
      write(*,*) '_____________________________________________________'
      write(*,*)

       call printSvnVersion()

c***********************************************************************
c  DEFINITION OF DRIVER VARIABLES                                      *
c      defines the switches internal to fortran                        *
c***********************************************************************
       call drivervalues()

c***********************************************************************
c  START OF INITIALIZATION                                             *
c     sp        = concentration array  (derived from Maple 'variables')*
c     spb       = array of boundary cond. derived from Maple 'bnddata' *
c     basic     = routine defining Physical Parameters (Maple)         *
c     molecular     = routine defining the molecular diffusion coeffs. *
c                  and temperature dependance (Maple)                  *
c     boundaries = routine prescribing upper b.c.+initial cond.(Maple) *
c     biogeo     = routine defining the Biogeochemical param. (Maple)  *
c     initialcond= routine computing the initial profiles              *
c                  see routine drivervalues.f for the different options*
c                                                                      *
c     initial calls to get parameters. This is required if any of them *
c     are being optimized. if not this might be a waste of time ;-)    *
c     depth dependent profiles (porosity, area, D, v) are calculated   *
c     in diagenesis.f, because they may depend of the parameters to be *
c     optimized (e.g. dispersivity aL)                                 *
c***********************************************************************
c  obtain depth dependent parameters                                   *
c  this could be done in the main program only if porosity, area,      *
c  diffusion and advection velocities are not part of the optimization *
c  because I don't want to exclude this it is repeated here            *
c***********************************************************************
      call basic()
      call biogeo()
      call molecular()
     
      call gridsetup()
      call porarea()
      call advdiffcoeff()
      call printdepth()

      call boundaries()
      call initialcond()

c***********************************************************************
c  FINAL RUN                                                           *
c       diagenesis: forward model (Maple/Fortran)                      *
c       tstart: start time of simulation                               *
c       tend:     end time of simulation                               *
c***********************************************************************
      tstart = 0.d0
      tend = endt
       ntopt2 = 1

c***********************************************************************
c  PRINT FINAL RUN                                                     *
c   limited to 100 components in format statement                      *
c***********************************************************************
       open(unit=5,file='conc.dat',status='unknown')
       open(unit=4,file='conc.txt',status='unknown')
! 3001     format(1x,f8.4, 1x, 200(1x,e14.7))
 3001     format(200(1x,e14.7))

      do ii=1,ntopt2
         if ((tend.eq.0.).and.(nsstate.ne.1)) write(*,*) 'tend=0!'
        if (nsstate.eq.1) write(5,*) 'steady state, C(x), R(x)'
        if (nsstate.eq.2) write(5,*) 'steady state &transient C(x),R(x)'
        if ((nsstate.ne.1).and.(nsstate.ne.2)) then
           write(5,*) 'conc & rates @ time: ', tend
        end if
c       if (nsstate.eq.1) then
c          write(5,*) 'steady state concentration profile & rates'
c       else
c          write(5,*) 'conc & rates @ time: ', tend
c       end if
         call diagenesis(tstart,tend)
        do j=1,nx,2
           call rates(j)
           if(nsstate.eq.1) call rates(j)
          write(5,3001) x(j), (sp(i,j), i=1,ncomp),(r(k,j),k=1,nreac)
           ! create input file that can be read in initialcond.f
           if ((nsstate.eq.1).or.(nsstate.eq.2).or.(ii.eq.ntopt2)) then
c          if ((nsstate.eq.1).or.(ii.eq.ntopt2)) then
            write(4,3001) x(j), (sp(i,j), i=1,ncomp)
           end if
        end do
      end do
      close (4)
      close (5)

!!c      elapsed_time = TIMEF()
!!c      write(*,*) 'elapsed time: ', elapsed_time

      stop
      end






