c      
c     SUBROUTINE biogeo
c      
      subroutine biogeo()
        include 'common_geo.inc'
        include 'common.inc'
          kdet = 0.1D1
          kmax1 = 0.1D1
          kl = 0.5444503676D-2
          Yo = 0.25D0
          st = 0.1D0
          ksodoc = 0.5D3
          ksox = 0.2D2
          ksamm = 0.2D2
          doxmin = 0.1D-1*kmax1
          Bfmax = 0.5D0
          Yn = 0.17D0
          ksndoc = 0.5D3
          kmax2 = 0.18D0
          ksno3 = 0.1D3
          kindox = 0.15D2
          no3min = 0.18D-2
          katt = 0.3D1
          kdeac = 0.1D1
          kmax3 = 0.3D-1
          ksso4 = 0.5D2
          kssdoc = 0.1D4
          kinno3 = 0.5D2
          so4min = 0.3D-3
          Ys = 0.4D-1
          kmax4 = 0.1D0
          kreac = 0.3D0
          km = 0.1D-1
          kpd = 0.15D-4
          ammin = 0.1D-2
          swAquaDiva = 1
          sw2 = 1
          ksnit = 0.2D2
          fdoco = 1.D0
          fdocn = 5.D0
          fdocs = 2.D0
          foo = 1.D0
          foa = 1.D0
          fnitra = 4.D0
          fsulph = 1.D0
          fcn = 0.1D0
          faa = 0.5D0
          Ya = 0.76D-1
          amming = 0.1D1
          kmax5 = 0.1D-1
          Ym = 0.38D-2
          ws = 0.1D1
          vel = 0
          khyd = 0.6D-4
          ksoxa = 0.2D1
          ktrans = 0.1D2
          faz = 0.3D2
          kmax6 = 0.3D-1
      end
