c      
c     SUBROUTINE notransport
c      
      subroutine notransport(k,itransp)
        include 'common_geo.inc'
        include 'common.inc'
        integer k,itransp
            if (k.eq.1) then
              itransp = 1
            endif
            if (k.eq.6) then
              itransp = 1
            endif
            if (k.eq.8) then
              itransp = 1
            endif
            if (k.eq.11) then
              itransp = 1
            endif
            if (k.eq.13) then
              itransp = 1
            endif
            if (k.eq.16) then
              itransp = 1
            endif
            if (k.eq.20) then
              itransp = 1
            endif
            if (k.eq.22) then
              itransp = 1
            endif
      end
