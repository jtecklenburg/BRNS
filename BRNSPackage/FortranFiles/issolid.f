c      
c     SUBROUTINE issolid
c      
      subroutine issolid(k,isolid)
        include 'common_geo.inc'
        include 'common.inc'
            if (k.eq.1) then
              isolid = 1
            endif
            if (k.eq.6) then
              isolid = 1
            endif
            if (k.eq.8) then
              isolid = 1
            endif
            if (k.eq.11) then
              isolid = 1
            endif
            if (k.eq.13) then
              isolid = 1
            endif
            if (k.eq.16) then
              isolid = 1
            endif
            if (k.eq.20) then
              isolid = 1
            endif
            if (k.eq.22) then
              isolid = 1
            endif
      end
