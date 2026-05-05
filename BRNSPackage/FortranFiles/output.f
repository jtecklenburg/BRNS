c      
c     SUBROUTINE out
c      
      subroutine out(j,nt,time,depth,v_out,v_int)
        include 'common_geo.inc'
        include 'common.inc'
        real*8 time
        if (nt.eq.1.and.j.eq.1) then
          open(unit=11,file='Bfo1.dat',status='replace')
          close(11)
          open(unit=11,file='Bmo1.dat',status='replace')
          close(11)
          open(unit=11,file='doc1.dat',status='replace')
          close(11)
          open(unit=11,file='dox1.dat',status='replace')
          close(11)
          open(unit=11,file='Amm1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifo1.dat',status='replace')
          close(11)
          open(unit=11,file='Bimo1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfn1.dat',status='replace')
          close(11)
          open(unit=11,file='Bmn1.dat',status='replace')
          close(11)
          open(unit=11,file='nitra1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifn1.dat',status='replace')
          close(11)
          open(unit=11,file='Bimn1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfs1.dat',status='replace')
          close(11)
          open(unit=11,file='Bms1.dat',status='replace')
          close(11)
          open(unit=11,file='sulpha1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifs1.dat',status='replace')
          close(11)
          open(unit=11,file='Bims1.dat',status='replace')
          close(11)
          open(unit=11,file='POM1.dat',status='replace')
          close(11)
          open(unit=11,file='tr1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfa1.dat',status='replace')
          close(11)
          open(unit=11,file='Bma1.dat',status='replace')
          close(11)
          open(unit=11,file='Bifa1.dat',status='replace')
          close(11)
          open(unit=11,file='Bima1.dat',status='replace')
          close(11)
          open(unit=11,file='Bfo1.inp',status='replace')
          close(11)
          open(unit=11,file='Bmo1.inp',status='replace')
          close(11)
          open(unit=11,file='doc1.inp',status='replace')
          close(11)
          open(unit=11,file='dox1.inp',status='replace')
          close(11)
          open(unit=11,file='Amm1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifo1.inp',status='replace')
          close(11)
          open(unit=11,file='Bimo1.inp',status='replace')
          close(11)
          open(unit=11,file='Bfn1.inp',status='replace')
          close(11)
          open(unit=11,file='Bmn1.inp',status='replace')
          close(11)
          open(unit=11,file='nitra1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifn1.inp',status='replace')
          close(11)
          open(unit=11,file='Bimn1.inp',status='replace')
          close(11)
          open(unit=11,file='Bfs1.inp',status='replace')
          close(11)
          open(unit=11,file='Bms1.inp',status='replace')
          close(11)
          open(unit=11,file='sulpha1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifs1.inp',status='replace')
          close(11)
          open(unit=11,file='Bims1.inp',status='replace')
          close(11)
          open(unit=11,file='POM1.inp',status='replace')
          close(11)
          open(unit=11,file='tr1.inp',status='replace')
          close(11)
          open(unit=11,file='Bfa1.inp',status='replace')
          close(11)
          open(unit=11,file='Bma1.inp',status='replace')
          close(11)
          open(unit=11,file='Bifa1.inp',status='replace')
          close(11)
          open(unit=11,file='Bima1.inp',status='replace')
          close(11)
        v_out = 10000.D0
        v_int = 10000.D0
        endif
        if (time.le.v_out.and.v_out.lt.time+delt) then
          open(unit=11,file='Bfo1.dat',
     +      status='old',access='append')
          write(11,2000) sp(1,j),depth
 2000     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmo1.dat',
     +      status='old',access='append')
          write(11,2001) sp(2,j),depth
 2001     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='doc1.dat',
     +      status='old',access='append')
          write(11,2002) sp(3,j),depth
 2002     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='dox1.dat',
     +      status='old',access='append')
          write(11,2003) sp(4,j),depth
 2003     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Amm1.dat',
     +      status='old',access='append')
          write(11,2004) sp(5,j),depth
 2004     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifo1.dat',
     +      status='old',access='append')
          write(11,2005) sp(6,j),depth
 2005     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimo1.dat',
     +      status='old',access='append')
          write(11,2006) sp(7,j),depth
 2006     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfn1.dat',
     +      status='old',access='append')
          write(11,2007) sp(8,j),depth
 2007     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmn1.dat',
     +      status='old',access='append')
          write(11,2008) sp(9,j),depth
 2008     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='nitra1.dat',
     +      status='old',access='append')
          write(11,2009) sp(10,j),depth
 2009     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifn1.dat',
     +      status='old',access='append')
          write(11,2010) sp(11,j),depth
 2010     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimn1.dat',
     +      status='old',access='append')
          write(11,2011) sp(12,j),depth
 2011     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfs1.dat',
     +      status='old',access='append')
          write(11,2012) sp(13,j),depth
 2012     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bms1.dat',
     +      status='old',access='append')
          write(11,2013) sp(14,j),depth
 2013     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='sulpha1.dat',
     +      status='old',access='append')
          write(11,2014) sp(15,j),depth
 2014     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifs1.dat',
     +      status='old',access='append')
          write(11,2015) sp(16,j),depth
 2015     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bims1.dat',
     +      status='old',access='append')
          write(11,2016) sp(17,j),depth
 2016     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='POM1.dat',
     +      status='old',access='append')
          write(11,2017) sp(18,j),depth
 2017     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='tr1.dat',
     +      status='old',access='append')
          write(11,2018) sp(19,j),depth
 2018     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfa1.dat',
     +      status='old',access='append')
          write(11,2019) sp(20,j),depth
 2019     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bma1.dat',
     +      status='old',access='append')
          write(11,2020) sp(21,j),depth
 2020     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifa1.dat',
     +      status='old',access='append')
          write(11,2021) sp(22,j),depth
 2021     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bima1.dat',
     +      status='old',access='append')
          write(11,2022) sp(23,j),depth
 2022     format(1x,e14.7,2x,f12.4)
          close(11)
        if (j.eq.nx) then
        v_out = v_out+v_int
        endif
        endif
        if (time.eq.endt) then
          open(unit=11,file='Bfo1.inp',
     +      status='old',access='append')
          write(11,2023) sp(1,j),depth
 2023     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmo1.inp',
     +      status='old',access='append')
          write(11,2024) sp(2,j),depth
 2024     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='doc1.inp',
     +      status='old',access='append')
          write(11,2025) sp(3,j),depth
 2025     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='dox1.inp',
     +      status='old',access='append')
          write(11,2026) sp(4,j),depth
 2026     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Amm1.inp',
     +      status='old',access='append')
          write(11,2027) sp(5,j),depth
 2027     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifo1.inp',
     +      status='old',access='append')
          write(11,2028) sp(6,j),depth
 2028     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimo1.inp',
     +      status='old',access='append')
          write(11,2029) sp(7,j),depth
 2029     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfn1.inp',
     +      status='old',access='append')
          write(11,2030) sp(8,j),depth
 2030     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bmn1.inp',
     +      status='old',access='append')
          write(11,2031) sp(9,j),depth
 2031     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='nitra1.inp',
     +      status='old',access='append')
          write(11,2032) sp(10,j),depth
 2032     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifn1.inp',
     +      status='old',access='append')
          write(11,2033) sp(11,j),depth
 2033     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bimn1.inp',
     +      status='old',access='append')
          write(11,2034) sp(12,j),depth
 2034     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfs1.inp',
     +      status='old',access='append')
          write(11,2035) sp(13,j),depth
 2035     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bms1.inp',
     +      status='old',access='append')
          write(11,2036) sp(14,j),depth
 2036     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='sulpha1.inp',
     +      status='old',access='append')
          write(11,2037) sp(15,j),depth
 2037     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifs1.inp',
     +      status='old',access='append')
          write(11,2038) sp(16,j),depth
 2038     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bims1.inp',
     +      status='old',access='append')
          write(11,2039) sp(17,j),depth
 2039     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='POM1.inp',
     +      status='old',access='append')
          write(11,2040) sp(18,j),depth
 2040     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='tr1.inp',
     +      status='old',access='append')
          write(11,2041) sp(19,j),depth
 2041     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bfa1.inp',
     +      status='old',access='append')
          write(11,2042) sp(20,j),depth
 2042     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bma1.inp',
     +      status='old',access='append')
          write(11,2043) sp(21,j),depth
 2043     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bifa1.inp',
     +      status='old',access='append')
          write(11,2044) sp(22,j),depth
 2044     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Bima1.inp',
     +      status='old',access='append')
          write(11,2045) sp(23,j),depth
 2045     format(1x,e14.7,2x,f12.4)
          close(11)
        endif
      end
