c      
c     SUBROUTINE initialcond
c      
      subroutine initialcond()
          include 'common_geo.inc'
          include 'common.inc'
          include 'common_drive.inc'
        real*8 spi(ncomp)
          if (ic.eq.1) then
          open(unit=85,file='initialconc.txt',status='old')
c           
          do 1000, j=1,nx,2
            read(85,*) de,(sp(i,j),i=1,ncomp) 
 1000     continue
c         
          close(85)
          endif
          if (ic.eq.2) then
              spi(1) = 0.2D-5
              spi(2) = 0.2D-5
              spi(3) = 625.D0
              spi(4) = 365.D0
              spi(5) = 0.33D1
              spi(6) = 0.2D-5
              spi(7) = 0.2D-5
              spi(8) = 0.2D-5
              spi(9) = 0.2D-5
              spi(10) = 119.D0
              spi(11) = 0.2D-5
              spi(12) = 0.2D-5
              spi(13) = 0.2D-5
              spi(14) = 0.2D-5
              spi(15) = 144.D0
              spi(16) = 0.2D-5
              spi(17) = 0.2D-5
              spi(18) = 150.D0
              spi(19) = 15.D0
              spi(20) = 0.2D-5
              spi(21) = 0.2D-5
              spi(22) = 0.2D-5
              spi(23) = 0.2D-5
c             
            do 1001, i=1,ncomp
c               
              do 1002, j=1,nx,2
                sp(i,j) = spi(i)
 1002         continue
c             
 1001       continue
c           
          endif
          if (ic.eq.3) then
              open(unit=101,file='Bfo1.inp',status='old')
              open(unit=102,file='Bmo1.inp',status='old')
              open(unit=103,file='doc1.inp',status='old')
              open(unit=104,file='dox1.inp',status='old')
              open(unit=105,file='Amm1.inp',status='old')
              open(unit=106,file='Bifo1.inp',status='old')
              open(unit=107,file='Bimo1.inp',status='old')
              open(unit=108,file='Bfn1.inp',status='old')
              open(unit=109,file='Bmn1.inp',status='old')
              open(unit=110,file='nitra1.inp',status='old')
              open(unit=111,file='Bifn1.inp',status='old')
              open(unit=112,file='Bimn1.inp',status='old')
              open(unit=113,file='Bfs1.inp',status='old')
              open(unit=114,file='Bms1.inp',status='old')
              open(unit=115,file='sulpha1.inp',status='old')
              open(unit=116,file='Bifs1.inp',status='old')
              open(unit=117,file='Bims1.inp',status='old')
              open(unit=118,file='POM1.inp',status='old')
              open(unit=119,file='tr1.inp',status='old')
              open(unit=120,file='Bfa1.inp',status='old')
              open(unit=121,file='Bma1.inp',status='old')
              open(unit=122,file='Bifa1.inp',status='old')
              open(unit=123,file='Bima1.inp',status='old')
c             
            do 1003, j=1,nx,2
                read(101,2000) sp(1,j),depth
 2000           format(1x,e14.7,2x,f8.4)
                read(102,2001) sp(2,j),depth
 2001           format(1x,e14.7,2x,f8.4)
                read(103,2002) sp(3,j),depth
 2002           format(1x,e14.7,2x,f8.4)
                read(104,2003) sp(4,j),depth
 2003           format(1x,e14.7,2x,f8.4)
                read(105,2004) sp(5,j),depth
 2004           format(1x,e14.7,2x,f8.4)
                read(106,2005) sp(6,j),depth
 2005           format(1x,e14.7,2x,f8.4)
                read(107,2006) sp(7,j),depth
 2006           format(1x,e14.7,2x,f8.4)
                read(108,2007) sp(8,j),depth
 2007           format(1x,e14.7,2x,f8.4)
                read(109,2008) sp(9,j),depth
 2008           format(1x,e14.7,2x,f8.4)
                read(110,2009) sp(10,j),depth
 2009           format(1x,e14.7,2x,f8.4)
                read(111,2010) sp(11,j),depth
 2010           format(1x,e14.7,2x,f8.4)
                read(112,2011) sp(12,j),depth
 2011           format(1x,e14.7,2x,f8.4)
                read(113,2012) sp(13,j),depth
 2012           format(1x,e14.7,2x,f8.4)
                read(114,2013) sp(14,j),depth
 2013           format(1x,e14.7,2x,f8.4)
                read(115,2014) sp(15,j),depth
 2014           format(1x,e14.7,2x,f8.4)
                read(116,2015) sp(16,j),depth
 2015           format(1x,e14.7,2x,f8.4)
                read(117,2016) sp(17,j),depth
 2016           format(1x,e14.7,2x,f8.4)
                read(118,2017) sp(18,j),depth
 2017           format(1x,e14.7,2x,f8.4)
                read(119,2018) sp(19,j),depth
 2018           format(1x,e14.7,2x,f8.4)
                read(120,2019) sp(20,j),depth
 2019           format(1x,e14.7,2x,f8.4)
                read(121,2020) sp(21,j),depth
 2020           format(1x,e14.7,2x,f8.4)
                read(122,2021) sp(22,j),depth
 2021           format(1x,e14.7,2x,f8.4)
                read(123,2022) sp(23,j),depth
 2022           format(1x,e14.7,2x,f8.4)
 1003       continue
c           
              close(101)
              close(102)
              close(103)
              close(104)
              close(105)
              close(106)
              close(107)
              close(108)
              close(109)
              close(110)
              close(111)
              close(112)
              close(113)
              close(114)
              close(115)
              close(116)
              close(117)
              close(118)
              close(119)
              close(120)
              close(121)
              close(122)
              close(123)
          endif
      end
