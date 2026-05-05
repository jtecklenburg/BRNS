c      
c     SUBROUTINE readbiogeo
c      
      subroutine readbiogeo()
        include 'common_geo.inc'
        include 'common.inc'
	  
	open(unit=999,file='column-param.in',status='old')

          read(999,*) k_exch_B1
          write(*,*) k_exch_B1
          read(999,*) k_exch_B2
          write(*,*) k_exch_B2
          read(999,*) k_exch_B3
          write(*,*) k_exch_B3
          read(999,*) dec_t_B1
          write(*,*) dec_t_B1
          read(999,*) dec_t_B2
          write(*,*) dec_t_B2
          read(999,*) dec_t_B3
          write(*,*) dec_t_B3
          read(999,*) alpha_tol_b
          write(*,*) alpha_tol_b
          read(999,*) dummy1
          kmax_l_B1=kmax_l_B1*dummy1
          write(*,*) kmax_l_B1
          read(999,*) dummy2
          kmax_l_B2=kmax_l_B2*dummy2
          write(*,*) kmax_l_B2
          read(999,*) dummy3
          kmax_l_B3=kmax_l_B3*dummy3
          write(*,*) kmax_l_B3

	close(999)

      end
