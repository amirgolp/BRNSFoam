c      
c     SUBROUTINE rates
c      
      subroutine rates(j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            r(1,j) = miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j)
     +)*sp(4,j)
            r(2,j) = nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO
     ++sp(3,j))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)
            r(3,j) = b_S*sp(4,j)
      end
