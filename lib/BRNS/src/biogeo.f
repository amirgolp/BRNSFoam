c      
c     SUBROUTINE biogeo
c      
      subroutine biogeo()
        include 'common_geo.inc'
        include 'common.inc'
          miu_max = 1.D0/28800.D0
          b_S = 0.7175925926D-5
          K_S = 0.3D0
          K_O2 = 0.625D-2
          K_NO = 0.357D-1
          K_i_O2 = 0.3125D-2
          nu_anox = 0.8D0
      end
