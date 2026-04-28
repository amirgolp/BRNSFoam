c      
c     SUBROUTINE boundaries
c      
      subroutine boundaries()
        include 'common_geo.inc'
        include 'common.inc'
        j = 1
          spb(1,1) = 0.1D1
          spb(2,1) = 1
          ibc(1,1) = 0
          ibc(2,1) = 0
        j = nx
          spb(1,2) = 0.D0
          spb(2,2) = 0
          ibc(1,2) = 1
          ibc(2,2) = 1
      end
