c      
c     SUBROUTINE boundaries
c      
      subroutine boundaries()
        include 'common_geo.inc'
        include 'common.inc'
        j = 1
          spb(1,1) = 1
          spb(2,1) = 1
          spb(3,1) = 1
          spb(4,1) = 1
          ibc(1,1) = 1
          ibc(2,1) = 1
          ibc(3,1) = 1
          ibc(4,1) = 1
        j = nx
          spb(1,2) = 0
          spb(2,2) = 0
          spb(3,2) = 0
          spb(4,2) = 0
          ibc(1,2) = 0
          ibc(2,2) = 0
          ibc(3,2) = 0
          ibc(4,2) = 0
      end
