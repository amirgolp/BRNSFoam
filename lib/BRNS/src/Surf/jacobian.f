c      
c     SUBROUTINE jacobian
c      
      subroutine jacobian(pd,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension pd(ncomp,ncomp)
        call switches(j)
         pd(2,2) = 1
         pd(1,2) = 0
         pd(2,1) = 1
         pd(1,1) = -k_ads-1/delt
      end
