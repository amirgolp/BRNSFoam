c      
c     SUBROUTINE ssrates
c      
      subroutine ssrates(rat,drdc,isp,j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            if (isp.eq.1) then
            rat = -k_ads*sp(1,j)
            drdc = -k_ads
            endif
            if (isp.eq.2) then
            rat = k_ads*sp(1,j)
            drdc = 0
            endif
      end
