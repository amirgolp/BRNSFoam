c      
c     SUBROUTINE residual
c      
      subroutine residual(funcs,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension funcs(ncomp)
        call switches(j)
          funcs(1) = -k_ads*sp(1,j)-sp(1,j)/delt+spold(1,j)/delt
          funcs(2) = sp(1,j)-spold(1,j)+sp(2,j)-spold(2,j)
      end
