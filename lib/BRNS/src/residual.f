c      
c     SUBROUTINE residual
c      
      subroutine residual(funcs,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension funcs(ncomp)
        call switches(j)
          funcs(1) = -miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j
     +))*sp(4,j)-sp(2,j)/delt+spold(2,j)/delt
          funcs(2) = nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO
     ++sp(3,j))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)-(-sp(1,j)+sp(2,j))/delt+
     +(-spold(1,j)+spold(2,j))/delt
          funcs(3) = -4.D0/5.D0*sp(1,j)+4.D0/5.D0*spold(1,j)+4.D0/5.D0*s
     +p(2,j)-4.D0/5.D0*spold(2,j)+sp(3,j)-spold(3,j)
          funcs(4) = -b_S*sp(4,j)-(sp(1,j)+sp(4,j))/delt+(spold(1,j)+spo
     +ld(4,j))/delt
      end
