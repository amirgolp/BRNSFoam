c      
c     SUBROUTINE jacobian
c      
      subroutine jacobian(pd,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension pd(ncomp,ncomp)
        call switches(j)
         pd(3,4) = 0
         pd(1,2) = -miu_max*sp(1,j)/(K_S+sp(1,j))/(K_O2+sp(2,j))*sp(4,j)
     ++miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))**2*sp(4,j)-1
     +/delt
         pd(4,4) = -b_S-1/delt
         pd(3,1) = -4.D0/5.D0
         pd(4,1) = -1/delt
         pd(3,3) = 1
         pd(1,4) = -miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))
         pd(1,1) = -miu_max/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))*sp(4,j)
     ++miu_max*sp(1,j)/(K_S+sp(1,j))**2*sp(2,j)/(K_O2+sp(2,j))*sp(4,j)
         pd(1,3) = 0
         pd(2,2) = -nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO+
     +sp(3,j))*K_i_O2/(K_i_O2+sp(2,j))**2*sp(4,j)-1/delt
         pd(2,4) = nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO+s
     +p(3,j))*K_i_O2/(K_i_O2+sp(2,j))
         pd(3,2) = 4.D0/5.D0
         pd(2,3) = nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))/(K_NO+sp(3,j))*
     +K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)-nu_anox*miu_max*sp(1,j)/(K_S+sp(1,
     +j))*sp(3,j)/(K_NO+sp(3,j))**2*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)
         pd(4,3) = 0
         pd(4,2) = 0
         pd(2,1) = nu_anox*miu_max/(K_S+sp(1,j))*sp(3,j)/(K_NO+sp(3,j))*
     +K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)-nu_anox*miu_max*sp(1,j)/(K_S+sp(1,
     +j))**2*sp(3,j)/(K_NO+sp(3,j))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)+1/de
     +lt
      end
