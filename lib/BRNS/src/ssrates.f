c      
c     SUBROUTINE ssrates
c      
      subroutine ssrates(rat,drdc,isp,j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            if (isp.eq.1) then
            rat = -miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))*
     +sp(4,j)-nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO+sp(3,j
     +))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)
            drdc = -miu_max/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))*sp(4,j)
     ++miu_max*sp(1,j)/(K_S+sp(1,j))**2*sp(2,j)/(K_O2+sp(2,j))*sp(4,j)-n
     +u_anox*miu_max/(K_S+sp(1,j))*sp(3,j)/(K_NO+sp(3,j))*K_i_O2/(K_i_O2
     ++sp(2,j))*sp(4,j)+nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))**2*sp(3,j)
     +/(K_NO+sp(3,j))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)
            endif
            if (isp.eq.2) then
            rat = -miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))*
     +sp(4,j)
            drdc = -miu_max*sp(1,j)/(K_S+sp(1,j))/(K_O2+sp(2,j))*sp(4,j)
     ++miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))**2*sp(4,j)
            endif
            if (isp.eq.3) then
            rat = -0.8D0*nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(
     +K_NO+sp(3,j))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)
            drdc = -0.8D0*nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))/(K_NO+sp
     +(3,j))*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)+0.8D0*nu_anox*miu_max*sp(1,
     +j)/(K_S+sp(1,j))*sp(3,j)/(K_NO+sp(3,j))**2*K_i_O2/(K_i_O2+sp(2,j))
     +*sp(4,j)
            endif
            if (isp.eq.4) then
            rat = miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))*s
     +p(4,j)+nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO+sp(3,j)
     +)*K_i_O2/(K_i_O2+sp(2,j))*sp(4,j)-b_S*sp(4,j)
            drdc = miu_max*sp(1,j)/(K_S+sp(1,j))*sp(2,j)/(K_O2+sp(2,j))+
     +nu_anox*miu_max*sp(1,j)/(K_S+sp(1,j))*sp(3,j)/(K_NO+sp(3,j))*K_i_O
     +2/(K_i_O2+sp(2,j))-b_S
            endif
      end
