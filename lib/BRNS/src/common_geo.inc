        implicit real*8(a-h,o-z) 
        parameter (nsolid=0,ndiss=4,ncomp=4,nreac=3,nx=101)
        common/kinetics/miu_max,b_S,K_S,K_O2,K_NO,K_i_O2,nu_anox
        real*8 miu_max,b_S,K_S,K_O2,K_NO,K_i_O2,nu_anox
        common/physics/al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,dept
     +hmax,endt
        real*8 al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,end
     +t
        common/physics2/iq,iw,iDb,ipor,igrid,iarea,ic
        integer iq,iw,iDb,ipor,igrid,iarea,ic
