        implicit real*8(a-h,o-z) 
        parameter (nsolid=0,ndiss=2,ncomp=2,nreac=1,nx=101)
        common/kinetics/k_ads
        real*8 k_ads
        common/physics/al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,dept
     +hmax,endt
        real*8 al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,end
     +t
        common/physics2/iq,iw,iDb,ipor,igrid,iarea,ic
        integer iq,iw,iDb,ipor,igrid,iarea,ic
