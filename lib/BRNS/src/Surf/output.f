c      
c     SUBROUTINE out
c      
      subroutine out(j,nt,time,depth,v_out,v_int)
        include 'common_geo.inc'
        include 'common.inc'
        real*8 time
        if (nt.eq.1.and.j.eq.1) then
          open(unit=11,file='Tracer.dat',status='replace')
          close(11)
          open(unit=11,file='Surf_wT.dat',status='replace')
          close(11)
          open(unit=11,file='rate_a.dat',status='replace')
          close(11)
          open(unit=11,file='Tracer.inp',status='replace')
          close(11)
          open(unit=11,file='Surf_wT.inp',status='replace')
          close(11)
        v_out = 50
        v_int = 500
        endif
        if (time.le.v_out.and.v_out.lt.time+delt) then
          open(unit=11,file='Tracer.dat',
     +      status='old',access='append')
          write(11,2000) sp(1,j),depth
 2000     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Surf_wT.dat',
     +      status='old',access='append')
          write(11,2001) sp(2,j),depth
 2001     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='rate_a.dat',
     +      status='old',access='append')
          write(11,2002) r(1,j),depth
 2002     format(1x,e14.7,2x,f12.4)
          close(11)
        if (j.eq.nx) then
        v_out = v_out+v_int
        endif
        endif
        if (time.eq.endt) then
          open(unit=11,file='Tracer.inp',
     +      status='old',access='append')
          write(11,2003) sp(1,j),depth
 2003     format(1x,e14.7,2x,f12.4)
          close(11)
          open(unit=11,file='Surf_wT.inp',
     +      status='old',access='append')
          write(11,2004) sp(2,j),depth
 2004     format(1x,e14.7,2x,f12.4)
          close(11)
        endif
      end
