      subroutine fpbspl(t,n,k,x,l,h)
!c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
!c  degree k at t(l) <= x < t(l+1) using the stable recurrence
!c  relation of de boor and cox.
!c  Travis Oliphant  2007
!c    changed so that weighting of 0 is used when knots with
!c      multiplicity are present.
!c    Also, notice that l+k <= n and 1 <= l+1-k
!c      or else the routine will be accessing memory outside t
!c      Thus it is imperative that that k <= l <= n-k but this
!c      is not checked.
!c  ..
!c  ..scalar arguments..
      real*8 x
      integer n,k,l
!c  ..array arguments..
      real*8 t(n),h(20)
!c  ..local scalars..
      real*8 f,one
      integer i,j,li,lj
!c  ..local arrays..
      real*8 hh(19)
!c  ..
      one = 0.1d+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.0d0
        do 20 i=1,j
          li = l+i
          lj = li-j
          if (t(li).ne.t(lj)) goto 15
          h(i+1) = 0.0d0 
          goto 20
  15      f = hh(i)/(t(li)-t(lj)) 
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end subroutine fpbspl


      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
!c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my
!c  ..array arguments..
      integer lx(mx),ly(my)
!ahsan      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
!     * wx(mx,kx+1),wy(my,ky+1)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my), wx(mx,kx+1),wy(my,ky+1)
!c  ..local scalars..
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      real*8 arg,sp,tb,te
!c  ..local arrays..
      real*8 h(6)
!c  ..subroutine references..
!c    fpbspl
!c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
  30    continue
  40  continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
  90    continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100        continue
            l1 = l1+nky1
 110      continue
          m = m+1
          z(m) = sp
 120    continue
 130  continue
      return
      end subroutine fpbisp

!*************************************************************************
!**************************************************************************
      subroutine splev(t,n,c,k,x,y,m,e,ier)
!c subroutine splev evaluates in a number of points x(i),i=1,2,...,m
!c a spline s(x) of degree k, given in its b-spline representation.
!c
!c calling sequence:
!c    call splev(t,n,c,k,x,y,m,e,ier)
!c
!c input parameters:
!c   t    : array,length n, which contains the position of the knots.
!c   n    : integer, giving the total number of knots of s(x).
!c   c    : array,length n, which contains the b-spline coefficients.
!c   k    : integer, giving the degree of s(x).
!c   x    : array,length m, which contains the points where s(x) must
!c          be evaluated.
!c   m    : integer, giving the number of points where s(x) must be
!c          evaluated.
!c   e    : integer, if 0 the spline is extrapolated from the end
!c          spans for points not in the support, if 1 the spline
!c          evaluates to zero for those points, if 2 ier is set to
!c          1 and the subroutine returns, and if 3 the spline evaluates
!c          to the value of the nearest boundary point.
!c
!c output parameter:
!c   y    : array,length m, giving the value of s(x) at the different
!c          points.
!c   ier  : error flag
!c     ier = 0 : normal return
!c     ier = 1 : argument out of bounds and e == 2
!c     ier =10 : invalid input data (see restrictions)
!c
!c restrictions:
!c   m >= 1
!c--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!c
!c other subroutines required: fpbspl.
!c
!c references :
!c   de boor c : on calculating with b-splines, j. approximation theory
!c                6 (1972) 50-62.
!c   cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!c                applics 10 (1972) 134-149.
!c   dierckx p. : curve and surface fitting with splines, monographs on
!c                numerical analysis, oxford university press, 1993.
!c
!c author :
!c    p.dierckx
!c    dept. computer science, k.u.leuven
!c    celestijnenlaan 200a, b-3001 heverlee, belgium.
!c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!c
!c  latest update : march 1987
!c
!c++ pearu: 11 aug 2003
!c++   - disabled cliping x values to interval [min(t),max(t)]
!c++   - removed the restriction of the orderness of x values
!c++   - fixed initialization of sp to double precision value
!c
!c  ..scalar arguments..
      integer n, k, m, e, ier
!c  ..array arguments..
      real*8 t(n), c(n), x(m), y(m)
!c  ..local scalars..
      integer i, j, k1, l, ll, l1, nk1
!c++..
      integer k2
!c..++
      real*8 arg, sp, tb, te
!c  ..local array..
      real*8 h(20)
!c  ..
!c  before starting computations a data check is made. if the input data
!c  are invalid control is immediately repassed to the calling program.
      ier = 10
!c--      if(m-1) 100,30,10
!c++..
      if (m .lt. 1) go to 100
!c..++
!c--  10  do 20 i=2,m
!c--        if(x(i).lt.x(i-1)) go to 100
!c--  20  continue
  30  ier = 0
!c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k + 1
!c++..
      k2 = k1 + 1
!c..++
      nk1 = n - k1
      tb = t(k1)
      te = t(nk1 + 1)
      l = k1
      l1 = l + 1
!c  main loop for the different points.
      do 80 i = 1, m
!c  fetch a new x-value arg.
        arg = x(i)
!c  check if arg is in the support
        if (arg .lt. tb .or. arg .gt. te) then
            if (e .eq. 0) then
                goto 35
            else if (e .eq. 1) then
                y(i) = 0
                goto 80
            else if (e .eq. 2) then
                ier = 1
                goto 100
            else if (e .eq. 3) then
                if (arg .lt. tb) then
                    arg = tb
                else
                    arg = te
                endif
            endif
        endif
!c  search for knot interval t(l) <= arg < t(l+1)
!c++..
 35     if (arg .ge. t(l) .or. l1 .eq. k2) go to 40
        l1 = l
        l = l - 1
        go to 35
!c..++
  40    if(arg .lt. t(l1) .or. l .eq. nk1) go to 50
        l = l1
        l1 = l + 1
        go to 40
!c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t, n, k, arg, l, h)
!c  find the value of s(x) at x=arg.
        sp = 0.0d0
        ll = l - k1
        do 60 j = 1, k1
          ll = ll + 1
          sp = sp + c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end subroutine splev

     subroutine splev_cpp(t,n,c,k,x,y,m,ier)
     implicit none 
       integer, intent(in):: n
       real*8, dimension(n),intent(in) :: t
       real*8, dimension(n),intent(in) :: c
       integer, intent(in) :: k
       integer, intent(in) :: m
       real*8, dimension(m),intent(in) :: x
       real*8, dimension(m),intent(out) :: y
       integer :: e
       integer, intent(out) :: ier
       
       e=2 !see splev for explanation
       call splev(t,n,c,k,x,y,m,e,ier)
     end subroutine splev_cpp


!*************************************************************************
!**************************************************************************
      subroutine bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk, ier)
!c  subroutine bispeu evaluates on a set of points (x(i),y(i)),i=1,...,m
!c  a bivariate spline s(x,y) of degrees kx and ky, given in the
!c  b-spline representation.
!c
!c  calling sequence:
!c     call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,
!c    * iwrk,kwrk,ier)
!c
!c  input parameters:
!c   tx    : real array, length nx, which contains the position of the
!c           knots in the x-direction.
!c   nx    : integer, giving the total number of knots in the x-direction
!c   ty    : real array, length ny, which contains the position of the
!c           knots in the y-direction.
!c   ny    : integer, giving the total number of knots in the y-direction
!c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!c           b-spline coefficients.
!c   kx,ky : integer values, giving the degrees of the spline.
!c   x     : real array of dimension (mx).
!c   y     : real array of dimension (my).
!c   m     : on entry m must specify the number points. m >= 1.
!c   wrk   : real array of dimension lwrk. used as workspace.
!c   lwrk  : integer, specifying the dimension of wrk.
!c           lwrk >= kx+ky+2
!c
!c  output parameters:
!c   z     : real array of dimension m.
!c           on succesful exit z(i) contains the value of s(x,y)
!c           at the point (x(i),y(i)), i=1,...,m.
!c   ier   : integer error flag
!c    ier=0 : normal return
!c    ier=10: invalid input data (see restrictions)
!c
!c  restrictions:
!c   m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
!c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!c
!c  other subroutines required:
!c    fpbisp,fpbspl
!c
!c  ..scalar arguments..
      integer nx,ny,kx,ky,m,lwrk,kwrk,ier
!c  ..array arguments..
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(m),y(m),z(m),wrk(lwrk)
!c  ..local scalars..
      integer iwrk(2)
      integer i,iw,lwest
!c  ..
!c  before starting computations a data check is made. if the input data
!c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = kx+ky+2
      if (lwrk.lt.lwest) go to 100
      if (m.lt.1) go to 100
      ier = 0
      do 10 i=1,m
         call fpbisp(tx,nx,ty,ny,c,kx,ky,x(i),1,y(i),1,z(i),wrk(1),&
                     wrk(kx+2),iwrk(1),iwrk(2))
 10   continue
 100  return
      end subroutine bispeu

     subroutine bispeu_cpp(tx,nx,ty,ny,c,nc,kx,ky,x,y,z,m,wrk,lwrk,ier)
     implicit none
 
       integer, intent(in) :: nx
       real*8, dimension(nx),intent(in) :: tx
       integer, intent(in) :: ny
       real*8, dimension(ny),intent(in) :: ty
       
       integer, intent(in) :: kx
       integer, intent(in) :: ky
       integer, intent(in) :: nc
       real*8, intent(in), dimension(nc):: c !(nx-kx-1)*(ny-ky-1)
       
       integer, intent(in) :: m !len(x)
       real*8, intent(in),dimension(m) :: x
       real*8, intent(in),dimension(m) :: y
       real*8, dimension(m),intent(out) :: z
       
       integer, intent(in) :: lwrk !=kx+ky+2
       real*8, dimension(lwrk) :: wrk !wrk is the workspace  required by bispeu
      
       integer, intent(out) :: ier
       
       call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,ier)
     end subroutine bispeu_cpp

!**************************************************************************
!**************************************************************************

! 
!  program myprog
! 	implicit none
! 	
! 	!REAL, PARAMETER :: Pi = 3.1415927
! 	integer, parameter :: n = 29
! 	integer, parameter ::  k = 3
! 	real*8, dimension(n):: t, c
! 	integer, parameter :: size_x = 3
! 	real*8, dimension(size_x):: x,y
! 	integer i, ier
! 	!!!!
! 	integer, parameter :: size_tx = 215
! 	integer, parameter :: size_ty = 303
! 	integer, parameter :: size_2dc = 63089
! 	real*8, dimension(size_tx) :: tx
! 	real*8, dimension(size_ty) :: ty
! 	real*8, dimension(size_2dc) :: c2d
! 	real*8, dimension(4) :: x2d_vals, y2d_vals
! 	real*8, dimension(4) :: z2d_vals
! 	real*8, dimension(10) :: wrk2d
! 	
! 	open (unit = 1, file = "product_tc.dat")
! 	do i = 1, n
! 		read (1,*) t(i), c(i)
! 	enddo
! 	close(1)
! 	
! 	x(:) = (/4.22,-2.3,0.8/)
! 	
! 	!splev(t,n,c,k,x,y,m,e,ier)
! 	call splev(t,n,c,k,x,y,size_x,2,ier)
! 	write(*,*) 'x_vals ',x
! 	write(*,*) 'y_vals ',y 
! 	write(*,*) 'ier ', ier
! 	write(*,*) t(3), c(5)
! 	write(*,*) "How you doing bro?"
! 	
! 	do i = 1, 900000
! 		call splev(t,n,c,k,x(1),y(1),1,2,ier)
! 		call splev(t,n,c,k,x(2),y(2),1,2,ier)
! 		call splev(t,n,c,k,x(3),y(3),1,2,ier)
! 	end do
! 	
! 	open (unit = 2, file = "product_2dtx.dat")
! 	do i = 1, size_tx
! 		read (2,*) tx(i)
! 	enddo
! 	close(2)
! 	
! 	open (unit = 3, file = "product_2dty.dat")
! 	do i = 1, size_ty
! 		read (3,*) ty(i)
! 	enddo
! 	close(3)
! 
! 	open (unit = 4, file = "product_2dc.dat")
! 	do i = 1, size_2dc
! 		read (4,*) c2d(i)
! 	enddo
! 	close(4)
! 	
! 	!call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk, ier)
! 	x2d_vals = (/-4.33,0.7,2.34,5./) 
! 	y2d_vals = (/-6.,3.5,6.89,-1./)
! 	call bispeu(tx,size_tx,ty,size_ty,c2d,3,3,x2d_vals,y2d_vals,z2d_vals,4,wrk2d,10, ier)
! 	write(*,*) " 2d ier ", ier
! 	write(*,*) "x, y, spline_val, real_val"
! 	do i=1,4
! 		write(*,*) x2d_vals(i), y2d_vals(i), z2d_vals(i)/dexp(-x2d_vals(i)**2+3.*y2d_vals(i))
! 	end do
! 	
! 		do i = 1, 900000
! call bispeu(tx,size_tx,ty,size_ty,c2d,3,3,x2d_vals,y2d_vals,z2d_vals,4,wrk2d,10, ier)
! 	end do
!  end program myprog
 
 
 
 
 
 
 
 
 
 
 
 
