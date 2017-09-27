
c  a front end to several numerical recipes routines
c  for least square fitting - straight line, errors
c  in 0, 1, or 2 coords

c lstsq.f is copyright 2006 Benjamin Weiner and released under
c the GNU General Public License.  You may modify, you may redistribute,
c you must retain the copyright notice.  Please contact me with
c any comments and suggestions.

c See Weiner et al. 2006 (ApJ, in press, astro-ph/------)
c for a discussion of fitting lines to data with intrinsic scatter.

c 4/23/00 - modified to ignore input lines beginning with '#',
c   and to allow fitting only a given region

c 12/20/05 - allows adding the intrinsic scatter in quadrature
c    Iterate until you find the value of intrinsic scatter that
c    gives chisq/N = 1.  At this value, the error bars on a and b
c    should be reliable.  If you _don't_ add intrinsic scatter,
c    the errors will be underestimates.

c Benjamin Weiner, bjw at astro.umd.edu (Astronomy Dept, University
c of Maryland; Steward Observatory as of Sept 2006)

      program lstsq

      parameter(NMAX=10000)

c      real x1(NMAX),x2(NMAX),x3(NMAX)x4(NMAX)
c      real x1err(NMAX),x2err(NMAX),x3err(NMAX)
      real x(NMAX), xerr(NMAX)
      real y(NMAX), yerr(NMAX)
      real xf(NMAX), xferr(NMAX)
      real yf(NMAX), yferr(NMAX)
      real xerrbot(NMAX),xerrtop(NMAX)
      real yerrbot(NMAX),yerrtop(NMAX)
      real xpl(2),ypl(2)

      character fname*40,dline*80

      big = 1.e20

      write(*,*) "0- plain least squares, 1- errors in y"
c      write(*,'(" 2- errors in x and y : "$)')
      write(*,'(" 2- errors in x and y, ",
     +     "3- x & y errors and y scatter: "$)')
c      write(*,'(" 3- y=y(x1,x2,...) : "$)')
      read(*,*) ichoose

c      if (ichoose .eq. 3) then
c         write(*,'("number of x'es [1-4]: "
c         read (*,*) nx
c      end if
      
 10   write (*,'("enter file: "$)')
      read (*,*) fname

      open (2,file=fname,status='old',err=10)
      npts = 1
 20   continue
         read (2,'(a)',end=30) dline
         if (dline(1:1) .eq. '#') go to 20
         if (ichoose .eq. 0) then
            xerr(npts) = 0.
            yerr(npts) = 0.
            read (dline,*,err=20) x(npts),y(npts)
         else if (ichoose .eq. 1) then
            xerr(npts) = 0.
            read (dline,*,err=20) x(npts),y(npts),yerr(npts)
         else if (ichoose .eq. 2 .or. ichoose .eq. 3) then
            read (dline,*,err=20) x(npts),xerr(npts),
     $           y(npts),yerr(npts)
         end if
         npts=npts+1
      go to 20
      
 30   npts = npts-1
      write(*,*) 'Read ', npts, ' points'
      close(2)
      if (npts .gt. NMAX) then
         write(*,*) "Currently dimensioned for only max ",NMAX," points"
      end if
      if (ichoose .ge. 2 .and. npts .gt. 20000) then
       write(*,*) 
     +   "Fitting w/both X and Y error will not work for >20000 points"
       write(*,*) "without redimensioning arrays in bigfitexy.f and"
       write(*,*) "bigchixy.f."
      end if

c for ichoose = 3, add the scatter in quadrature to the error
c and then continue as if we were doing ichoose=2
      if (ichoose .eq. 3) then
         write(*,'("intrinsic scatter in y-coord: ",$)')
         read (*,*) yintscatt
         yintscattsq = yintscatt**2
         do i=1,npts
            yerr(i) = sqrt(yerr(i)**2 + yintscattsq)
         end do
         ichoose = 2
      end if

c  moved plotting points to here

c this is sloppy and works because the errors are 0 if
c they are not supplied

      xmin = big
      xmax = -big
      ymin = big
      ymax = -big
      do i=1,npts
         xerrbot(i) = x(i)-xerr(i)
         xerrtop(i) = x(i)+xerr(i)
         yerrbot(i) = y(i)-yerr(i)
         yerrtop(i) = y(i)+yerr(i)
         xmin = min(xmin,xerrbot(i))
         xmax = max(xmax,xerrtop(i))
         ymin = min(ymin,yerrbot(i))
         ymax = max(ymax,yerrtop(i))
      end do

      xmargin = 0.1*(xmax-xmin)
      ymargin = 0.1*(ymax-ymin)
      wxmin = xmin-xmargin
      wxmax = xmax+xmargin
      wymin = ymin-ymargin
      wymax = ymax+ymargin

      call pgbegin (0,"?",1,1)
      call pgask(.false.)
      call pgscf(2)
      call pgsch(1.1)

      call pgenv(wxmin,wxmax,wymin,wymax,0,1)
      call pglabel ("x","y","Least squares fit")
      call pgpoint (npts,x,y,17)
      if(ichoose .eq. 1 .or. ichoose .eq. 2) then
         call pgerrb(6,npts,x,y,yerr,1.0)
         if (ichoose .eq. 2) then
            call pgerrb(5,npts,x,y,xerr,1.0)
         end if
      end if

c  select fitting window
      xfmin = wxmin
      xfmax = wxmax
      write (*,'("Enter xmin, xmax to fit [use all]: ",$)')
      read (*,'(a)') dline
      if (dline(1:5) .ne. '     ') then
         read(dline,*,err=40) xfmin,xfmax
      end if
 40   continue

      nf = 0
      do i=1,npts
         if(x(i) .gt. xfmin .and. x(i) .lt. xfmax) then
            nf = nf+1
            xf(nf) = x(i)
            yf(nf) = y(i)
            xferr(nf) = xerr(i)
            yferr(nf) = yerr(i)
         end if
      end do
      write(*,*) "Including ",nf," points"
      if(nf .lt.2) go to 666

c  do fit
      
      if (ichoose .eq. 0) then
c         call fit(x,y,npts,yerr,0,a,b,siga,sigb,chisq,q)
         call fit(xf,yf,nf,yferr,0,a,b,siga,sigb,chisq,q)
      else if (ichoose .eq. 1) then
c         call fit(x,y,npts,yerr,1,a,b,siga,sigb,chisq,q)
         call fit(xf,yf,nf,yferr,1,a,b,siga,sigb,chisq,q)
      else if (ichoose .eq. 2) then
c         call fitexy(x,y,npts,xerr,yerr,a,b,siga,sigb,chisq,q)
c         call fitexy_big(xf,yf,nf,xferr,yferr,a,b,siga,sigb,chisq,q)
         call bigfitexy(xf,yf,nf,xferr,yferr,a,b,siga,sigb,chisq,q)
      end if

c  plot fit

      xpl(1) = wxmin
      ypl(1) = a + b*wxmin
      xpl(2) = wxmax
      ypl(2) = a + b*wxmax

      call pgline(2,xpl,ypl)

      write(*,*) 'y = a + bx'
      write(*,*) 'a = ',a,' +/- ',siga
      write(*,*) 'b = ',b,' +/- ',sigb
      write(*,*) 'chisq = ',chisq,' prob. = ',q

      sumsqerr = 0.0
c      do i=1,npts
      do i=1,nf
         sumsqerr = sumsqerr + (y(i)-a-b*x(i))**2
      end do
c      write (*,*) 'rms error = ',sqrt(sumsqerr/npts)
      write (*,*) 'rms error = ',sqrt(sumsqerr/nf)

c  write data and fit to file

      open(3,file='lstsq.out',status='unknown')
      if (ichoose .eq. 0) then
         write (3,*) 'xdata       ydata       yfit'
         do i=1,npts
            write(3,*) x(i),y(i),a + b*x(i)
         end do
      else if (ichoose .eq. 1) then
         write (3,*) 'xdata       ydata        yerror       yfit'
         do i=1,npts
            write(3,*) x(i),y(i),yerr(i),a + b*x(i)
         end do
      else if (ichoose .eq. 2) then
         write (3,*) 
     $        'xdata       xerror       ydata       yerror   yfit'
         do i=1,npts
            write(3,*) x(i),xerr(i),y(i),yerr(i),a + b*x(i)
         end do
      end if
      close(3)

 666  continue
      call pgend()

      end

