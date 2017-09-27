
c  a front end to several numerical recipes routines
c  for least square fitting - straight line, errors
c  in 0, 1, or 2 coords

c lstsq_bootstrap.f is copyright 2006 Benjamin Weiner and released under
c the GNU General Public License.  You may modify, you may redistribute,
c you must retain the copyright notice.  Please contact me with
c any comments and suggestions.

c See Weiner et al. 2006 (ApJ, in press, astro-ph/------)
c for a discussion of fitting lines to data with intrinsic scatter.

c 4/23/00 - modified to ignore input lines beginning with '#',
c   and to allow fitting only a given region

c 3/18/05 - added bootstrap estimation of the parameters/errors
c   by sampling with replacement and refitting

c 12/20/05 - allows adding the intrinsic scatter in quadrature.
c    following lstsq.f.  Note this is redundant with bootstrapping
c    for error estimation - once you have the correct intrinsic
c    scatter that makes chisq/N=1, the original fit should return
c    an error that is very close to the RMS of the bootstraps.
c    However, the bootstraps are still useful for looking at the
c    covariance of a and b.

c  8/30/06 - added the second plot of the distribution of a vs b
c    as fit to the bootstrap resamples.

c Benjamin Weiner, bjw at astro.umd.edu (Astronomy Dept, University
c of Maryland; Steward Observatory as of Sept 2006)

      program lstsq_bootstrap

      parameter(NMAX=100000)
c      parameter(NMAX_XY=1000)
      parameter(NMAX_XY=10000)
      parameter(NBOOTMAX=10000)

c      real x1(NMAX),x2(NMAX),x3(NMAX)x4(NMAX)
c      real x1err(NMAX),x2err(NMAX),x3err(NMAX)
      real x(NMAX), xerr(NMAX)
      real y(NMAX), yerr(NMAX)
      real xf(NMAX), xferr(NMAX)
      real yf(NMAX), yferr(NMAX)
      real xerrbot(NMAX),xerrtop(NMAX)
      real yerrbot(NMAX),yerrtop(NMAX)
      real xpl(2),ypl(2)

      real xftry(NMAX),yftry(NMAX),xftryerr(NMAX),yftryerr(NMAX)
      real aboot(NBOOTMAX),bboot(NBOOTMAX)

      character fname*40,dline*80

      big = 1.e20

      write(*,*) "0- plain least squares, 1- errors in y"
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
      if (ichoose .ge. 2 .and. npts .gt. NMAX_XY) then
       write(*,*) 
     +   "Fitting w/both X and Y error will not work for >",NMAX_XY
       write(*,*) "points without redimensioning arrays in fitexy.f and"
       write(*,*) "other Numerical Recipes subroutines: chixy.f."
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
         call big_fitexy(xf,yf,nf,xferr,yferr,a,b,siga,sigb,chisq,q)
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

c save the best fit values
      abest = a
      bbest = b
      abesterr = siga
      bbesterr = sigb

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

c initialize random
      idum = -1
 200  continue
      write(*,'("Number of resamples for bootstrap: ",$)')
      read (*,*) nboot
      if (nboot .gt. NBOOTMAX) then
         write(*,*) "max samples ",NBOOTMAX," setting to that"
         nboot = NBOOTMAX
      end if

      asum = 0.0
      asumsq = 0.0
      bsum = 0.0
      bsumsq = 0.0
      amax = -1.e10
      amin =  1.e10
      bmax = -1.e10
      bmin =  1.e10
      open(3,file='lstsq_bootstrap.out',status='unknown')
c For nboot resamples
      do i=1,nboot
c         write(*,*) 'Doing resample ',i
c make random resample
         do j=1,nf
            jrand = int(nf*ran1(idum)) +1
            xftry(j) = xf(jrand)
            yftry(j) = yf(jrand)
            xftryerr(j) = xferr(jrand)
            yftryerr(j) = yferr(jrand)
         end do
c do fit
         if (ichoose .eq. 0) then
c            call fit(x,y,npts,yerr,0,a,b,siga,sigb,chisq,q)
            call fit(xftry,yftry,nf,yftryerr,0,a,b,siga,sigb,chisq,q)
         else if (ichoose .eq. 1) then
c            call fit(x,y,npts,yerr,1,a,b,siga,sigb,chisq,q)
            call fit(xftry,yftry,nf,yftryerr,1,a,b,siga,sigb,chisq,q)
         else if (ichoose .eq. 2) then
c            call fitexy(x,y,npts,xerr,yerr,a,b,siga,sigb,chisq,q)
c            call fitexy_big(xf,yf,nf,xferr,yferr,a,b,siga,sigb,chisq,q)
            call big_fitexy(xftry,yftry,nf,xftryerr,yftryerr,
     +           a,b,siga,sigb,chisq,q)
         end if
c record parameters
         aboot(i) = a
         bboot(i) = b
         asum = asum + a
         asumsq = asumsq + a*a
         bsum = bsum + b
         bsumsq = bsumsq + b*b
         write(3,*) i,a,b
         amax = max(amax,a)
         amin = min(amin,a)
         bmax = max(bmax,b)
         bmin = min(bmin,b)
      end do
      close(3)

c  compute mean and rms of parameters
      amean = asum /nboot
      arms = sqrt(asumsq/nboot - amean**2)
      bmean = bsum /nboot
      brms = sqrt(bsumsq/nboot - bmean**2)

c Plot the bootstrap realizations' values of a and b - 
c useful for seeing covariance

      amarg = 0.2*(amax-amin)
      bmarg = 0.2*(bmax-bmin)
      call pgenv(amin-amarg,amax+amarg,bmin-bmarg,bmax+bmarg,0,1)
      call pglabel("a, intercept","b, slope",
     $     "y = a + bx, best fit and bootstrap realizations")
      call pgpoint(nboot,aboot,bboot,17)
c plot the original best fit point with the error bar from the
c rms of the bootstraps.  It might be better to use the 68% range,
c although they should be similar
      call pgsch(2.0)
      call pgsci(2)
      call pgpt1(abest,bbest,12)
      call pgerrb(5,1,abest,bbest,arms,1.0)
      call pgerrb(6,1,abest,bbest,brms,1.0)
      call pgsch(1.1)
      call pgsci(1)

c compute median and conf levels.  This sorts the aboot and bboot
c arrays so we need to do it last, after plotting and output.
      call sort(nboot,aboot)
      call sort(nboot,bboot)
      if (mod(nboot,2) .eq. 1) then
         indx = int(nboot/2) + 1
         amed = aboot(indx)
         bmed = bboot(indx)
      else
         indx=nboot/2
         amed = (aboot(indx) + aboot(indx+1))/2.0
         bmed = (bboot(indx) + bboot(indx+1))/2.0
      end if
      indx = nint(nboot*(0.5-0.34))
      a68lo = aboot(indx)
      b68lo = bboot(indx)
      indx = nint(nboot*(0.5+0.34))
      a68hi = aboot(indx)
      b68hi = bboot(indx)
      indx = nint(nboot*(0.5-0.477))
      a95lo = aboot(indx)
      b95lo = bboot(indx)
      indx = nint(nboot*(0.5+0.477))
      a95hi = aboot(indx)
      b95hi = bboot(indx)
      
      write(*,*) "Pars from bootstrap with ",nboot," resamples of ",
     +     nf," points"
      write(*,*) "mean, rms a = ",amean," +- ",arms
      write(*,*) "mean, rms b = ",bmean," +- ",brms

      write(*,*) "median and 68%, 95.4% ranges: "
      write(*,*) "a: ",amed,a68lo,a68hi,a95lo,a95hi
      write(*,*) "b: ",bmed,b68lo,b68hi,b95lo,b95hi

      write(*,*) "sigma estimates from 68% and 95%:"
      write(*,*) "a sigma: ",(a68hi-a68lo)/2.0,(a95hi-a95lo)/4.0
      write(*,*) "b sigma: ",(b68hi-b68lo)/2.0,(b95hi-b95lo)/4.0

 666  continue
      call pgend()

      end

