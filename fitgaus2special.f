
c fitgaus2special has the parameters in order inty,sigma,cont,vel
c to allow fitting just inty and sigma.  This shouldn't be necessary,
c but I think there is a problem in mrqfit?

c  subroutine to fit a gaussian to some data
c    uses modified mrqmin, from Numerical Recipes but returns 
c      error flag instead of irritating pausing.
c  x, y, yerr are the data and errors
c  n is number of data points
c  pars are the parameters of the gaussian - cont, inty, mean, sigma
c  errpars are the errors on fitted params
c  npars is the number of params
c  ifit is an array that specifies whether to fit a param or not
c  iguess tells if we should make our own initial guesses.
c  if iguess=2, program is allowed to guess whether em/abs line
c  if iguess=1, program makes its own initial guesses for em line
c  if iguess=0, the initial guesses are passed in.
c  if iguess=-1, make our own initial guesses but fit an absorption line.

c      subroutine fitgaus(x,y,yerr,np,pars,errpars,npars,ifit,iguess)
      subroutine fitgaus2special
     $     (x,y,yerr,n,pars,errpars,ifit,iguess,chisq,ier)

c number of params to fit
      parameter (NPARAMS=4)
c 1-iguess controls whether to fit em or abs line
c 0-try to figure out em or abs from the data
      parameter (IFEMASSUME=1)
      real x(n),y(n),yerr(n)
c      real pars(npars)
c      integer ifit(npars)
      real pars(NPARAMS),errpars(NPARAMS)
      integer ifit(NPARAMS)
      real covar(NPARAMS,NPARAMS),alpha(NPARAMS,NPARAMS)
      logical badfit

      external gausfunc
      
      npars = NPARAMS
      ifitabs = 0

c  generate our own initial guesses?
      if (iguess .ne. 0) then
         xmin = 1.e10
         xmax = -1.e10
         ymin = 1.e10
         ymax = -1.e10
         do i=1,n
            xmin=min(x(i),xmin)
            xmax=max(x(i),xmax)
            if (y(i).gt.ymax) imax = i
            if (y(i).lt.ymin) imin = i
            ymin=min(y(i),ymin)
            ymax=max(y(i),ymax)
         end do
         if (IFEMASSUME .eq. 0 .or. iguess .eq. 2) then
c don't assume whether we fit em or abs
c a guess for the central location should be passed in as pars(3) ?
c or maybe we should just take the midpoint of the data?
c            xguess = (x(1)+x(n))/2.
            if (pars(4) .lt. x(1)) then
               iclose=1
            else if (pars(4) .gt. x(n)) then
               iclose=n
            else
c   find the point ~closest to the guess
               iclose=2
 110           continue
               if (pars(4) .gt. x(iclose) .and. iclose .lt. n) then
                  iclose=iclose+1
                  go to 110
               end if
            end if
            pars(2) = (x(n)-x(1))/20.
c which tercentile are we in?
            ithird = int(((iclose-1)*3.0)/n)+1
            if (ithird .eq. 1) then
               pars(3) = (y(n) + y(n-1) + y(n-2))/3.0
            else if (ithird .eq. 3) then
               pars(3) = (y(1) + y(2) + y(3))/3.0
            else 
               pars(3) = (y(1)+y(2)+y(n-1)+y(n))/4.0
            end if
            pars(1) = (y(iclose) - pars(3))*2.5*pars(2)
         else
            if (iguess .eq. -1) then
c  we are fitting an absorption line
               pars(3) = ymax
               pars(4) = x(imin)
c  this is completely bogus
c               pars(2) = (xmax-xmin)/20.
               pars(2) = (x(n)-x(1))/20.
c  2.5 is about sqrt(2pi)
               pars(1) = (ymin-ymax)*2.5*pars(2)
            else
c  we are fitting an emission line            
               pars(3) = ymin
               pars(4) = x(imax)
c               pars(2) = (xmax-xmin)/20.
               pars(2) = (x(n)-x(1))/20.
               pars(1) = (ymax-ymin)*2.5*pars(2)         
            end if
         end if
      end if


      nfit=0
      do i=1,npars
         if (ifit(i) .ne. 0) nfit = nfit+1
      end do

      alamda = -1.
      itmax = 200
      conv = 1.e-3
      ochisq = 1000.
      
c      do iter = 1,itmax
c         oldalam = alamda
c         call mrqmin(x,y,yerr,n,pars,ifit,npars,
c     $        covar,alpha,npars,chisq,gausfunc,alamda)
c         dchisq = (ochisq - chisq) / chisq
c         ochisq = chisq
c         if((dchisq .lt. conv) .and. (oldalam .gt. alamda)) go to 300
c      end do
      call mrqfit(x,y,yerr,n,pars,ifit,npars,
     $     covar,alpha,npars,chisq,gausfunc,ifconv)

      if (ifconv .eq. 0) then
         badfit = .true.
c      write (*,'("failed to converge -  "$)') 
         ier = 1
      else
c  jump here if we converged
c 300  continue
c      alamda = 0.
c      call mrqmin(x,y,yerr,n,pars,ifit,npars,
c     $     covar,alpha,npars,chisq,gausfunc,alamda)
c  reduce chisq
c      chisq = chisq / real(n-nfit)
c change this
c         do i=1,nfit
         do i=1,npars
            if (ifit(i) .ne. 0) then
               errpars(i) = sqrt(covar(i,i))
            else
               errpars(i) = 0.0
            end if
         end do
         ier = 0
      end if

      return
      end
    

      subroutine gausfunc(x,pars,y,dydpar,npar)

      real pars(npar)
      real dydpar(npar)
      real y0,atot,x0,sig

      root2pi = sqrt(2.0*3.1415927)

      y0 = pars(3)
      atot = pars(1)
      x0 = pars(4)
      sig = pars(2)

      u = (x-x0)**2/2.0/sig**2
      expu = exp(-u)
      r2pisig = root2pi*sig

      y = y0 + atot / r2pisig * expu

      dydpar(3) = 1.0
      dydpar(1) = 1.0 / r2pisig * expu
      dydpar(4) = atot / r2pisig * expu * (x-x0) / sig**2
      dydpar(2) = atot / r2pisig * expu / sig * (2.0*u - 1.0)

      return
      end

