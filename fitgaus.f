
c  subroutine to fit a gaussian to some data
c    uses MRQMIN from Numerical Recipes
c  x, y, yerr are the data and errors
c  n is number of data points
c  pars are the parameters of the gaussian - cont, inty, mean, sigma
c  errpars are the errors on fitted params
c  npars is the number of params
c  ifit is an array that specifies whether to fit a param or not
c  iguess tells if we should make our own initial guesses.
c  if iguess=0, the initial guesses are passed in.
c  if iguess=-1, make our own initial guesses but fit an absorption line.

c      subroutine fitgaus(x,y,yerr,np,pars,errpars,npars,ifit,iguess)
      subroutine fitgaus(x,y,yerr,n,pars,errpars,ifit,iguess,chisq)

      parameter (NPARAMS=4)
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
         if (iguess .eq. -1) then
c  we are fitting an absorption line
            pars(1) = ymax
            pars(3) = x(imin)
c  this is completely bogus
            pars(4) = (xmax-xmin)/20.
c  2.5 is about sqrt(2pi)
            pars(2) = (ymin-ymax)*2.5*pars(4)
         else
c  we are fitting an emission line            
            pars(1) = ymin
            pars(3) = x(imax)
            pars(4) = (xmax-xmin)/20.
            pars(2) = (ymax-ymin)*2.5*pars(4)         
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
      
      do iter = 1,itmax
         oldalam = alamda

         call mrqmin(x,y,yerr,n,pars,ifit,npars,
     $        covar,alpha,npars,chisq,gausfunc,alamda)

         dchisq = (ochisq - chisq) / chisq
         ochisq = chisq
         if((dchisq .lt. conv) .and. (oldalam .gt. alamda)) go to 300
      end do

      badfit = .true.
c      write (*,'("failed to converge -  "$)') 


c  jump here if we converged
 300  continue

      alamda = 0.
      call mrqmin(x,y,yerr,n,pars,ifit,npars,
     $     covar,alpha,npars,chisq,gausfunc,alamda)
c  reduce chisq
      chisq = chisq / real(n-nfit)
      do i=1,nfit
         errpars(i) = sqrt(covar(i,i))
      end do

      return
      end
    

      subroutine gausfunc(x,pars,y,dydpar,npar)

      real pars(npar)
      real dydpar(npar)
      real y0,atot,x0,sig

      root2pi = sqrt(2.0*3.1415927)

      y0 = pars(1)
      atot = pars(2)
      x0 = pars(3)
      sig = pars(4)

      u = (x-x0)**2/2.0/sig**2
      expu = exp(-u)
      r2pisig = root2pi*sig

      y = y0 + atot / r2pisig * expu

      dydpar(1) = 1.0
      dydpar(2) = 1.0 / r2pisig * expu
      dydpar(3) = atot / r2pisig * expu * (x-x0) / sig**2
      dydpar(4) = atot / r2pisig * expu / sig * (2.0*u - 1.0)

      return
      end

