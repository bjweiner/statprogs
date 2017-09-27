
c this program fits a straight line relation in a maximum likelihood way.

c mlsfit.f is derived from linewidthfit2.f
c   allows looping over best value ofscatter

c input is a file with
c x, ex, y, ey

c See Weiner et al. 2006 (ApJ, in press, astro-ph/------)
c for a discussion of fitting lines to data with intrinsic scatter.

c MLSFIT is copyright 2006 Benjamin Weiner and is licensed under the 
c GNU GPL v2.  You may modify it and redistribute it freely; you must
c retain the copyright notice.  If you sell or distribute binaries of it
c (hah) you must also freely provide the source.

c Because this program is derived from a more complex program that
c handled certain obscure and non-Gaussian parts of a specific problem
c (Tully-Fisher relation) there is a fair amount of commented-out
c stuff and some subroutines that aren't used.

c Benjamin Weiner, bjw at astro.umd.edu (Astronomy Dept, University
c of Maryland; Steward Observatory as of Sept 2006)

      program mlsfit

c      parameter (NSIMMAX=1000,NMAX=1000)
c      parameter (NBINMAX=200,NMAGSIM=50)
      parameter(NMAX=12000)
      parameter(NPARMAX=4000,NCONTMAX=41)
      parameter(N3PARMAX=301)
      parameter(IFPOINTS=0,IFCONLABEL=0)
      parameter(IFGRAYSCALE=0)
      parameter(IFMAGERR=1,IFINDIVMAGERR=1)
      parameter(IFPLOTPTS=1)
c recontouring loop?
      parameter(IFRECONT=0)

      real xdat(NMAX),exdat(NMAX),ydat(NMAX),eydat(NMAX)
c      real siglog(NMAX),esiglog(NMAX)
c      real winst0(NMAX),restw(NMAX)
      real probprod(NPARMAX,NPARMAX),problike(NPARMAX,NPARMAX)
      real arraya(NPARMAX),arrayb(NPARMAX)
      real trans(6)
      real clevel(NCONTMAX)
      real xpl(2),ypl(2)

c      real wsigsim(NSIMMAX)
c      integer nsigbin(NBINMAX,NMAGSIM)
c      real rnsigbin(NBINMAX,NMAGSIM)
c      real wbincenter(NBINMAX),wgaus(NBINMAX)
c      real xpl(2),ypl(2)
      character fname*80,toplabel*80,xlabel*80,clabel*20,cline*80

c      external gasdev

      root2pi = sqrt(2.0*3.14159265)
      ln10 = 2.302585
      ncont = 21
      probfloor = 1.0e-40
c can use zpx to make fit around a zeropoint, e.g. if the range
c of x is -18 to -22, set zeropoint of x to -20 to fit as function
c of (x- -20) and reduce covariance between intercept and slope.
      zpx = 0.0
      write(*,'("input a zeropoint for the x range: ",$)')
      read(*,*) zpx

      call pgbegin(0,'?',2,2)
      call pgsch(1.5)
      call pgscf(2)

 100  continue
      write(*,'("file with x,err_x,y,err_y: ",$)')
      read (*,'(a)') fname
      open(2,file=fname,status='old',err=100)
      nobj=1
 110  continue
c if each point has its own x err
      if (IFMAGERR .eq. 1 .and. IFINDIVMAGERR .eq. 1) then
            read(2,*,err=120,end=120) xdat(nobj),exdat(nobj),ydat(nobj),
     $           eydat(nobj)
      else
         read(2,*,err=120,end=120) xdat(nobj),ydat(nobj),
     $        eydat(nobj)
      end if
      nobj=nobj+1
      if (nobj .gt. NMAX) write(*,*) "error, nobj > nmax =",NMAX
      go to 110
 120  continue
      close(2)
      nobj=nobj-1
      write (*,*) "Read ",nobj," data points"

c hack to set an x err if we didn't have one for each point
      if (IFMAGERR .eq. 1 .and. IFINDIVMAGERR .ne. 1) then
         write(*,'("x error for all points: ",$)')
         read (*,*) exconst
         do i=1,nobj
            exdat(i) = exconst
         end do
      end if

 200  continue
c      write(*,'("a,b for y = a + b x: ",$)')
c      read(*,*) a,b

c choose which parameters to vary/loop over when computing likelihood
c use x and y to switch the ones that can be fixed
c For mlsfit, don't fit z evolution, so this is redundant
c      write(*,990)
c 990  format("Fix z evolution (0), TF slope (1), neither (2) [quit]: ",
c     $     $)
c      read (*,'(a)') cline
c      if (cline(1:3) .eq. '   ') go to 666
c      read(cline,*) ifixslope
      write(*,'("test y = a + b*(x- ",f5.2,")")') zpx
      write(*,'("intercept amin, amax, da: ",$)')
      read(*,*) amin,amax,da
      na = int((amax-amin)/da)+1
c      nx = na
      write(*,'("slope bmin, bmax, db: ",$)')
      read(*,*) bmin,bmax,db
      nb = int((bmax-bmin)/db)+1
      ny = nb
      ymin = bmin
      ymax = bmax
      dy = db
c      write(*,'("TF rms scatter in logsig: ",$)')
c      read (*,*) tfsig
      write(*,'("rms intrinsic scatter in y, dmin, dmax, dd: ",$)')
      read(*,*) dmin,dmax,dd
      nd = int((dmax-dmin)/dd)+1
c The minimum error fudge isn't really needed when you are adding
c intrinsic scatter.  Otherwise, it is useful to keep points with
c unrealistically small errors from dominating the fit.
c      write(*,'("Minimum error fudge in err(y): ",$)')
c      read (*,*) eyfudge
      eyfudge = 0.0

c if disp in A was input, compute the logsig.
c this is just for plotting the relation.    
c      do i=1,nobj
c         siglog(i) = wobs(i)
c         esiglog(i) = ewobs(i)
c      end do

c Note that y is for plotting the fit params, so they
c are really a and b, not x and y the input data points
c      na = int((amax-amin)/da)+1
c      nx = int((xmax-xmin)/dx)+1
c      nb = int((bmax-ymin)/dy)+1
         
      if (na .gt. NPARMAX .or. nb .gt. NPARMAX
     $     .or. nd .gt. NPARMAX) then
         write(*,*) "Error, na,nb,nd,NPARMAX : ",na,nb,nd,NPARMAX
      else
         write(*,*) "na, nb, nd  = ",na,nb,nd
      end if

      pmaxoverall = -1.0e20

c loop over scatter values
      do id = 1,nd
         d = dmin + (id-1)*dd
         open(3,file="mlsfit.out",status='unknown')
         plikemin = 1.0e10
         plikemax = -1.0e20
c loop over slope
         do ib=1,nb
            b = bmin + (ib-1)*db
c loop over intercept
            do ia = 1,na
               a = amin + (ia-1)*da
               probprod(ia,ib) = 1.0
               problike(ia,ib) = 0.0
               chisq = 0.0
c loop over objects
               do i=1,nobj

c Some of this stuff is irrelevant to the real MLS-fit method
c that uses a convolution in a subroutine, e.g. calcprob1()

c What to do if working in Angstroms
c               if (IFKMS .eq. 0) then
c               if (IFMAGERR .eq. 0) then
c what to do if working in km/s
c               else
c predicted y given a,b,x
               ypred = a + b*(xdat(i)-zpx)
c scatter in y
               scatt = d
c fudged error bar on observed y
               yerrtotsq = eydat(i)**2 + eyfudge**2
c probability of getting the observed dispersion
c we have two gaussians, one from the estimated error of the 
c observation, the other from the scatter in TF
c  just convolve the two
                  yscatttotsq = yerrtotsq + scatt**2
                  ydiff = ydat(i) - ypred
                  devsq = ydiff**2/yscatttotsq
c note, devsq just calculates the deviation from the y-displacement.
c it doesn't take into account the x.  there is a more complex
c method of finding the orthogonal deviation:
c Try calculating orthogonal deviation. 
c distance from point to line is D^2 = (xp-x(i))**2 + (yp-y(i))**2
c where xp,yp is the closest point on the line, which is given by
c minimizing D^2 wrt x.  solving d(D^2)/dx=0 gives 
c xp = (xi + b*yi - ab)/(1 +b**2), then calculate deviation
c                  xp = (x(i) + b*y(i) - a*b)/(1. + b**2)
c                  yp = a + b*xp
c                  devsq = (xp-x(i))**2/ex(i)**2 + 
c     +                 (yp-y(i))**2/(ey(i)**2 + dsq)
c However this still isn't really right, because we should calculate
c the min distance in terms x and y scaled by the error bars and 
c intrinsic scatter.  Apportion the intrinsic scatter S (=d) into x and y
c components Sx and Sy: Sy/Sx = 1/b, so (1+b**2)Sy**2 = S**2
c then deviation D^2 = (xp-xi)**2/(ex**2 + Sx**2) + (yp-yi)**2/(ey**2 + Sy**2)
c and we need to minimize D^2 to find xp,yp
                  scattysq = dsq / (1. + b**2)
                  scattxsq = dsq - scattysq
                  wxsq = 1./(exdat(i)**2 + scattxsq)
                  wysq = 1./(yscatttotsq + scattysq)
                  xp = (wxsq*xdat(i) + wysq*b*ydat(i) - wysq*a) / 
     +                 (wxsq + wysq * b**2)
                  yp = a + b*xp
                 orthdevsq = wxsq*(xp-xdat(i))**2 + wysq*(yp-ydat(i))**2
               if (IFMAGERR .eq. 0) then
                  prob = 1.0/root2pi/sqrt(yscatttotsq) * 
     $                 exp(-devsq/2.0)
               else
c or calc. max likelihood in subroutine. This is better and how
c to take care of mag errors
c                  aadj = a + c*z(i)
                  aadj = a
                  badj = b
                  call calcprob1(xdat(i),ydat(i),exdat(i),
     $                 sqrt(yerrtotsq),zpx,aadj,badj,
     $                 scatt,prob)
               end if

c end the (IFKMS .eq. 0) clause
c               end if

c to avoid underflow
               if (prob .lt. probfloor) prob = probfloor
               probprod(ia,ib) = probprod(ia,ib) * prob
               problike(ia,ib) = problike(ia,ib) + log(prob)
               chisq = chisq + devsq
               orthchisq = orthchisq + orthdevsq
c end looping over objects
            end do
            
            plikemin=min(plikemin,problike(ia,ib))
            if (problike(ia,ib) .gt. plikemax) then
               plikemax = problike(ia,ib)
               iamax = ia
               ibmax = ib
               pamax = a
               pbmax = b
               chimin = chisq
               orthchimin = orthchisq
            end if
            if (problike(ia,ib) .gt. pmaxoverall) then
               pmaxoverall = problike(ia,ib)
               iaall = ia
               iball = ib
               aall = a
               ball = b
c               call = c
               dall = d
               chiminall = chisq
               orthchiminall = orthchisq
            end if
            write(3,1030) ia,ib,a,b,d,probprod(ia,ib),problike(ia,ib),
     $           nobj,chisq,orthchisq
 1030      format(i5,1x,i5,3(2x,f6.3),2(2x,1pe10.3),2x,i5,2(2x,1pe10.3))
c end looping over slopes and intercepts
            end do
         end do

         close(3)

c make a plot, for this particular value of the scatter
c (still inside the scatter loop)
      trans(1) = amin-da
      trans(2) = da
      trans(3) = 0.0
      trans(4) = bmin-db
      trans(5) = 0.0
      trans(6) = db

c linear: iflog=0 or log: iflog=1 contours?
      iflog = 0
      if (iflog .eq. 0) then
         dcont = (plikemax-plikemin) / ncont
         do i=1,ncont
            clevel(i) = plikemin + i*dcont
         end do
      else
         dcont = (log(plikemax-plikemin)) /ncont
         do i=1,ncont
            clevel(i) = plikemin + exp(i*dcont)
         end do
      end if

      write(*,1040) plikemax,iamax,ibmax,pamax,pbmax,d,nobj,
     $     chimin,orthchimin
 1040 format("Max log likelihood: ",1pe10.3," at a,b,d",2(2x,i4),
     $     3(2x,1pe10.3),", N= ",i5," chisq= ",1pe10.3,2x,1pe10.3)
c      write(*,*) plikemax,iamax,ixmax,pamax,pxmax

c compute error bars on the values (a,x) that are being contoured
c over by finding the 68% probability ranges.  assume that the
c probability outside the grid range is small.  actually computing
c the +-34% range from the peak value, not from the expectation 
c value (sigh).  Need to sample the par space finely or you
c can have a situation where there isn't 34% of probability on
c one side of the peak.

c loop through and accumulate total probability, and prob
c as function of a, and x (x=b)

c Try finding the 2 sigma error range, because 1 sigma is sometimes
c within 1 grid cell?

c      write(*,*) "call finderr,",na,nx,iamax,ixmax,NPARMAX
      pthresh = 0.34
c      pthresh = 0.477
      call finderrrange(problike,na,nb,NPARMAX,NPARMAX,iamax,ibmax,
     +     pthresh,erra1index,erra2index,errb1index,errb2index)
      erra1 = amin + (erra1index-1)*da
      erra2 = amin + (erra2index-1)*da
      errb1 = bmin + (errb1index-1)*db
      errb2 = bmin + (errb2index-1)*db

      write(*,'("Error range including ",f6.4)') 2.0*pthresh
      write(*,1042) pamax,pamax-erra1,erra2-pamax,
     +     pbmax,pbmax-errb1,errb2-pbmax
 1042 format("Best a,b and err ",1pe10.3," - ",1pe10.3," + ",1pe10.3,
     +     3x,1pe10.3," - ",1pe10.3," + ",1pe10.3)

 500  continue
      write(toplabel,1050) pamax,pbmax,zpx,d
 1050 format("Best y= ",f5.2," + ",f6.3,"(x- ",f6.2,
     $     "), scatt",f6.3)
c 1060 format("Best logsig = ",f5.2," + ",f6.3,"(M_B- ",f5.1,
c     $     ") + ",f6.3," z")
      if (IFGRAYSCALE .eq. 1) then
         call pgenv(amin-da,amax+da,bmin-db,bmax+db,0,0)
         call pglabel("a","b",toplabel)
         call pggray(problike,NPARMAX,NPARMAX,1,na,1,nb,
     $     clevel(ncont),clevel(1),trans)
         call pgqci(ioldci)
         call pgsci(2)
         call pgpt1(pamax,pbmax,17)
         call pgsci(ioldci)
      end if

      call pgenv(amin-da,amax+da,bmin-db,bmax+db,0,0)

      if (IFPOINTS .eq. 1) then
c show points where the likelihood was calculated
         do ia=1,na
            a = amin + (ia-1)*da
            do ib = 1,nb
               arraya(ib) = a 
               arrayb(ib) = bmin + (ib-1)*db
            end do
            call pgqch(cheight)
            call pgsch(0.5)
            call pgpt(nx,arraya,arrayb,17)
            call pgsch(cheight)
         end do
      end if

c contour plot the likelihood
      if (iflog .eq. 0) then
         write(*,*) "linear contours ",clevel(1),clevel(ncont),dcont
      else
         write(*,*) "log contours ",clevel(1),clevel(ncont),dcont
      end if
      write(xlabel,1070) clevel(1),clevel(ncont),dcont
 1070 format("a, contours ",1pe9.2," to ",1pe9.2," by ",1pe9.2)
      call pglabel(xlabel,"b",toplabel)
c      call pgcont(problike,NPARMAX,NPARMAX,1,na,1,nx,
c     $        clevel,ncont,trans)
      call pgcont(problike,NPARMAX,NPARMAX,1,na,1,nb,
     $     clevel,-ncont,trans)
      if (IFCONLABEL .eq. 1) then
         call pgqch(cheight)
         call pgsch(0.9)
         do i=1,ncont
            write(clabel,'(1pe8.1)') clevel(i)
            call pgconl(problike,NPARMAX,NPARMAX,1,na,1,nb,
     $           clevel(i),trans,clabel,70,10)
         end do
         call pgsch(cheight)
      end if
      call pgqci(ioldci)
      call pgsci(2)
      call pgpt1(pamax,pbmax,17)
      call pgsci(ioldci)


c plot lines for the error ranges.
      xpl(1) = amin-da
      xpl(2) = amax+da
      ypl(1) = errb1
      ypl(2) = errb1
      call pgline(2,xpl,ypl)
      ypl(1) = errb2
      ypl(2) = errb2
      call pgline(2,xpl,ypl)
      xpl(1) = erra1
      xpl(2) = erra1
      ypl(1) = bmin-db
      ypl(2) = bmax+db
      call pgline(2,xpl,ypl)
      xpl(1) = erra2
      xpl(2) = erra2
      call pgline(2,xpl,ypl)

c plot points and best fit
      if (IFPLOTPTS .ne. 0) then
         xpmin = 1.0e10
         xpmax = -1.0e10
         ypmin = 1.0e10
         ypmax = -1.0e10
         do i=1,nobj
            xpmin = min(xpmin,xdat(i))
            xpmax = max(xpmax,xdat(i))
            ypmin = min(ypmin,ydat(i))
            ypmax = max(ypmax,ydat(i))
         end do
         xmarg = 0.05 * (xpmax-xpmin)
         ymarg = 0.05 * (ypmax-ypmin)
         write(*,*) xpmin,xpmax, ypmin,ypmax
         call pgenv(xpmin-xmarg,xpmax+xmarg,ypmin-ymarg,ypmax+ymarg,0,0)
c         call pgenv(4.0,-3.0,0.8,2.7,0,0)
         call pglabel("x","y","data and best fit")
         call pgpt(nobj,xdat,ydat,17)
c         call pgerry(nobj,bmag,siglog,esiglog,0.0)
c         if (IFMAGERR .ne. 0 .and. IFKMS .ne. 0)
c     $        call pgerrx(nobj,bmag,siglog,ebmag,0.0)
c draw lines for best relation and its 1-sigma scatter
         xpl(1) = xpmax+xmarg
         xpl(2) = xpmin-xmarg
c         if (ifixslope .eq. 0) then
            ypl(1) = pamax + pbmax * xpl(1)
            ypl(2) = pamax + pbmax * xpl(2)
c         else
c            ypl(1) = pamax + b * xpl(1)
c            ypl(2) = pamax + b * xpl(2)
c         end if
         call pgline(2,xpl,ypl)
         ypl(1) = ypl(1) + d
         ypl(2) = ypl(2) + d
         call pgline(2,xpl,ypl)
         ypl(1) = ypl(1) - 2.0*d
         ypl(2) = ypl(2) - 2.0*d
         call pgline(2,xpl,ypl)
      end if

      
      if (IFRECONT .ne. 0) then
c go through the recontouring loop

 510     continue
         write(*,'("New contmin,contmax,ncont (- for log) [next]: ",$)')
         read(*,'(a)') cline
         if (cline(1:4) .ne. '    ') then
            read(cline,*,err=510) contmin,contmax,ncont
            if (ncont .gt. 0) then
               iflog = 0
               ncont=min(ncont,NCONTMAX)
               dcont = (contmax-contmin)/(ncont-1.0)
               do i=1,ncont
                  clevel(i) = contmin + (i-1.0)*dcont
               end do
               go to 500
            else if (ncont .lt. 0) then
               iflog = 1
               ncont = -ncont
               ncont=min(ncont,NCONTMAX)
               contmin = max(contmin,1.0-e6)
               contminlog = log(contmin)
               dcont = (log(contmax)-contminlog)/(ncont-1.0)
               do i=1,ncont
                  clevel(i) = exp(contminlog + (i-1.0)*dcont)
               end do
               go to 500
            end if
         end if
      end if      


c end the do iy loop
c      end do
c end the do id loop - looping over scatter
      end do
c This ends the main loops.

      write(*,*) "Overall minimum: "
c      write (*,*) iaall,iball,aall,ball,dall,
c     $     chiall,nobj,chiminall
      write (*,*) aall,ball,dall,nobj,chiminall

c      go to 200

 666  continue
      call pgend()

      end

c find error ranges by checking a threshold eg. 68/2 percent
c compute error bars on the values (a,x) that are being contoured
c over by finding the 68% probability ranges.  assume that the
c probability outside the grid range is small.  actually computing
c the +-34% range from the peak value, not from the expectation 
c value (sigh).  This is a bit tedious...
c returns fractional INDEXES, not actual a,b values

c I had to make some variables double precision because the
c probability = exponential can produce some underflows

      subroutine finderrrange(probarr,na,nb,naphys,nbphys,
     +     iamax,ibmax,pthresh,ea1,ea2,eb1,eb2)

c is our array probability, or likelihood = log P?
      parameter(IFLIKE=1)

      real probarr(naphys,nbphys)
c      real proba(na), probb(nb)
      double precision proba(na), probb(nb)
      double precision prob1,sumprobtot

c because the probability can often be quite small, scale it up
      pgridmax=-1.e10
      do ia=1,na
         do ib=1,nb
            pgridmax = max(pgridmax,probarr(ia,ib))
         end do
      end do

      sumprobtot = 0.d0
      do ia=1,na
         proba(ia) = 0.d0
      end do
      do ib=1,nb
         probb(ib) = 0.d0
         do ia=1,na
            if (IFLIKE .eq. 1) then
c               prob1 = exp(dble(probarr(ia,ib)))
               prob1 = exp(dble(probarr(ia,ib)-pgridmax))
            else
c               prob1 = dble(probarr(ia,ib))
               prob1 = dble(probarr(ia,ib) / pgridmax)
            end if
            sumprobtot = sumprobtot + prob1
            proba(ia) = proba(ia) + prob1
            probb(ib) = probb(ib) + prob1
         end do
      end do
c  normalize
c      suma = 0.0
c      sumb = 0.0
      open(11,file='mlsfit.err.out',status='unknown')
      do ia=1,na
         proba(ia) = proba(ia)/sumprobtot
c         suma = suma + proba(ia)
         write(11,*) "ia ",ia,proba(ia)
      end do
      do ib=1,nb
         probb(ib) = probb(ib)/sumprobtot
c         sumb = sumb + probb(ib)
         write(11,*) "ib ",ib,probb(ib)
      end do
      close(11)
      write(*,*) "finderr: sumtot, P(amax), P(bmax): ",sumprobtot,
     +     proba(iamax),probb(ibmax)
c      write(*,*) suma,sumb

c  start at max point iamax,ixmax and go up looking to accumulate 34 percent
c      pthresh = 0.34
      oldpsum = 0.0
      psum = proba(iamax) / 2.0
      ia = iamax
 610  continue
c      write (*,*) ia,oldpsum,psum,pthresh
      if (psum .gt. pthresh) then
         frac = (pthresh-oldpsum)/(psum-oldpsum)
         ea2 = (ia-0.5) + frac
      else
         ia = ia+1
         if (ia .le. na) then
            oldpsum = psum
            psum = psum + proba(ia)
            go to 610
         else
            ea2 = real(na)
         end if
      end if
c go down in a
      oldpsum = 0.0
      psum = proba(iamax)/2.0
      ia = iamax
 615  continue
c      write (*,*) ia,oldpsum,psum,pthresh
      if (psum .gt. pthresh) then
         frac = (pthresh-oldpsum)/(psum-oldpsum)
         ea1 = (ia+0.5) - frac
      else
         ia = ia-1
         if (ia .ge. 1) then
            oldpsum = psum
            psum = psum + proba(ia)
            go to 615
         else
            ea1 = 1.0
         end if
      end if
c start at max point iamax,ixmax and go up looking to accumulate 34 percent
      oldpsum = 0.0
      psum = probb(ibmax)/2.0
      ib = ibmax
 620  continue
c      write (*,*) ib,oldpsum,psum,pthresh
      if (psum .gt. pthresh) then
         frac = (pthresh-oldpsum)/(psum-oldpsum)
         eb2 = (ib-0.5) + frac
      else
         ib = ib+1
         if (ib .le. nb) then
            oldpsum = psum
            psum = psum + probb(ib)
            go to 620
         else
            eb2 = real(nb)
         end if
      end if
c go down in b
      oldpsum = 0.0
      psum = probb(ibmax)/2.0
      ib = ibmax
 625  continue
c      write (*,*) ib,oldpsum,psum,pthresh
      if (psum .gt. pthresh) then
         frac = (pthresh-oldpsum)/(psum-oldpsum)
         eb1 = (ib+0.5) - frac
      else
         ib = ib-1
         if (ib .ge. 1) then
            oldpsum = psum
            psum = psum + probb(ib)
            go to 625
         else
            eb1 = 1.0
         end if
      end if

      return 
      end


c in order to do mag errors and sigma errors and TF scatter all
c at the same time, need to find the probability a point has of
c occurring given a TF relation ypred = a + bx
c with scatter tfsig
c and a data point xi, yi, with errors ex, ey
c this means convolving the TF distribution (probability contours
c parallel to TF line) with prob. ellipses around the data point,
c I think.

c This was the wrong way to do that:
c define deviations as a fn of location in (x,y):
c dtfsq = (y - ypred)**2 / tfsig**2 = (y - a - bx)**2 / tfsig**2
c dobssq = (x - xi)**2 / ex**2 + (y - yi)**2 / ey**2
c probability P = 1/sqrt(2pi)/something * exp(-(dtfsq + dobsssq)/2)
c max in probability is when S = dtfsq + dobssq is minimized
c 0 = dS/dx = ex**2 (-ab + b**2 * x - by) + tfsig**2 (-xi + x)
c 0 = dS/dy = ey**2 (-a  - bx       + y)  + tfsig**2 (-yi + y)
c so
c 1) 0 = ex**2 (-ab) + tfsig**2 (-xi) + x(ex**2 b**2 + tfsig**2) + y(-ex**2 b)
c 2) 0 = ey**2 (-a)  + tfsig**2 (-yi) + x(-ey**2 b)              + y(ey**2 + tfsig**2)
c    let K = (ex**2 b**2 + tfsig**2) / (ey**2 b), and compute 1+2*K to elim. x
c  0 = ex**2(-ab) + tfsig**2(-xi) + ey**2(-a K) +tfsig**2(-yi K) + 
c        y (-ex**2 b + K (ey**2 + tfsig**2))
c  so y = (ex**2(-ab) + ey**2(-a K) + tfsig**2(-xi - yi * K)) /
c            (ex**2 b + K * (ey**2 + tfsig**2))
c Knowing y, then
c   x = (y (ex**2 + tfsig**2) + ey**2(-a) + tfsig**2(-yi)) / (ey**2 b)
c Then can calculate dtfsq, dobssq, and P, assuming we can
c normalize P properly ...
c  It seems like a use of Bayes's theorem might simplify this mess?

c actually I think this is the wrong approach.  I need to find the
c probability by convolving the distributions: the error ellipse
c with the TF relation distrib.  fortunately, since the elliptical
c gaussian is separable, the problem should be _much_ easier.

c Actually this is the wrong way to go about it.  Rather, one wants to
c convolve the probability distribs of the TF and the observation to
c find a total prob
c  dP_tf  = 1/(root2pi tfsig) exp( -(y-a-bx)**2/ (2 tfsig**2) )
c  dP_obs = 1/(2pi ex ey) exp( -(x-xi)**2/(2 ex**2) - (y-yi)**2/(2 ey**2) )
c  P = integral( dP_tf dP_obs dx dy)
c  P = 1 / (2pi^(3/2) tfsig ex ey) *
c     integral( -(y-a-bx)**2/ (2 tfsig**2) -(x-xi)**2/(2 ex**2) - (y-yi)**2/(2 ey**2) )
c this has a cross term; might be able to approximate it by 
c evaluating at a small number of x samples

c calculate probability by doing convolution of:
c gaussian error ellipse at x,y with sigmas ex, ey
c and TF line y = f + g(x-zp) with scatter yscatt
c note that f here can = a + c*z(i) in main program
c slope g could also change with z

c This convolution is kind of cheesy and will work better
c when the slope of the best fit line is relatively shallow.
c if the slope b(or g) is very steep, so that the line changes a lot
c per x-sample, the convolution won't be good, but this is a sign
c that the user needs to change variables.

      subroutine calcprob1(x,y,ex,ey,zpx,f,g,yscatt,prob)
c number of sigma in x to go to and number of samples per sigma
c      parameter(NSIG=4,NSAMP=4)
      parameter(NSIG=4,NSAMP=4)

      root2 = sqrt(2.0)
      root2pi = sqrt(2.0*3.14159265)
      nx = 2*NSIG*NSAMP
      dx = ex / NSAMP
      x1 = x - (NSIG*NSAMP-0.5)*dx
      prob = 0.0
      xnormtot = 0.0
      do i=1,nx
         xcen = x1 + (i-1.0)*dx
         xmin = xcen - dx/2.0
         xmax = xmin + dx
c erf is defined as integral exp(-x^2), so need to add a sqrt(2) factor???
c         xnorm = abs( erf((xmin-x)/ex) - erf((xmax-x)/ex) ) / 2.0
         xnorm = abs( erf((xmin-x)/ex/root2) - 
     $        erf((xmax-x)/ex/root2) ) / 2.0
c         write (*,*) erf((xmin-x)/ex),erf((xmax-x)/ex)
         ysigtotsq = ey**2 + yscatt**2
         ydiff = y - (f + g*(xcen-zpx))
         p1 = xnorm * 1./root2pi/sqrt(ysigtotsq) * 
     $        exp(-ydiff**2/2.0/ysigtotsq)
         xnormtot = xnormtot + xnorm
         prob = prob + p1
      end do
c xnormtot should come out slightly less than 1.
      prob = prob / xnormtot
c      write(*,*) f,g,x,y,xnormtot,prob,log(prob)
      return
      end

cccccccccccccccccccccccccccccc
c  calc prob with mag error and doing it in angstrom space
c  much of this is similar to calcprob1
      subroutine calcprob2ang(x,ex,wdisp,ewdisp,zpx,f,g,winst,
     $     winstsig,waveobs,tfscatt,prob)

      parameter(NSIG=4,NSAMP=4)

      root2 = sqrt(2.0)
      root2pi = sqrt(2.0*3.14159265)
      if (ex .lt. 1.e-4) then
c mag err negligible. x-input is on bin center
c with 1 bin, xnorm is about 1 sigma and should cancel out.
         nx = 1
         dx = ex
         x1 = x
      else
c even number of bins, so x-input is on bin edge
         nx = 2*NSIG*NSAMP
         dx = ex / NSAMP
         x1 = x - (NSIG*NSAMP-0.5)*dx
      end if
      prob = 0.0
      xnormtot = 0.0
      do i=1,nx
         xcen = x1 + (i-1.0)*dx
         xmin = xcen - dx/2.0
         xmax = xmin + dx
c erf is defined as integral exp(-x^2), so need to add a sqrt(2) factor???
c         xnorm = abs( erf((xmin-x)/ex) - erf((xmax-x)/ex) ) / 2.0
         xnorm = abs( erf((xmin-x)/ex/root2) - 
     $        erf((xmax-x)/ex/root2) ) / 2.0
c predicted wdisp from TF
         ypredlogkms = f + g*(xcen-zpx)
         ypredkms = 10**ypredlogkms
         ypredang = ypredkms * waveobs/3.0e5
         predwdispsq = ypredang**2 + winst**2
         predwdisp = sqrt(predwdispsq)
         wdiffang = predwdisp - wdisp
c tfscatt in A space is effectively an error on y and
c transforms as such, so need to multiply by d(predwdisp)/dy
c here we're neglecting, unfortunately, that the scatter, if normal
c in log kms, is not normal in wdisp.
c It might be a good idea to fix this.
c         ddispsqdy = 2.0*ypredang *waveobs/3.0e5 * ypredkms * 2.3026
         ddispsqdy = 2.0*ypredang**2 * 2.3026
c tfscatt in A = d(disp)/dy * tfscatt in log sigma
         scattang = (0.5/predwdisp * ddispsqdy) * tfscatt
         wsigtotsq = (ewdisp**2 + scattang**2 + winstsig**2)
         p1 = xnorm * 1./root2pi/sqrt(wsigtotsq) * 
     $        exp(-wdiffang**2/2.0/wsigtotsq)
         xnormtot = xnormtot + xnorm
         prob = prob + p1
      end do
c xnormtot should come out slightly less than 1, except for n=1 case
      prob = prob / xnormtot
c      write(*,*) f,g,x,xnormtot,prob,log(prob),p1errsum,p1scattsum
      return
      end

cccccccccccccccccccccccccccccc
c  calc prob with mag error and doing it in angstrom space
c  much of this is similar to calcprob1
c  calcprob3ang is like calcprob2ang but tries to treat the
c  non-gaussian tf scatter properly, which requires another
c  integral
      subroutine calcprob3ang(x,ex,wdisp,ewdisp,zpx,f,g,winst,
     $     winstsig,waveobs,tfscatt,prob)

      parameter(NSIG=4,NSAMP=4)
      parameter(NWSIG=4,NWSAMP=2)

      root2 = sqrt(2.0)
      root2pi = sqrt(2.0*3.14159265)
      winstsq = winst**2
      tfscattsq = tfscatt**2
      if (ex .lt. 1.e-4) then
c mag err negligible. x-input is on bin center
c with 1 bin, xnorm is about 1 sigma and should cancel out.
         nx = 1
         dx = ex
         x1 = x
      else
c even number of bins, so x-input is on bin edge
         nx = 2*NSIG*NSAMP
         dx = ex / NSAMP
         x1 = x - (NSIG*NSAMP-0.5)*dx
      end if
      prob = 0.0
      xnormtot = 0.0
      do i=1,nx
         xcen = x1 + (i-1.0)*dx
         xmin = xcen - dx/2.0
         xmax = xmin + dx
c erf is defined as integral exp(-x^2), so need to add a sqrt(2) factor???
c         xnorm = abs( erf((xmin-x)/ex) - erf((xmax-x)/ex) ) / 2.0
         xnorm = abs( erf((xmin-x)/ex/root2) - 
     $        erf((xmax-x)/ex/root2) ) / 2.0
c predicted wdisp from TF
         ypredlogkms = f + g*(xcen-zpx)
         ypredkms = 10**ypredlogkms
         ypredang = ypredkms * waveobs/3.0e5
         predwdispsq = ypredang**2 + winst**2
         predwdisp = sqrt(predwdispsq)
         wdiffang = predwdisp - wdisp
c tfscatt in A space is effectively an error on y and
c transforms as such, so need to multiply by d(predwdisp)/dy
c here we're neglecting, unfortunately, that the scatter, if normal
c in log kms, is not normal in wdisp.
c It might be a good idea to fix this.
c         ddispsqdy = 2.0*ypredang *waveobs/3.0e5 * ypredkms * 2.3026
         ddispsqdy = 2.0*ypredang**2 * 2.3026
c tfscatt in A = d(disp)/dy * tfscatt in log sigma
         scattang = (0.5/predwdisp * ddispsqdy) * tfscatt
         wsigtotsq = (ewdisp**2 + scattang**2 + winstsig**2)
         p1 = xnorm * 1./root2pi/sqrt(wsigtotsq) * 
     $        exp(-wdiffang**2/2.0/wsigtotsq)
c here we do a loop to integrate the convolution of 
c the error and TF scatter distributions in wavelength space.
c figure the range of the loop by the (min,max) of the NWSIG sigma
c ranges in error and scatter, and sample by NWSAMP per 
c whichever sigma is smaller
         w1 = min(wdisp-NWSIG*ewdisp,predwdisp-NWSIG*scattang)
         w2 = max(wdisp+NWSIG*ewdisp,predwdisp+NWSIG*scattang)
         dw = min(ewdisp,scattang)/NWSAMP
         nw = int((w2-w1)/dw)+1
c formally, can't have a dispersion lower than winst so the 
c convolution is reliant on the part of p1err that is above winst.
c if none, then the observation is theoretically impossible.
c in practice, i want to force it to have at least some
c probability above winst.  Try: if the wdisp is more than 1 sigma (ewdisp)
c below winst, force the error larger to bring it to 1 sigma.
c use a tmp variable to avoid redefining ewdisp
         if (wdisp+ewdisp .lt. winst) then
            ewdisptmp = winst-wdisp
         else
            ewdisptmp = ewdisp
         end if
c allow for variations in siginst
         ewdisptmpsq = ewdisptmp**2 + winstsig**2
         ewdisptmp=sqrt(ewdisptmpsq)
         p1 = 0.0
         p1errsum = 0.0
         p1scattsum = 0.0
c I need to think about these in terms of probability density over w.
c p1err has dimensions 1/w
c p1scatt has dimensions 1/log(kms), I need to multiply 
c by d(log y)/d(predwdisp).  If I get this right, p1errsum and p1scattsum,
c which are only for debugging, should come out around 0.999 (4 sigma range)
         do j=1,nw
            wcen = w1 + (j-1)*dw
            p1err = 1./root2pi/ewdisptmp * 
     $           exp(-(wdisp-wcen)**2/2.0/ewdisptmpsq)
            p1err = p1err * dw
            p1errsum = p1errsum + p1err
c compute the y (log vel) that corresponds to wcen to get tfscatt
            if (wcen .le. winst) then
c formally impossible, just set prob very small
               p1scatt = 1.e-6
            else
               ycenang = sqrt(wcen**2 - winstsq)
               ycenkms = 3.0e5/waveobs * ycenang
               ycenlogkms = log10(ycenkms)
               p1scatt = 1./root2pi/tfscatt * 
     $              exp(-(ypredlogkms-ycenlogkms)**2/2.0/tfscattsq)
               dlogydw = 1./ycenkms/2.3026 * 3.0e5/waveobs * 
     $              wcen/ycenang
               p1scatt = p1scatt * dw * dlogydw
               p1scattsum = p1scattsum + p1scatt
            end if
            p1 = p1 + p1err*p1scatt
         end do
         xnormtot = xnormtot + xnorm
         prob = prob + xnorm*p1
      end do
c xnormtot should come out slightly less than 1, except for n=1 case
      prob = prob / xnormtot
c      write(*,*) f,g,x,xnormtot,prob,log(prob),p1errsum,p1scattsum
      return
      end


         
cccccccccccccccccccccccccccccc
c I now think the probability can be computed without doing a
c convolution
c It seems that I should be able to convolve the y-error
c with the y-scatter by adding in quadrature, and to
c convert the x-error into an additional component of
c y-error by a factor in the slope.

c x,y,ex,ey are the observed point and errors.
c the relation is y = f + g*x with scatter yscatt
      subroutine calcprob2(x,y,ex,ey,zpx,f,g,yscatt,prob)
      parameter (IFDEBUG=0)

      root2pi = sqrt(2.0*3.14159265)

c converting x error into y: if slope g is shallow, x error changes
c probability of finding this point little.  if slope is steep,
c and the point is (for example) 3*yerr off in y but only 1*xerr
c off in x, then it converts to a reasonable y-error
      yerrfromx = abs(g) * ex
      ydiff = y - (f + g*(x-zpx))
      ysigtotsq = yscatt**2 + ey**2 + yerrfromx**2
c not entirely confident of the normalization ?  normalizes to 1,
c appears to be fine.
      prob = 1./root2pi/sqrt(ysigtotsq) * exp(-ydiff**2/2.0/ysigtotsq)

      if (IFDEBUG .ge. 1) then
c for historical compatibility with calcprob1 write statement
         xnormtot = 1.0
         write(*,1020) f,g,x,y,xnormtot,prob,log(prob)
      end if
 1020 format(4(f6.3,2x),f6.4,2(2x,f9.4))
      return
      end


c You need an erf(x) function, which appears to be built in
c in Linux Fortran.  Or find one from NETLIB or someplace.

