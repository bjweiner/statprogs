
c read data from a text file and fit a gaussian to it

      program fitgausfile

      parameter (NMAX=10000)

      real xd(NMAX),yd(NMAX),yerr(NMAX),yfit(NMAX)
      real pars(4),errpars(4),dydp(4)
      integer ifit(4)
      character fname*80

c file is assumed sorted with x increasing

      call pgbeg(0,"?",1,1)
      call pgscf(2)
      call pgsch(1.5)
      
 100  continue
      write(*,'("File with x,y,yerr [quit]: ",$)')
      read(*,'(a)') fname
      if (fname(1:3) .eq. '   ') go to 666
      open(2,file=fname,status='old',err=100)

      np = 1
 200  continue
      read(2,*,err=220,end=220) xd(np),yd(np),yerr(np)
      np = np + 1
      go to 200
 220  continue
      np = np-1
      close(2)
      write(*,'("Read ",i5," points")') np

c plot data
c      wxmin = 1.e10
c      wxmax = -1.e10
      wxmin = xd(1)
      wxmax = xd(np)
      wymin = 1.e10
      wymax = -1.e10
      do i=1,np
         wymin = min(wymin,yd(i))
         wymax = max(wymax,yd(i))
      end do
      wxmarg = (wxmax-wxmin) / 100.
      wymarg = (wymax-wymin) / 40.
      call pgenv(wxmin-wxmarg,wxmax+wxmarg,wymin-wymarg,
     $     wymax+wymarg,0,1)
      call pgpt(np,xd,yd,17)
      call pgerrb(6,np,xd,yd,yerr,0.0)

      write(*,'("x range to fit: ",$)')
      read(*,*) xmin,xmax
      call locate(xd,np,xmin,ixmin)
      ixmin = max(ixmin,1)
      call locate(xd,np,xmax,ixmax)
      ixmax = min(ixmax,np)
      nfit = ixmax-ixmin+1

      iguess=1
      do i=1,4
         ifit(i) = 1
      end do
      call fitgaus(xd(ixmin),yd(ixmin),yerr(ixmin),nfit,
     $     pars,errpars,ifit,iguess,chisq)
      do i=1,4
         write (*,*) pars(i)," +- ", errpars(i)
      end do
      write(*,*) "chisq = ",chisq

c compute and plot fit
      do i=ixmin,ixmax
         call gausfunc(xd(i),pars,yfit(i),dydp,4)
      end do
      call pgqci(indexc)
      call pgsci(2)
      call pgline(nfit,xd(ixmin),yfit(ixmin))
      call pgsci(indexc)

c write out data, fit, subtracted

      open(3,file='fitgausfile.out',status='unknown')
      do i=1,np
c         call gausfunc(xd(i),pars,gfitval,dydp,4)
c         write(3,*) i,xd(i),yd(i),yerr(i),gfitval,yd(i)-gfitval
         call gausfunc(xd(i),pars,yfit(i),dydp,4)
         write(3,*) i,xd(i),yd(i),yerr(i),yfit(i),yd(i)-yfit(i)
      end do
      close(3)

      go to 100

 666  continue

      end

      
