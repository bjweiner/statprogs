# statprogs
Statistics and fitting programs, mostly least squares / max likelihood, mostly Fortran.

The lstsq and mlsfit programs can be used for fitting lines to data with errors in both coordinates and intrinsic scatter.  See http://mingus.mmto.arizona.edu/~bjw/tkrs_kinematics/fitprogs.html
and the discussion in the appendix of B.J. Weiner et al, 2006, ApJ, 653, 1049, http://adsabs.harvard.edu/abs/2006ApJ...653.1049W  . If you find the discussion useful or use these programs to fit your data, please cite that paper.

The lstsq program uses PGPLOT to plot the output. It will still run without pgplot, but plotting is nice. PGPLOT is free. If you have difficulty compiling it, eg on OS X, instructions are in another of my Github repositories, https://github.com/bjweiner/PGPLOT-Mac-OSX-config


