
include $(HOME)/programs/Make.defs

mlsfit: mlsfit.f
	g77 -O -o mlsfit mlsfit.f $(PGPLOT) 
	