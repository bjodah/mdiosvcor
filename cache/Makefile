all: *.f90 *.h
	for FILE in *.f90; do \
	MODNAME=`echo $${FILE} | cut -d. -f1`; \
	f2py -m $${MODNAME} -c $${MODNAME}.f90; \
	done ;

.PHONY: all