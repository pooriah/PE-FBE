include par.make
#m=160 32 64
#bs=40 8 64
#np=16
#tl=600

#build=build
hook=hook/build/libhook.so

$(hook): hook
	make -C hook

#conf: $(hook)

default: cleanall

body.dat par.make par.conf: setup par.py
	./setup

par.py:
	touch $@

#setup.copy: setup
#	cp $< $@

conf: $(hook) par.conf add.conf #setup.copy

#bc.vtk: body.dat Makefile std.conf add.conf
#	make run

#a.jpg: bc.vtk plot.py
#	which pvbatch && pvbatch ./plot.py || touch a.jpg

#.PHONY: default

include $(shell ap.makesim)

clean::
	rm -f a.jpg body.dat par.make par.conf body.dat bc.dat gen.copy bub_generation_hist.txt bubgen.conf *.xmf *.raw *.vtk *.csv electro_hist.txt elec.conf
