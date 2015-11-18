# VERSION NUMBER
vers=3.5.6
# Build type 32 or 64 (for linux versions)
TOBUILD=64
# contact email address to go into exe error messages
email=arsf-processing@pml.ac.uk

CPPFLAGS=-Wall -O4  -D VERSION='"$(vers)"' -D CONTACTEMAIL='"$(email)"'
# rather than `pkg-config --cflags blitz` set up path to version 0.9 of blitz as APL is incompatible with later versions
CPPFLAGS += -Iexternal_code/blitz-0.9/
LDFLAGS=

# don't actually need to link to blitz because we're only using the template functions defined in the .h files
# LDFLAGS=$(LDFLAGS) `pkg-config --libs blitz`

# just for apltran (note you need to yum install proj-devel.i686 for 32bit builds)
transform_ldflags=-lproj

ifeq ($(TOBUILD),32)
    # 32 bit build needs large file support
	CPPFLAGS += -m32 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -ffloat-store
endif

# Compiler to use
CC=g++
AR=ar -cq

# Object file directory
obj=objectfiles

# source dir
src=src

# lib dir
libs=libs

# exe dir
bin=bin

# file to save dependencies to
dependfile = .depend

common_libs = $(libs)/libcommandline.a $(libs)/liblogger.a $(libs)/libbinaryreader.a

# Make all the exes
all: $(dependfile) aplcal aplnav aplcorr apltran aplmap aplmask aplshift  

# make individual exes
aplcal:  $(bin)/aplcal
aplnav:  $(bin)/aplnav
aplcorr: $(bin)/aplcorr
apltran: $(bin)/apltran
aplmap: $(bin)/aplmap
aplmask: $(bin)/aplmask
aplshift: $(bin)/aplshift

# make command line library
$(libs)/libcommandline.a: $(obj)/commandline.o $(obj)/commonfunctions.o
	rm -f $@
	$(AR) $@ $^

# make logger library
$(libs)/liblogger.a: $(obj)/logger.o $(obj)/commonfunctions.o
	rm -f $@
	$(AR) $@ $^

# make binary reader library
$(libs)/libbinaryreader.a: $(obj)/binaryreader.o $(obj)/bil.o $(obj)/bsq.o $(obj)/binfile.o $(obj)/multifile.o $(obj)/commonfunctions.o
	rm -f $@
	$(AR) $@ $^

$(bin)/aplcal: $(obj)/radcal.o $(obj)/specimsensors.o $(obj)/calibration.o $(obj)/mainworker.o $(obj)/commonfunctions.o $(obj)/bilwriter.o $(obj)/os_dependant.o $(common_libs)
	$(CC)  $(CPPFLAGS) -o $@ $^

$(bin)/aplnav: $(obj)/navigation.o $(obj)/navfileclasses.o $(obj)/datahandler.o $(obj)/navigationsyncer.o $(obj)/navigationinterpolator.o $(obj)/interpolationfunctions.o $(obj)/leverbore.o $(obj)/transformations.o $(obj)/conversions.o $(obj)/commonfunctions.o $(obj)/bilwriter.o  $(obj)/os_dependant.o $(common_libs)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

$(bin)/aplcorr: $(obj)/geolocation.o $(obj)/geodesics.o $(obj)/cartesianvector.o $(obj)/dems.o $(obj)/viewvectors.o $(obj)/navbaseclass.o $(obj)/conversions.o $(obj)/planarsurface.o $(obj)/transformations.o $(obj)/leverbore.o $(obj)/commonfunctions.o $(obj)/bilwriter.o  $(common_libs)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

$(bin)/apltran: $(obj)/bilwriter.o $(obj)/commonfunctions.o $(obj)/transform.o $(obj)/basic_igm_worker.o $(common_libs)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(transform_ldflags) -o $@ $^

$(bin)/aplmap: $(obj)/bilwriter.o $(obj)/commonfunctions.o $(obj)/level3grid.o $(obj)/basic_igm_worker.o $(obj)/map_main.o $(obj)/TreeGrid.o $(obj)/treegrid_support.o $(obj)/os_dependant.o $(obj)/geodesics.o $(obj)/conversions.o $(common_libs)  
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

$(bin)/aplmask: $(obj)/bilwriter.o $(obj)/applymask.o $(common_libs)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

$(bin)/aplshift: $(obj)/bilwriter.o $(obj)/navshift.o $(obj)/datahandler.o $(obj)/interpolationfunctions.o $(obj)/navbaseclass.o $(common_libs)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

clean: 
	rm -f $(obj)/*.o 
	rm -f $(dependfile) 

cleanall: clean
	rm -f $(libs)/*.a
	rm -f $(bin)/* 

depend: $(dependfile)

# Get the dependencies for each cpp file and store in dependfile
$(dependfile): $(src)/*.cpp
	rm -f $(dependfile)

	for fname in $(wildcard $(src)/*.cpp) ; do \
		fnameo=`echo $$fname | sed 's/.cpp/.o/' | sed 's/$(src)\///'`;\
		$(CC) -MM -MT $(obj)/$$fnameo $(CPPFLAGS) $$fname -MF $(dependfile) ; \
   done

# These are targets that dont actually exist after the build
.PHONY: clean cleanall all depend

.SUFFIXES:.o .cpp

# make the object file using the equivalent named cpp file and dependencies
# $< means first dependency, $@ means file on the left of :
$(obj)/%.o: $(src)/%.cpp 
	$(CC) $(CPPFLAGS) -c $< -o $@

-include $(dependfile) 
