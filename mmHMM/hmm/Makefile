# 
#  Copyright (C) 2008 Cold Spring Harbor Laboratory
#  Authors: Andrew D. Smith
# 
#  This file is part of rmap.
# 
#  rmap is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  rmap is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with rmap; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
### ==========================

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

COMMON_DIR = ../common
INCLUDEDIRS =  $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))
LIBS = -lgsl -lgslcblas

PROGS = hmr-cthmm hmr-vdhmm

CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
CXXFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
INCLUDEARGS += -I/usr/local/include
LIBS = -L/usr/local/lib
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
	               then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif



all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, \
	smithlab_os.o smithlab_utils.o GenomicRegion.o OptionParser.o) \

hmr-cthmm: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o \
	TwoStateCTHMM.o Distro.o BetaBin.o false_discovery_rate.o numerical-utils.o) \
	$(addprefix $(SMITHLAB_CPP)/, Smoothing.o)

hmr-vdhmm: $(addprefix $(COMMON_DIR)/, MethpipeFiles.o \
  TwoStateVDHMM.o Distro.o BetaBin.o false_discovery_rate.o numerical-utils.o) \
  $(addprefix $(SMITHLAB_CPP)/, Smoothing.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean

