#    Makefile for subHMR
#
#    Copyright (C) 2015-2016 Andrew D. Smith
#
#    Authors: Andrew D. Smith
#
#    This code is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

SUBHMR_LIBDIR = ./lib

PROGS = SubunitFinder UnitSignal CollapseIntervals BinarizeCpG AssignCpGs \
				DetectBreaks
SOURCES = $(wildcard *.cpp)

INCLUDEDIRS = $(SMITHLAB_CPP) $(SUBHMR_LIBDIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -L/usr/local/lib -lgsl -lgslcblas # -lefence


CXX = g++
CXXFLAGS = -Wall -Wextra -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

# Flags passed to the C++ compiler.
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
CXXFLAGS += 
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

$(PROGS): $(addprefix $(SUBHMR_LIBDIR)/, ProcSubunit.o)


UnitSignal: $(addprefix $(SUBHMR_LIBDIR)/, UnitCluster.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
