if DEBUG
AM_CFLAGS = -g -O0
AM_CXXFLAGS = -g -O0
else
AM_CFLAGS = -O2
AM_CXXFLAGS = -O2
endif

if PROFILING
AM_CFLAGS += -pg
AM_CXXFLAGS += -pg
endif

bin_PROGRAMS=gtenricher

gtenricher_SOURCES=main.cc
