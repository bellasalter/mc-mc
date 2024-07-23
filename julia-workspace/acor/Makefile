##--------------------------------
## targets

TARGETS = acor

ifeq ($(HOSTTYPE),intel-mac)
  ALL = $(patsubst %,%.dylib,$(TARGETS))
else
  # assume linux
  ALL = $(patsubst %,%.so,$(TARGETS))
endif

##--------------------------------
## commands

ifeq ($(HOSTTYPE),intel-mac)
  CC = cc
else
  # assume linux
  CC = gcc
endif


CC += -O3\
      -Wall\

RM = /bin/rm

##--------------------------------
## rules

all: $(ALL)

clean:
	$(RM) $(ALL)

again: clean all

%.dylib: %.c %.h
	$(CC) -dynamiclib -o $@ $<

%.so: %.c %.h
	$(CC) -shared -o $@ $<
