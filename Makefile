CC=g++
CC+=-DDEBUG -g
CFLAGS=-c -Wall -m64
LDFLAGS=-fPIC

DIR_SRC = ./src
DIR_BIN = ./bin

SOURCES=read_oscillation_v01.cxx $(wildcard $(DIR_SRC)/*.cxx)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=$(DIR_BIN)/read_oscillation_v01

# -------- ROOT --------
# Remote
ROOTSYS=/home/jpmendez/code/root_install
# Local
# ROOTSYS=/usr

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs)

CFLAGS += -I./inc/ -I$(ROOTSYS)/include/

# -------- Remote Config --------
REMOTE = pg8
REMOTE_DIR = ~/analysis/test-remote-analysis

# -------- Remote Build --------

# Sync only changed source + headers
push:
	rsync -ravz \
		--include="*/" \
		--include="*.cxx" \
		--include="*.h" \
		--include="*.hpp" \
        --include="*.icc" \
		--exclude="*" \
		./ $(REMOTE):$(REMOTE_DIR) \
    && scp Makefile.remote $(REMOTE):$(REMOTE_DIR)/Makefile

# Always sync Makefile explicitly
# rsync -avz Makefile $(REMOTE):$(REMOTE_DIR)

# Compile on remote machine
remote-build:
	ssh $(REMOTE) "cd $(REMOTE_DIR) && make"

# Combined
deploy: push remote-build

# -------- Cleanup --------

clean:
	rm -f *.o $(DIR_SRC)/*.o
	rm -f $(EXECUTABLE)
	rm -f *.pcm *.d *.so

canv:
	rm -f canv*

rroot:
	rm -f *.root
