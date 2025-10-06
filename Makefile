CC=g++
CC+=-DDEBUG -g  
CFLAGS=-c -Wall -m64
LDFLAGS=-fPIC
DIR_SRC = ./src
DIR_BIN = ./bin
SOURCES=read_oscillation_v01.cxx $(wildcard $(DIR_SRC)/*.cxx)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=$(DIR_BIN)/read_oscillation_v01

ROOTSYS=/home/jpmendez/code/root_install

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) 

#CFLAGS += -std=c++11 -I./inc/ -I$(ROOTSYS)/include/
CFLAGS += -I./inc/ -I$(ROOTSYS)/include/

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE):$(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lMinuit2

#.C.o:
%.o:%.cxx
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(DIR_SRC)/*.o; rm $(EXECUTABLE); rm -f *.pcm *.d *.so
canv:
	rm -f canv*
rroot:
	rm -f *.root
