
INCLUDE = eigen/ -Ihelpers.h
BASE = shapeedit

all: $(BASE)

OS := $(shell uname -s)

ifeq ($(OS), Linux) # Science Center Linux Boxes
  CPPFLAGS = -I/home/l/i/lib175/usr/glew/include
  LDFLAGS += -L/home/l/i/lib175/usr/glew/lib -L/usr/X11R6/lib
  LIBS += -lGL -lGLU -lglut
endif

ifeq ($(OS), Darwin) # Assume OS X
  CPPFLAGS += -D__MAC__ -stdlib=libstdc++
  LDFLAGS += -framework GLUT -framework OpenGL
endif

ifdef OPT 
  #turn on optimization
  CXXFLAGS += -O3
else 
  #turn on debugging
  CXXFLAGS += -g
endif

CXX = g++ 

OBJ = $(BASE).o ppm.o glsupport.o poisson.o

$(BASE): $(OBJ)
	$(LINK.cpp) -I$(INCLUDE) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJ) $(BASE)
