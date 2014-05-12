EXE = PyOCP_problem
OBJS = PyOCP_problem.o\
      PyOCP_common.o\
       Double_Integrator__0.o\
       Double_Integrator__1.o\
       Double_Integrator__2.o\

ADDLIBS =
ADDINCFLAGS =

CC = cc
CFLAGS = -std=c99 -O3 -pipe -DNDEBUG -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD
CLINKFLAGS =

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=/opt/local/lib64/pkgconfig:/opt/local/lib/pkgconfig:/opt/local/share/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/opt/local/lib64/pkgconfig:/opt/local/lib/pkgconfig:/opt/local/share/pkgconfig: pkg-config --libs ipopt` -lstdc++ -lm

CYGPATH_W = echo


SRC=$(wildcard *.c)

all: $(EXE)

.SUFFIXES: .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CC) $(CFLAGS) $(CLINKFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

subarc: $(SRC)
	$(CC) $(CFLAGS) $(CLINKFLAGS) -o $@ $^ $(LIBS) $(ADDLIBS)

clean:
	rm -f $(EXE) $(OBJS) ipopt.out

.c.o:
	$(CC) $(CFLAGS) $(INCL) -c -o $@ $<

