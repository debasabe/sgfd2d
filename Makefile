LIBS =

INCL_DIRS =

## GNU
CC = icpx
#CFLAGS = -O3 -qopenmp ${INCL_DIRS} -std=c++11 -diag-disable=10441
CFLAGS = -O3 -qopenmp

## Intel
# CC = icpc
# CFLAGS = -O3 -fopenmp ${INCL_DIRS} -std=c++11 -diag-disable=10441
# CFLAGS = -g ${INCL_DIRS}
# CFLAGS = -O3 -axAVX -fp-model precise -qopenmp -qopt-report=0 -qopt-report-phase=openmp
# CFLAGS = -Ofast -fp-model precise -qopenmp

O = OBJ/

HEADERS = macros.h cvector.h cmatrix.h

default: sgfdm2d

all: sgfdm2d fdm1d

# Compile all the files
sgfdm2d: $(O)main.o $(O)sgfdm2d.o
	@echo "Linking " $@
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

fdm1d: $(O)fdm1d.o
	@echo "Linking " $@
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# Compile CPP files that have a specific header file
$(O)%.o: %.cpp %.h Makefile $(HEADERS)
	@echo "Creating obj file from" $<
	$(CC) -o $@ -c $< $(CFLAGS)

# Compile CPP files
$(O)%.o: %.cpp Makefile $(HEADERS)
	@echo "Creating obj file from" $<
	$(CC) -o $@ -c $< $(CFLAGS)

# Compile CPP files
$(O)%.o: %.f Makefile $(HEADERS)
	@echo "Creating obj file from" $<
	$(CF) -o $@ -c $< $(CFLAGS)

# Delete all the output files
clrout:
	rm -f -v OUTPUT/*

# Delete the temporay files and executable
clean:
	rm -f sgfdm2d fdm1d $(O)*

remake: clean default

