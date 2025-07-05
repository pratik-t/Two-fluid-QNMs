FC = gfortran
SRC = $(wildcard *.f90)
OBJDIR = BUILD
MATHDIR = MATH
FFLAGS = -march=native -flto -O3 -floop-nest-optimize -fopenmp  -fno-strict-aliasing -Wall -I./$(OBJDIR) -J./$(OBJDIR)
OBJ = $(patsubst %.f90, $(OBJDIR)/%.o,$(SRC))
EXEC = Main
MODULES = $(patsubst %.f90,$(OBJDIR)/%.mod,$(SRC))

.SECONDARY: $(OBJ) $(MODULES)

.PHONY: all clean

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(EXEC): $(OBJ)	
	$(FC) $(FFLAGS) $(MATHDIR)/*.o $(OBJ) -o $@ -L$(MATHDIR) -lodepack -llapack -lblas

clean:
	del /Q /F /S $(OBJDIR)\*.o Main
