# Compilators

ifdef _USE_CLANG
  CC = clang++
  COMPILER_SPECIFIC = -stdlib=libstdc++ -I/usr/include/c++/8 -I/usr/include/x86_64-linux-gnu/c++/8
  LINKING_SPECIFIC  = -stdlib=libstdc++ -L/usr/lib/gcc/x86_64-linux-gnu/8 -L/usr/lib/clang
else
  CC = g++
  COMPILER_SPECIFIC = 
endif

CFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp

# Number of sites
#
ifdef NUM 
  SITES = -D_SITESNUMBER=${NUM} -D_SILENCE
else 
  SITES = -D_SITESNUMBER=96 -D_SILENCE
endif

# Verbose
# If verbose mode is needed please use:
# make V=-D[verbose flag] instead of make
# Usually I use _VERBOSE as a flag /andrzej

# Linking
LFLAGS = -std=c++17 -lm -lgsl -lgslcblas -fopenmp

 
# File names
EXEC = test.x

# Sources
CPP_SOURCES = $(wildcard *.cpp) $(wildcard ../../../asa/*.cpp) $(wildcard ../../../asa/solver/*.cpp)
CPP_OBJECTS = $(CPP_SOURCES:.cpp=.o)
 
#+++++++++++++++++++#
# COMPILATION RULES #
#+++++++++++++++++++#

# Main target
$(EXEC):  $(CPP_OBJECTS)
	$(CC) $(LINKING_SPECIFIC) $(CPP_OBJECTS) -o $(EXEC) $(LFLAGS)
 
# Other targets 
%.o: %.cpp
	$(CC) $(COMPILER_SPECIFIC) $(CFLAGS) $(SITES) $(VERBOSE) -c $< -o $@

# Cleaning
clean:
	rm -f $(EXEC) $(CPP_OBJECTS) best.flips 

# Force compilation
all:	clean $(EXEC)

# Force compilation and run valgind
test:	all
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 0 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 1 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 2 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 3 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 4 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 5 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 6 1
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC) basis.dat supercell.dat flippable.dat 7 1
	rm -f $(EXEC) $(CPP_OBJECTS) best.flips 
