# Compilators

ifdef _USE_CLANG
  CPP = clang++
  COMPILER_SPECIFIC = -stdlib=libstdc++ -I/usr/include/c++/8 -I/usr/include/x86_64-linux-gnu/c++/8
  LINKING_SPECIFIC  = -stdlib=libstdc++ -L/usr/lib/gcc/x86_64-linux-gnu/8 -L/usr/lib/clang
else
  CPP = g++
  COMPILER_SPECIFIC = 
endif

CFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp

# Number of sites
#
# SITES

# Verbose
# If verbose mode is needed please use:
# make V=-D[verbose flag] instead of make
# Usually I use _VERBOSE as a flag /andrzej

# Linking
LFLAGS = -std=c++17 -lm -lgsl -lgslcblas 

 
# File names
EXEC = start

# Sources
CPP_SOURCES = $(wildcard *.cpp) $(wildcard ../*.cpp)
CPP_OBJECTS = $(CPP_SOURCES:.cpp=.o)
 
#+++++++++++++++++++#
# COMPILATION RULES #
#+++++++++++++++++++#

# Main target
$(EXEC):  $(CPP_OBJECTS)
	$(CPP) $(LINKING_SPECIFIC) $(CPP_OBJECTS) -o $(EXEC) $(LFLAGS)
 
# Other targets 
%.o: %.cpp
	$(CPP) $(COMPILER_SPECIFIC) $(CFLAGS) $(SITES) $(VERBOSE) -c $< -o $@

# Cleaning
clean:
	rm -f $(EXEC) $(CPP_OBJECTS)  

# Spacing
spacing:
	@echo "#"
	@echo "#"
	@echo "#"
	@echo "#"
	@echo "#"

# Force compilation
all:	spacing clean $(EXEC)

# Force compilation and run valgind
check:	all
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC)
