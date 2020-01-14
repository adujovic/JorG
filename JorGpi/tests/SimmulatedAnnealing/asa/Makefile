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

# Verbose
# If verbose mode is needed please use:
# make V=-D[verbose flag] instead of make
# Usually I use _VERBOSE as a flag /andrzej
V =

# Linking
LFLAGS = -std=c++17 -lm -lgsl -lgslcblas 

 
# File names
EXEC = start

# Sources
CPP_SOURCES = $(wildcard *.cpp) $(wildcard ../../../asa/*.cpp)
CPP_OBJECTS = $(CPP_SOURCES:.cpp=.o)
 
#+++++++++++++++++++#
# COMPILATION RULES #
#+++++++++++++++++++#

# Main target
$(EXEC):  $(CPP_OBJECTS)
	$(CC) $(LINKING_SPECIFIC) $(CPP_OBJECTS) -o $(EXEC) $(LFLAGS)
 
# Other targets 
%.o: %.cpp
	$(CC) $(COMPILER_SPECIFIC) $(CFLAGS) -c $< -o $@

# Cleaning
clean:
	rm -f $(EXEC) $(CPP_OBJECTS)  

# Force compilation
all:	clean $(EXEC)

# Force compilation and run valgind
test:	all
	valgrind --leak-check=full --show-leak-kinds=all ./$(EXEC)
