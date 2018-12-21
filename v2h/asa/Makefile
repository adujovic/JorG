# Compilators
#CC = clang
CC = g++
CFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp

# Verbose
# If verbose mode is needed please use:
# make V=-D[verbose flag] instead of make
# Usually I use _VERBOSE as a flag /andrzej
V =

# Linking
LFLAGS = -lm -lgsl -lgslcblas 
 
# File names
EXEC = start

# Sources
C_SOURCES = $(wildcard *.cpp)
C_OBJECTS = $(C_SOURCES:.cpp=.o)
 
#+++++++++++++++++++#
# COMPILATION RULES #
#+++++++++++++++++++#

# Main target
$(EXEC):  $(C_OBJECTS)
	$(CC) $(C_OBJECTS) -o $(EXEC) $(LFLAGS)
 
# Other targets 
%.o: %.cpp
	$(CC) -c $(CFLAGS) -c $< -o $@

# Cleaning
clean:
	rm -f $(EXEC) $(C_OBJECTS)  

