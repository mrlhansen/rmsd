# Martin Hansen (martin.hansen@roma2.infn.it)
# March, 2018
# Copyright (C), all rights reserved.

# Variables
CC = gcc
LDFLAGS = -lm -lmpfr -larb -lflint
CCFLAGS = -Wall -Wextra -O2 -I.

# Do not change the following
SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)
TARGETS = fakecorr realcorr
MODULES = $(filter-out $(addsuffix .o,$(TARGETS)),$(OBJECTS))

all: $(OBJECTS) $(TARGETS)

$(TARGETS): $(OBJECTS)
	$(CC) $(MODULES) $@.o $(LDFLAGS) -o $@

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) $(TARGETS)
