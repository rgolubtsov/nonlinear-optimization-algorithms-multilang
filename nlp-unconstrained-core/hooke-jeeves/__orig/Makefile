# nlp-unconstrained-core/hooke-jeeves/__orig/Makefile
# =====================================================
#  Project: nonlinear-optimization-algorithms-multilang
#           nlp-unconstrained-core
#           hooke-jeeves
#           __orig
# Platform: Win32
# =====================================================
# Copyright (C) 2008 Radislav (Radic) Golubtsov

PROJ0 = hooke
PROJ1 = $(PROJ0)-woods

OBJS0 = $(PROJ0).o
OBJS1 = $(PROJ1).o

CC = gcc
RM = del

$(PROJ0).exe: $(OBJS0)
	$(CC) -o $@ $<

$(OBJS1): CFLAGS = -DWoods
$(OBJS1): %-woods.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(PROJ1).exe: $(OBJS1)
	$(CC) -o $@ $<

.PHONY: rosenbrock woods all clean

rosenbrock: $(PROJ0).exe

woods: $(PROJ1).exe

all: rosenbrock woods

clean:
	$(RM) $(PROJ0).exe $(OBJS0) $(PROJ1).exe $(OBJS1)
