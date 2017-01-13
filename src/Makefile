CC = g++
COPTS = -Wall -Wno-deprecated -std=c++11
LKOPTS = -Wall -Wno-deprecated -std=c++11

PRODUCT_NAME = Ximulator

OBJS =\
	hopcroft.o\
	edmond.o\
	sunflowSolver.o\
	util.o\
	global.o\
	scheduler.o\
	schedulerVarys.o\
	schedulerAalo.o\
	schedulerOptc.o\
	schedulerTMS.o\
	schedulerEdmond.o\
	schedulerSunflow.o\
	schedulerSolstice.o\
	schedulerBvN.o\
    trafficGen.o\
	coflow.o\
	events.o\
	main.o

HEADRES =\
	hopcroft.h\
	edmond.h\
	sunflowSolver.h\
	global.h\
	util.h\
	scheduler.h\
	trafficGen.h\
	coflow.h\
	events.h

%.o: %.cc
	$(CC) $(COPTS) -c $< -o $@

all:    simulator

simulator: $(OBJS)
	$(CC) $(LKOPTS) -o $(PRODUCT_NAME) $(OBJS)

$(OBJS): $(HEADRES)

clean:
	rm -f *.o $(PRODUCT_NAME)
	rm -f ./audit.txt
	rm -f ./output.txt