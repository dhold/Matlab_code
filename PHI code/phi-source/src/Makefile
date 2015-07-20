include ../make.config
PROJECT= phi
PLATFORM= mac
SRCS= HierarchyIntegrator.cpp PhiParameters.cpp ComplexMatrix.cpp HierarchyNode.cpp Barrier.cpp main.cpp
OBJS= ${SRCS:.cpp=.o}

.cpp.o: 
	$(CC) $(CFLAGS) $(OPTIMISE) -c $<


all: $(PROJECT)


$(PROJECT): $(OBJS)
	$(LD) $(LDFLAGS) $(OPTIMISE) -o $@ $(OBJS)	

debug:
	$(CC) $(CFLAGS) -g $(LDFLAGS) $(SRCS) -o $(PROJECT)

static:
	$(CC) $(CFLAGS) $(LDFLAGS_LIN_STATIC) $(SRCS) $(CBLAS_STATIC) -o $(PROJECT)
	

clean:
	rm -f *.o
	rm -f phi

