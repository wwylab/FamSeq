CC=g++
CFLAGS=-O3 -g3 -c -Wall
LDFLAGS=
SOURCES=FamSeq.cpp family.cpp normal.cpp checkInput.cpp file.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=FamSeq

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $< -o $@