CXXFLAGS =	-O3 -g -Wall -fmessage-length=0 -std=c++11 -Wno-unused-variable

OBJS =		SBStask.o kmer_hash.o

TESTDIR = tests

LIBS =

TARGET =	SBStask

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)
		(cd tests; make all)

clean:
	rm -f $(OBJS) $(TARGET)
