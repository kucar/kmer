CXXFLAGS =	-fopenmp -O3 -g -Wall -fmessage-length=0 -std=c++11 -Wno-unused-variable

OBJS =		SBStask.o kmer_hash.o file_op.o
LIBS =

TARGET =	SBStask

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) -fopenmp

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
