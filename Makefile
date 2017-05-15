CXXFLAGS =	-O3 -g -Wall -fmessage-length=0 -std=c++11 -Wno-unused-variable #-DDEBUG_M

OBJS =		SBStask.o kmer_hash.o
LIBS =

TARGET =	SBStask

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
