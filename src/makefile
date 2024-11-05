# .KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -g   -m64  $(shell pkg-config --cflags jellyfish-2.0 sdsl-lite)
# LIB		= -lz -lm -std=c++11 -pthread -lc -lbz2 -llzma -lstdc++  -Wl,--rpath=/home/b8402/22_liangjialang/lib   -Wl,--dynamic-linker=/home/b8402/22_liangjialang/lib/ld-linux-x86-64.so.2
LDFLAGS=$(shell pkg-config --libs sdsl-lite) $(shell pkg-config --libs jellyfish-2.0) -lpthread -lz $(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L/-Wl,-rpath,/')

LIB		= -lz -lm -std=c++11 -pthread -lc -lbz2 -llzma -lstdc++  

SOURCE		= main.cpp mapping.cpp bwt_fm_index.cpp KmerAnalysis.cpp nw_alignment.cpp tool.cpp alignment.cpp  mappingpe.cpp  SBT/Kmers.cpp SBT/BF.cpp SBT/BloomTree.cpp SBT/util.cpp SBT/Query.cpp
HEADER		= util.h task_queue/concurrentqueue.h task_queue/readerwriterqueue.h task_queue/atomicops.h SBT/Kmers.h SBT/BF.h SBT/BloomTree.h SBT/util.h SBT/Query.h
OBJECT		= $(SOURCE:%.cpp=%.o)
%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $< $(LDFLAGS)

main:		$(OBJECT) 
			$(CXX) $(FLAGS) $(OBJECT) -o preAloc BWT_Index/libbwa.a  htslib/libhts.a $(LIB) $(LDFLAGS)

clean:
		rm -f *.o *~