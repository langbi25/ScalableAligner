# .KEEP_STAT:

all: main

CXX		= g++
FLAGS		= -Wall -D NDEBUG -O3 -m64  
# LIB		= -lz -lm -std=c++11 -pthread -lc -lbz2 -llzma -lstdc++  -Wl,--rpath=/home/b8402/22_liangjialang/lib   -Wl,--dynamic-linker=/home/b8402/22_liangjialang/lib/ld-linux-x86-64.so.2
LIB		= -lz -lm -std=c++11 -pthread -lc -lbz2 -llzma -lstdc++  

SOURCE		= main.cpp mapping.cpp bwt_fm_index.cpp KmerAnalysis.cpp nw_alignment.cpp tool.cpp alignment.cpp  mappingpe.cpp
HEADER		= util.h task_queue/concurrentqueue.h task_queue/readerwriterqueue.h task_queue/atomicops.h
OBJECT		= $(SOURCE:%.cpp=%.o)
%.o:		%.cpp $(HEADER)
			$(CXX) $(FLAGS) -c $<

main:		$(OBJECT) 
			$(CXX) $(FLAGS) $(OBJECT) -o preAloc BWT_Index/libbwa.a  htslib/libhts.a $(LIB)

clean:
		rm -f *.o *~