# .KEEP_STAT:

all: preAloc

htslib:
		mkdir -p bin/
		$(MAKE) -C src/htslib libhts.a

bwt_index:
		$(MAKE) -C src/BWT_Index && mv -f src/BWT_Index/$@ bin/

# kart: htslib bwt_index
# 		$(MAKE) -C src && mv -f src/$@ bin/

preAloc:  bwt_index
		$(MAKE) -C src && mv -f src/$@ bin/		 

clean:
		rm -f preAloc bwt_index
		$(MAKE) clean -C src
		# $(MAKE) clean -C src/htslib
		$(MAKE) clean -C src/BWT_Index
