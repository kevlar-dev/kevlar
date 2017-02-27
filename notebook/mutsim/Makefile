KHMERROOT=khmer
CXXFLAGS=-Wall -O3 --std=c++11 \
		 -I src -I $(KHMERROOT)/lib \
		 -I $(KHMERROOT)/third-party/seqan/core/include \
		 -I $(KHMERROOT)/third-party/smhasher

OBJS=obj/log.o obj/hist.o obj/mut.o obj/snv.o obj/del.o
BINS=bin/mut-hist

all: $(BINS)

$(KHMERROOT)/lib/liboxli.a:
	cd $(KHMERROOT)/lib && make

bin/%: $(OBJS) src/%.cpp $(KHMERROOT)/lib/liboxli.a
	@- echo "[compile $*]"
	@- mkdir -p bin/
	@ g++ $(CXXFLAGS) -o $@ $^

$(OBJS): obj/%.o : src/%.cpp src/%.hpp
	@- echo "[compile $*]"
	@- mkdir -p obj/
	@ g++ $(CXXFLAGS) -c -o $@ $<

clean:
	@- echo "[cleanup]"
	@ rm -f $(OBJS) $(BINS)
