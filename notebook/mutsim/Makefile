KHMERROOT=khmer
CXXFLAGS=-Wall -O3 --std=c++11 \
		 -I src -I $(KHMERROOT)/include/oxli \
		 -I src -I $(KHMERROOT)/include/khmer \
		 -I $(KHMERROOT)/third-party/seqan/core/include \
		 -I $(KHMERROOT)/third-party/smhasher \
		 -I $(KHMERROOT)/third-party/cqf \
		 -I $(KHMERROOT)/third-party/rollinghash

OBJS=obj/log.o obj/hist.o obj/mut.o obj/snv.o obj/del.o
BINS=bin/mut-hist

all: $(BINS)

$(KHMERROOT)/src/oxli/liboxli.a:
	cd $(KHMERROOT)/src/oxli && make

bin/%: $(OBJS) src/%.cpp $(KHMERROOT)/src/oxli/liboxli.a
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
