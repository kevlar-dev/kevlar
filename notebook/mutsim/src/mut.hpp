#ifndef KEVLAR_CPP_MUT
#define KEVLAR_CPP_MUT

#include <string>
#include "hist.hpp"
#include "log.hpp"
#include "hashgraph.hh"
#include "khmer.hh"

using namespace khmer;

class Mutator
{
    protected:
        uint k;
        ulong limit;
        ulong nuclcount;
        Histogram abund_hist;
        Histogram unique_hist;
        Logger& logger;

    public:
        Mutator(uint ksize, Logger& l, uint maxabund = 16, ulong lim = 0);
        virtual unsigned long process(std::string& sequence, Countgraph& countgraph) = 0;
        std::ostream& print(std::ostream& stream) const;
};

std::ostream& operator<<(std::ostream& stream, const Mutator& m);

#endif // KEVLAR_CPP_MUT
