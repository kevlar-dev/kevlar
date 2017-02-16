#ifndef KEVLAR_CPP_MUT
#define KEVLAR_CPP_MUT

#include <random>
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

        double sampling_rate;
        std::uniform_real_distribution<double> dist;
        std::random_device rand;
        std::mt19937 prng;

    public:
        Mutator(uint ksize, Logger& l, uint maxabund = 16, ulong lim = 0);
        virtual unsigned long process(std::string& sequence, Countgraph& countgraph) = 0;
        std::ostream& print(std::ostream& stream) const;
        bool skip_nucl();
        void set_sampling_rate(double rate, int seed);
        virtual ulong get_mut_count();
};

std::ostream& operator<<(std::ostream& stream, const Mutator& m);

#endif // KEVLAR_CPP_MUT
