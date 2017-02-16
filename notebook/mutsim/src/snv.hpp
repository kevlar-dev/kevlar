#ifndef KEVLAR_CPP_SNV
#define KEVLAR_CPP_SNV

#include <array>
#include "mut.hpp"

using namespace khmer;

class MutatorSNV : public Mutator
{
    private:
        ulong nuclcount;
        std::string nucl;
        class SingleNucleotideVariant
        {
            public:
                std::string& sequence;
                MutatorSNV& mut;
                std::vector<uint> abunds;
                SingleNucleotideVariant(std::string& seq, MutatorSNV& m, Countgraph& countgraph);
                void print(std::ostream& stream);
        };

    public:
        MutatorSNV(uint ksize, Logger& l, uint maxabund = 16, ulong lim = 0);
        ulong process(std::string& sequence, Countgraph& countgraph);
        ulong get_mut_count();
};

#endif // KEVLAR_CPP_SNV
