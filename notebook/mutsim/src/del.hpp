#ifndef KEVLAR_CPP_DEL
#define KEVLAR_CPP_DEL

#include "mut.hpp"

using namespace oxli;

class MutatorDel : public Mutator
{
    private:
        ulong delcount;
        uint delsize;
        class Deletion
        {
            public:
                std::string& sequence;
                MutatorDel& mut;
                std::vector<uint> abunds;
                Deletion(std::string& seq, MutatorDel& m, Counttable& counttable);
                void print(std::ostream& stream);
        };

    public:
        MutatorDel(uint ksize, uint delsize, Logger& l, uint maxabund = 16, ulong lim = 0);
        ulong process(std::string& sequence, Counttable& counttable);
        ulong get_mut_count();
};

#endif // KEVLAR_CPP_DEL
