#ifndef KEVLAR_CPP_HIST
#define KEVLAR_CPP_HIST

#include <vector>
typedef unsigned long ulong;
typedef unsigned int uint;

class Histogram
{
    protected:
        std::vector<ulong> hist;
        uint maxbin;

    public:
        Histogram(uint max = 16);
        void increment(uint bin, uint by = 1);
        ulong get(uint bin);
        std::ostream& print(std::ostream& stream) const;
};

std::ostream& operator<<(std::ostream& stream, const Histogram& h);

#endif // KEVLAR_CPP_HIST
