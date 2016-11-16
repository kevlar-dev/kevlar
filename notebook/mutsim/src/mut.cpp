#include "mut.hpp"

Mutator::Mutator(uint ksize, Logger& l, uint maxabund, ulong lim)
    : k(ksize), limit(lim), nuclcount(0), abund_hist(maxabund),
      unique_hist(ksize), logger(l)
{

}

std::ostream& Mutator::print(std::ostream& stream) const
{
    stream << abund_hist << '\n' << unique_hist << '\n';
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Mutator& m)
{
    return m.print(stream);
}
