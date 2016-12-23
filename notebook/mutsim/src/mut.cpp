#include "mut.hpp"

Mutator::Mutator(uint ksize, Logger& l, uint maxabund, ulong lim)
    : k(ksize), limit(lim), nuclcount(0), abund_hist(maxabund),
      unique_hist(ksize), logger(l), sampling_rate(1.0), dist(0.0, 1.0), prng(0)
{

}

std::ostream& Mutator::print(std::ostream& stream) const
{
    stream << abund_hist << '\n' << unique_hist << '\n';
    return stream;
}

bool Mutator::skip_nucl()
{
    if (fabs(1.0 - sampling_rate) < 0.0001) {
        return false;
    }

    double random = dist(prng);
    if (random < sampling_rate) {
        return false;
    }
    return true;
}

void Mutator::set_sampling_rate(double rate, int seed)
{
    prng.seed(seed);
    dist.reset();
    sampling_rate = rate;
}

std::ostream& operator<<(std::ostream& stream, const Mutator& m)
{
    return m.print(stream);
}
