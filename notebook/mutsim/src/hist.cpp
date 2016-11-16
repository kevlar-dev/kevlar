#include <assert.h>
#include <iostream>
#include "hist.hpp"

Histogram::Histogram(uint max) : maxbin(max)
{
    while (hist.size() < maxbin + 1) {
        hist.push_back(0);
    }
}

void Histogram::increment(uint bin, uint by)
{
    if (bin > maxbin) {
        bin = maxbin;
    }
    hist[bin] += by;
}

ulong Histogram::get(uint bin)
{
    assert(bin <= maxbin);
    return hist[bin];
}

std::ostream& Histogram::print(std::ostream& stream) const
{
    bool first = true;
    stream << '[';
    for (auto value : hist) {
        if (first) {
            first = false;
        }
        else {
            stream << ", ";
        }
        stream << value;
    }
    stream << ']';
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Histogram& h)
{
    return h.print(stream);
}
