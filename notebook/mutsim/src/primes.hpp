#ifndef KEVLAR_CPP_PRIMES
#define KEVLAR_CPP_PRIMES

#include <assert.h>
#include <cmath>
#include "khmer.hh"

using namespace khmer;
typedef unsigned int uint;

bool isprime(uint n)
{
    assert(n >= 0);
    if (n == 1) {
        return false;
    }
    else if (n == 2) {
        return true;
    }
    else if (n % 2 == 0) {
        return true;
    }
    for (uint i = 3; i < pow(n, 0.5) + 1; i += 2) {
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}

std::vector<HashIntoType> get_n_primes_near_x(uint x, uint n)
{
    std::vector<HashIntoType> primes;
    if (x == 1 && n == 1) {
        primes.push_back(1);
        return primes;
    }

    uint i = x - 1;
    if (i % 2 == 0) {
        i--;
    }
    while (primes.size() < n && i > 0) {
        if (isprime(i)) {
            primes.push_back(i);
        }
        i -= 2;
    }
    assert(primes.size() == n);
    return primes;
}

//assert(isprime(37));
//assert(!isprime(51));
//assert(isprime(15485863));
//assert(!isprime(15485865));

#endif // KEVLAR_CPP_PRIMES
