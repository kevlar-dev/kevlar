#include "log.hpp"

Logger::Logger()
    : interval(1000000), threshold(1000000), counter(0), stream(std::cerr)
{

}

Logger::Logger(unsigned int intvl, std::ostream& os)
    : interval(intvl), threshold(intvl), counter(0), stream(os)
{

}

bool Logger::increment(unsigned int by)
{
    counter += by;
    if (threshold > 0 && counter > threshold) {
        threshold += interval;
        return true;
    }
    return false;
}
