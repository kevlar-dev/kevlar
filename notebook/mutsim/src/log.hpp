#ifndef KEVLAR_CPP_LOGGER
#define KEVLAR_CPP_LOGGER

#include <iostream>

class Logger
{
    public:
        unsigned int interval;
        unsigned long threshold;
        unsigned long counter;
        std::ostream& stream;

        Logger();
        Logger(unsigned int intvl, std::ostream& os);
        bool increment(unsigned int by = 1);
};

#endif // KEVLAR_CPP_LOGGER
