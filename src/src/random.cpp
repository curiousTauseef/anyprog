#include "random.hpp"

namespace anyprog {
random::random()
    : engine(std::chrono::system_clock::now().time_since_epoch().count())
    , distribution(0, 1)
{
}
random::random(double l, double u)
    : engine(std::chrono::system_clock::now().time_since_epoch().count())
    , distribution(l, u)
{
}

double random::generate()
{
    return this->distribution(this->engine);
}
}