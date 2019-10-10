#ifndef ANYPROG_RANDOM_HPP
#define ANYPROG_RANDOM_HPP

#include <chrono>
#include <random>

namespace anyprog {

class random {
private:
    std::default_random_engine engine;
    std::uniform_real_distribution<> distribution;

public:
    random();
    random(double l, double u);
    virtual ~random() = default;

public:
    double generate();
};
}

#endif