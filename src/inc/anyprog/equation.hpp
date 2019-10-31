#ifndef ANYPROG_EQUATION
#define ANYPROG_EQUATION

#include "block.hpp"
#include "optimization.hpp"
#include <functional>
#include <memory>
#include <vector>

namespace anyprog {
class equation {
public:
    typedef std::function<double(const real_block&)> funcation_t;

public:
    equation(const std::vector<equation::funcation_t>&, const real_block&);
    virtual ~equation() = default;

public:
    real_block solve(double eps = 1e-8,size_t max_iter = 1000);
    bool is_ok() const;

private:
    bool ok;
    std::shared_ptr<optimization> opt;
};
}

#endif