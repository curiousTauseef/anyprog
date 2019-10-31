#include "equation.hpp"

namespace anyprog {
equation::equation(const std::vector<equation::funcation_t>& eq, const real_block& p)
    : ok(false)
    , opt()
{
    this->opt = std::make_shared<optimization>([](const real_block& x) {
        return x.sum();
    },
        p);
    this->opt->set_equation_condition(eq);
}
bool equation::is_ok() const
{
    return this->ok;
}
real_block equation::solve(double eps, size_t max_iter)
{
    real_block ret =  this->opt->solve(optimization::method::LN_COBYLA, eps, max_iter);
    this->ok = this->opt->is_ok();
    return ret;
}
}