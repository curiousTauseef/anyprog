#ifndef ANYPROG_FIT_HPP
#define ANYPROG_FIT_HPP

#include "block.hpp"
#include "optimization.hpp"
#include <functional>
#include <vector>

namespace anyprog {
class fit {
public:
    typedef std::function<double(const real_block&, const real_block&)> function_t;
    typedef std::function<void(real_block&)> filter_function_t;
    typedef std::function<real_block(const real_block&)> gradient_function_t;

private:
    optimization::solver_t solver;
    real_block dat, X, point;
    function_t cb;
    filter_function_t filter_cb;
    gradient_function_t grad_cb;
    std::vector<optimization::equation_condition_function_t> eq_fun;
    std::vector<optimization::inequation_condition_function_t> ineq_fun;

public:
    fit() = delete;
    fit(const real_block&, size_t);
    fit(const real_block&, const std::vector<function_t>&, const real_block&);
    virtual ~fit() = default;
    fit& set_equation_condition(const std::vector<optimization::equation_condition_function_t>&);
    fit& set_inequation_condition(const std::vector<optimization::inequation_condition_function_t>&);
    fit& set_filter_function(const filter_function_t&);
    fit& set_gradient_function(const gradient_function_t&);
    fit& set_solver(optimization::solver_t);

public:
    const real_block& solve(const real_block&);
    real_block fitting(const real_block&) const;
    const real_block& lssolve(const real_block&, optimization::method = optimization::method::LN_COBYLA, double eps = 1e-5, size_t = 1000);
    const real_block& lssolve(const real_block&, const std::vector<optimization::range_t>& range, optimization::method = optimization::method::LN_COBYLA, double eps = 1e-5, size_t = 1000);
    const real_block& lssearch(const real_block&, const std::vector<optimization::range_t>& range, size_t = 100, size_t = 30, double = 0.382, optimization::method = optimization::method::LN_COBYLA, double eps = 1e-5, size_t = 1000);
    double r_squared(const real_block&,const real_block&) const;
};
}

#endif