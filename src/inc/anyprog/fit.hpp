#ifndef ANYPROG_FIT_HPP
#define ANYPROG_FIT_HPP

#include "block.hpp"
#include "optimization.hpp"
#include <functional>
#include <vector>

namespace anyprog {
class fit {
public:
    typedef std::function<double(const real_block&, const real_block&)> funcation_t;
    typedef std::function<void(real_block&)> filter_funcation_t;

private:
    optimization::solver_t solver;
    real_block dat, X, point;
    funcation_t cb;
    filter_funcation_t filter_cb;
    std::vector<optimization::equation_condition_funcation_t> eq_fun;
    std::vector<optimization::inequation_condition_funcation_t> ineq_fun;

public:
    fit() = delete;
    fit(const real_block&,size_t);
    fit(const real_block&, const std::vector<funcation_t>&, const real_block&);
    virtual ~fit() = default;
    fit& set_equation_condition(const std::vector<optimization::equation_condition_funcation_t>&);
    fit& set_inequation_condition(const std::vector<optimization::inequation_condition_funcation_t>&);
    fit& set_filter_function(const filter_funcation_t&);
    fit& set_solver(optimization::solver_t);

public:
    const real_block& solve(const real_block&);
    real_block fitting(const real_block&) const;
    const real_block& lssolve(const real_block&, optimization::method = optimization::method::LN_COBYLA, double eps = 1e-5, size_t = 1000);
    const real_block& lssolve(const real_block&, const std::vector<optimization::range_t>& range, optimization::method = optimization::method::LN_COBYLA, double eps = 1e-5, size_t = 1000);
    const real_block& lssearch(const real_block&, const std::vector<optimization::range_t>& range, size_t = 100, size_t = 30, double = 0.328, optimization::method = optimization::method::LN_COBYLA, double eps = 1e-5, size_t = 1000);
};
}

#endif