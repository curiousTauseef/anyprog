#ifndef ANYPROG_OPTIMIZATION
#define ANYPROG_OPTIMIZATION

#include "block.hpp"
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace anyprog {
class optimization {
public:
    typedef std::function<double(const real_block&)> funcation_t;
    typedef funcation_t equation_condition_funcation_t;
    typedef funcation_t inequation_condition_funcation_t;
    typedef std::function<void(real_block&)> filter_funcation_t;
    typedef std::function<real_block(const real_block&)> gradient_function_t;
    typedef std::pair<double, double> range_t;
    typedef std::vector<std::pair<double, real_block>> history_t;
    enum method {
        LN_COBYLA = 0,
        LN_NEWUOA,
        LN_NEWUOA_BOUND,
        LN_NELDERMEAD,
        LN_SBPLX,
        LN_AUGLAG,
        LN_AUGLAG_EQ,
        LN_BOBYQA,
        LN_PRAXIS,
        LD_MMA,
        LD_SLSQP,
        LD_LBFGS,
        GN_DIRECT,
        GN_DIRECT_L,
        GN_DIRECT_L_RAND,
        GN_DIRECT_NOSCAL,
        GN_DIRECT_L_NOSCAL,
        GN_DIRECT_L_RAND_NOSCAL,
        GN_ORIG_DIRECT,
        GN_ORIG_DIRECT_L,
        GN_MLSL,
        GN_MLSL_LDS,
        GN_ISRES,
        AUGLAG,
        AUGLAG_EQ,
        G_MLSL,
        G_MLSL_LDS,
        GN_ESCH,
        GN_CRS2_LM,
        GN_AGS
    };
    enum solver_t {
        NLOPT = 0
    };

private:
    class help_t {
    public:
        help_t()
            : fun(0)
            , filter(0)
            , grad(0)
        {
        }
        virtual ~help_t() = default;
        optimization::funcation_t* fun;
        optimization::filter_funcation_t* filter;
        optimization::gradient_function_t* grad;
    };
    solver_t solver;
    double fval;
    bool ok;
    real_block point;
    funcation_t cb;
    filter_funcation_t filter_cb;
    gradient_function_t grad_cb;
    std::vector<equation_condition_funcation_t> eq_fun;
    std::vector<inequation_condition_funcation_t> ineq_fun;
    std::vector<gradient_function_t> eq_grad_fun, ineq_grad_fun;
    std::vector<range_t> range;
    history_t history;
    bool check(const real_block&, double) const;
    int select_nlopt_method(optimization::method) const;

public:
    optimization() = delete;
    optimization(const funcation_t&, const real_block&);
    optimization(const funcation_t&, const std::vector<range_t>& range);
    optimization(const funcation_t&, const real_block&, const std::vector<range_t>& range);
    optimization(const real_block&, const real_block&);
    optimization(const real_block&, const std::vector<range_t>& range);
    optimization(const real_block&, const real_block&, const std::vector<range_t>& range);
    virtual ~optimization() = default;

public:
    optimization& set_equation_condition(const std::vector<equation_condition_funcation_t>&);
    optimization& set_inequation_condition(const std::vector<inequation_condition_funcation_t>&);
    optimization& set_equation_condition(const real_block&, const real_block&);
    optimization& set_inequation_condition(const real_block&, const real_block&);
    optimization& set_filter_function(const filter_funcation_t&);
    optimization& set_gradient_function(const gradient_function_t&);
    optimization& set_equation_gradient_function(const std::vector<gradient_function_t>&);
    optimization& set_inequation_gradient_function(const std::vector<gradient_function_t>&);
    optimization& set_enable_integer_filter();
    optimization& set_enable_binary_filter();
    optimization& set_solver(optimization::solver_t);
    const history_t& get_history() const;
    bool is_ok() const;

public:
    const real_block& solve(optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000);
    const real_block& search(size_t = 100, size_t = 30, double = 0.328, optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000);
    double obj(const real_block&) const;

private:
    const real_block& nlopt_solve(optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000);

private:
    static double instance_fun(unsigned n, const double* x, double* grad, void* my_func_data);
    static double instance_eq_fun(unsigned n, const double* x, double* grad, void* my_func_data);
    static double instance_ineq_fun(unsigned n, const double* x, double* grad, void* my_func_data);

public:
    static optimization::method default_local_method;
    static bool enable_auto_eq_objector, enable_auto_ineq_objector, enable_default_bound_step;
    static double default_bound_step;
    static size_t default_population;
    static size_t max_reloop_iter;

public:
    static real_block fminunc(const optimization::funcation_t&, const real_block&, double = 1e-5, size_t = 1000);
    static real_block fminunc(const optimization::funcation_t&, const optimization::gradient_function_t&, const real_block&, double = 1e-5, size_t = 1000);
    static real_block fminbnd(const optimization::funcation_t&, const std::vector<range_t>&, double = 1e-5, size_t = 1000);
    static real_block assignment(const real_block& c, bool& ok, double& fval, size_t = 20, size_t = 5, double = 0.328, optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000)

    class tsp {
    private:
        const real_block& bk;
        real_block data;
        double max_value, sum;
        std::vector<size_t> path;
        void find(size_t i);

    public:
        tsp() = delete;
        tsp(const real_block&);
        virtual ~tsp() = default;
        const std::vector<size_t>& solve() const;
        double obj() const;
    };
};
}

#endif