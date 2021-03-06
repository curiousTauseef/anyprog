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
    typedef std::function<double(const real_block&)> function_t;
    typedef function_t equation_condition_function_t;
    typedef function_t inequation_condition_function_t;
    typedef std::function<void(real_block&)> filter_function_t;
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
        optimization::function_t* fun;
        optimization::filter_function_t* filter;
        optimization::gradient_function_t* grad;
    };
    solver_t solver;
    double fval;
    bool ok;
    real_block point;
    function_t cb;
    filter_function_t filter_cb;
    gradient_function_t grad_cb;
    std::vector<equation_condition_function_t> eq_fun;
    std::vector<inequation_condition_function_t> ineq_fun;
    std::vector<gradient_function_t> eq_grad_fun, ineq_grad_fun;
    std::vector<range_t> range;
    history_t history;
    bool check(const real_block&, double) const;
    int select_nlopt_method(optimization::method) const;

public:
    optimization() = delete;
    optimization(const function_t&, const real_block&);
    optimization(const function_t&, const std::vector<range_t>& range);
    optimization(const function_t&, const range_t& range, size_t dim);
    optimization(const function_t&, const real_block&, const std::vector<range_t>& range);
    optimization(const real_block&, const real_block&);
    optimization(const real_block&, const std::vector<range_t>& range);
    optimization(const real_block&, const range_t& range);
    optimization(const real_block&, const real_block&, const std::vector<range_t>& range);
    virtual ~optimization() = default;

public:
    optimization& set_equation_condition(const std::vector<equation_condition_function_t>&);
    optimization& set_inequation_condition(const std::vector<inequation_condition_function_t>&);
    optimization& set_equation_condition(const real_block&, const real_block&);
    optimization& set_inequation_condition(const real_block&, const real_block&);
    optimization& set_filter_function(const filter_function_t&);
    optimization& set_gradient_function(const gradient_function_t&);
    optimization& set_equation_gradient_function(const std::vector<gradient_function_t>&);
    optimization& set_inequation_gradient_function(const std::vector<gradient_function_t>&);
    optimization& set_enable_integer_filter();
    optimization& set_enable_binary_filter();
    optimization& set_enable_integer_filter(const std::vector<size_t>&);
    optimization& set_enable_binary_filter(const std::vector<size_t>&);
    optimization& set_solver(optimization::solver_t);
    const history_t& get_history() const;
    bool is_ok() const;

public:
    const real_block& solve(optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000);
    const real_block& search(size_t = 100, size_t = 30, double = 0.382, optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000);
    double obj(const real_block&) const;

private:
    const real_block& nlopt_solve(optimization::method = optimization::method::LN_COBYLA, double = 1e-5, size_t = 1000);

private:
    static double instance_fun(unsigned n, const double* x, double* grad, void* my_func_data);
    static double instance_eq_fun(unsigned n, const double* x, double* grad, void* my_func_data);
    static double instance_ineq_fun(unsigned n, const double* x, double* grad, void* my_func_data);

public:
    static optimization::method default_local_method;
    static bool enable_default_bound_step;
    static double default_bound_step;
    static size_t default_population;
    static size_t max_reloop_iter;

public:
    static real_block fminunc(const optimization::function_t&, const real_block&, bool&, double = 1e-5, size_t = 1000);
    static real_block fminunc(const optimization::function_t&, const optimization::gradient_function_t&, const real_block&, bool&, double = 1e-5, size_t = 1000);
    static real_block fminbnd(const optimization::function_t&, const std::vector<range_t>&, bool&, double = 1e-5, size_t = 1000);
    static real_block fminbnd(const optimization::function_t&, const range_t&, size_t, bool&, double = 1e-5, size_t = 1000);

    static void print(bool ok, const real_block& ret, const optimization::function_t& obj);
    static void print(bool ok, const real_block& ret, const real_block& obj);

    class assignment {
    private:
        const real_block& bk;
        real_block data;
        double max_value, sum;
        std::vector<std::pair<size_t, size_t>> path;
        void find(size_t i);

    public:
        assignment() = delete;
        assignment(const real_block& c, double inf = 1e10);
        virtual ~assignment() = default;
        const std::vector<std::pair<size_t, size_t>>& solve() const;
        double obj() const;
    };

    class tsp {
    private:
        size_t start;
        const real_block& bk;
        real_block data;
        double max_value, sum;
        std::vector<size_t> path;
        void find(size_t i);

    public:
        typedef std::pair<double, double> point2d_t;
        typedef std::function<double(const point2d_t&, const point2d_t)> point2d_distance_function_t;
        tsp() = delete;
        tsp(const real_block&, size_t = 0, double inf = 1e10);
        virtual ~tsp() = default;
        const std::vector<size_t>& solve() const;
        double obj() const;

    public:
        static real_block distance(const std::vector<point2d_t>&, const point2d_distance_function_t&, double inf = 1e10);
    };
};
}
#endif