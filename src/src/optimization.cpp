#include "optimization.hpp"
#include "random.hpp"
#include <iostream>
#include <nlopt.h>

namespace anyprog {

optimization::method optimization::default_local_method = optimization::method::LN_COBYLA;
bool optimization::enable_auto_eq_objector = false;
bool optimization::enable_auto_ineq_objector = false;
bool optimization::enable_default_bound_step = true;
double optimization::default_bound_step = 10;
size_t optimization::default_population = 200;
size_t optimization::max_reloop_iter = 3;

double optimization::instance_fun(unsigned n, const double* x, double* grad, void* my_func_data)
{
    real_block ret(n, 1);
    for (size_t i = 0; i < n; ++i) {
        ret(i, 0) = x[i];
    }
    optimization::help_t* help = (optimization::help_t*)(my_func_data);
    if (*help->filter) {
        (*help->filter)(ret);
    }
    if (grad && *help->grad) {
        real_block tmp = (*help->grad)(ret);
        for (size_t i = 0; i < n; ++i) {
            grad[i] = tmp(i, 0);
        }
    }

    return (*help->fun)(ret);
}

double optimization::instance_eq_fun(unsigned n, const double* x, double* grad, void* my_func_data)
{
    real_block ret(n, 1);
    for (size_t i = 0; i < n; ++i) {
        ret(i, 0) = x[i];
    }
    optimization::help_t* help = (optimization::help_t*)(my_func_data);
    if (*help->filter) {
        (*help->filter)(ret);
    }
    return (*help->fun)(ret);
}

double optimization::instance_ineq_fun(unsigned n, const double* x, double* grad, void* my_func_data)
{
    real_block ret(n, 1);
    for (size_t i = 0; i < n; ++i) {
        ret(i, 0) = x[i];
    }
    optimization::help_t* help = (optimization::help_t*)(my_func_data);
    if (*help->filter) {
        (*help->filter)(ret);
    }
    return (*help->fun)(ret);
}

real_block optimization::fminunc(const optimization::funcation_t& obj, const real_block& p, double eps, size_t max_iter)
{
    optimization opt(obj, p);
    return opt.solve(optimization::method::LN_NEWUOA, eps, max_iter);
}
real_block optimization::fminunc(const optimization::funcation_t& obj, const optimization::gradient_function_t& grad, const real_block& p, double eps, size_t max_iter)
{
    optimization opt(obj, p);
    opt.set_gradient_function(grad);
    return opt.solve(optimization::method::LD_LBFGS, eps, max_iter);
}
real_block optimization::fminbnd(const optimization::funcation_t& obj, const std::vector<range_t>& range, double eps, size_t max_iter)
{
    optimization opt(obj, range);
    return opt.solve(optimization::method::LN_NEWUOA_BOUND, eps, max_iter);
}

optimization::optimization(const funcation_t& fun, const real_block& p)
    : solver(optimization::solver_t::NLOPT)
    , fval(0)
    , ok(false)
    , point(p)
    , cb(fun)
    , filter_cb()
    , grad_cb()
    , eq_fun()
    , ineq_fun()
    , range()
    , history()
{
    if (optimization::enable_default_bound_step) {
        double bound_step = fabs(optimization::default_bound_step);
        for (size_t i = 0; i < this->point.rows(); ++i) {
            this->range.push_back({ p(i, 0) - bound_step, p(i, 0) + bound_step });
        }
    }
}

optimization::optimization(const funcation_t& fun, const real_block& p, const std::vector<optimization::range_t>& range)
    : solver(optimization::solver_t::NLOPT)
    , fval(0)
    , ok(false)
    , point(p)
    , cb(fun)
    , filter_cb()
    , grad_cb()
    , eq_fun()
    , ineq_fun()
    , range(range)
    , history()
{
}

optimization::optimization(const funcation_t& fun, const std::vector<optimization::range_t>& range)
    : solver(optimization::solver_t::NLOPT)
    , fval(0)
    , ok(false)
    , point(range.size(), 1)
    , cb(fun)
    , filter_cb()
    , grad_cb()
    , eq_fun()
    , ineq_fun()
    , range(range)
    , history()
{
    random rng(0, 1);
    for (size_t i = 0; i < this->point.rows(); ++i) {
        this->point(i, 0) = this->range[i].first + (this->range[i].second - this->range[i].first) * rng.generate();
    }
}

optimization::optimization(const real_block& v, const real_block& p)
    : solver(optimization::solver_t::NLOPT)
    , fval(0)
    , ok(false)
    , point(p)
    , cb()
    , filter_cb()
    , grad_cb()
    , eq_fun()
    , ineq_fun()
    , range()
    , history()
{
    if (optimization::enable_default_bound_step) {
        double bound_step = fabs(optimization::default_bound_step);
        for (size_t i = 0; i < this->point.rows(); ++i) {
            this->range.push_back({ p(i, 0) - bound_step, p(i, 0) + bound_step });
        }
    }
    this->cb = [&](const real_block& x) {
        size_t m = x.rows();
        double sum = 0.0;
        for (size_t i = 0; i < m; ++i) {
            sum += v(i, 0) * x(i, 0);
        }
        return sum;
    };
}
optimization::optimization(const real_block& v, const std::vector<range_t>& range)
    : solver(optimization::solver_t::NLOPT)
    , fval(0)
    , ok(false)
    , point(range.size(), 1)
    , cb()
    , filter_cb()
    , grad_cb()
    , eq_fun()
    , ineq_fun()
    , range(range)
    , history()
{
    this->cb = [&](const real_block& x) {
        size_t m = x.rows();
        double sum = 0.0;
        for (size_t i = 0; i < m; ++i) {
            sum += v(i, 0) * x(i, 0);
        }
        return sum;
    };
    random rng(0, 1);
    for (size_t i = 0; i < this->point.rows(); ++i) {
        this->point(i, 0) = this->range[i].first + (this->range[i].second - this->range[i].first) * rng.generate();
    }
}
optimization::optimization(const real_block& v, const real_block& p, const std::vector<range_t>& range)
    : solver(optimization::solver_t::NLOPT)
    , fval(0)
    , ok(false)
    , point(p)
    , cb()
    , filter_cb()
    , grad_cb()
    , eq_fun()
    , ineq_fun()
    , range(range)
    , history()
{
    this->cb = [&](const real_block& x) {
        size_t m = x.rows();
        double sum = 0.0;
        for (size_t i = 0; i < m; ++i) {
            sum += x(i, 0) * v(i, 0);
        }
        return sum;
    };
}

optimization& optimization::set_equation_condition(const std::vector<equation_condition_funcation_t>& eq_cond)
{
    this->eq_fun = eq_cond;
    return *this;
}
optimization& optimization::set_inequation_condition(const std::vector<inequation_condition_funcation_t>& ineq_cond)
{
    this->ineq_fun = ineq_cond;
    return *this;
}

optimization& optimization::set_equation_condition(const real_block& A, const real_block& b)
{
    size_t m = A.rows();
    for (size_t i = 0; i < m; ++i) {
        this->eq_fun.emplace_back([&, i](const real_block& x) {
            double v = 0;
            size_t n = x.rows();
            for (size_t j = 0; j < n; ++j) {
                v += A(i, j) * x(j, 0);
            }
            return v - b(i, 0);
        });
    }
    return *this;
}
optimization& optimization::set_inequation_condition(const real_block& A, const real_block& b)
{
    size_t m = A.rows();
    for (size_t i = 0; i < m; ++i) {
        this->ineq_fun.emplace_back([&, i](const real_block& x) {
            double v = 0;
            size_t n = x.rows();
            for (size_t j = 0; j < n; ++j) {
                v += A(i, j) * x(j, 0);
            }
            return v - b(i, 0);
        });
    }
    return *this;
}

optimization& optimization::set_solver(optimization::solver_t s)
{
    this->solver = s;
    return *this;
}

bool optimization::is_ok() const
{
    return this->ok;
}

double optimization::obj(const real_block& ret) const
{
    return this->cb(ret);
}
bool optimization::check(const real_block& p, double eps) const
{
    bool eq_check = true, ineq_check = true;
    size_t m = this->eq_fun.size();
    for (size_t i = 0; eq_check && i < m; ++i) {
        eq_check = eq_check && fabs(this->eq_fun[i](p)) <= eps;
    }
    if (eq_check) {
        double v = 0;
        m = this->ineq_fun.size();
        for (size_t i = 0; ineq_check && i < m; ++i) {
            v = this->ineq_fun[i](p);
            if (v > 0) {
                ineq_check = ineq_check && v <= eps;
            }
        }
    } else {
        return eq_check;
    }
    return eq_check && ineq_check;
}

optimization& optimization::set_filter_function(const optimization::filter_funcation_t& cb)
{
    this->filter_cb = cb;
    return *this;
}

optimization& optimization::set_gradient_function(const optimization::gradient_function_t& cb)
{
    this->grad_cb = cb;
    return *this;
}

const real_block& optimization::search(size_t max_random_iter, size_t max_not_changed, double s, optimization::method m, double eps, size_t max_iter)
{
    if (!this->range.empty()) {
        size_t dim = this->range.size();
        double obj_value, global_obj_value;
        if (this->filter_cb) {
            this->filter_cb(this->point);
        }
        obj_value = this->cb(this->point);
        global_obj_value = obj_value;
        real_block global_point = this->point;
        size_t not_changed = 0, global_max_random_iter = 0, reloop_iter = 0;
        std::vector<std::shared_ptr<random>> rng;
        random rg(0.0 - eps, fabs(s) + eps);
        std::vector<range_t> range_bk = this->range;
    reloop:
        for (size_t i = 0; i < dim; ++i) {
            rng.emplace_back(std::make_shared<random>(this->range[i].first - eps, this->range[i].second + eps));
        }
    loop:
        for (size_t i = 0; i < max_random_iter; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                this->point(j, 0) = rng[j]->generate();
            }
            this->point = this->solve(m, eps, max_iter);
            obj_value = this->fval;
            bool gcheck = this->check(global_point, eps), lcheck = this->check(this->point, eps);
            bool case1 = !gcheck && lcheck, case2 = lcheck && (global_obj_value - obj_value) >= eps;
            if (this->ok && (case1 || case2)) {
                global_point = this->point;
                global_obj_value = obj_value;
                not_changed = 0;
                this->history.push_back({ global_obj_value, global_point });
            } else if (++not_changed > max_not_changed) {
                not_changed = 0;
                break;
            }
        }

        rng.clear();
        not_changed = 0;
        this->point = global_point;
        this->ok = !this->history.empty();
        for (size_t i = 0; i < dim; ++i) {
            range_t& p = range_bk[i];
            double c = p.second - p.first;
            if (c >= eps) {
                double best = this->point(i, 0);
                if (best > c / 2.0) {
                    p.first += (best - p.first) * rg.generate();
                } else {
                    p.second -= (p.second - best) * rg.generate();
                }
            } else {
                ++not_changed;
            }
            rng.emplace_back(std::make_shared<random>(p.first - eps, p.second + eps));
        }
        if ((not_changed < dim - 1) && ++global_max_random_iter <= max_not_changed) {
            not_changed = 0;
            goto loop;
        }
        if (!this->ok && ++reloop_iter <= optimization::max_reloop_iter) {
            global_max_random_iter = 0;
            not_changed = 0;
            range_bk = this->range;
            goto reloop;
        }
        return this->point;
    }
    return this->solve(m, eps, max_iter);
}
int optimization::select_nlopt_method(optimization::method m) const
{
    nlopt_algorithm method;
    switch (m) {
    case optimization::method::LN_COBYLA:
        method = NLOPT_LN_COBYLA;
        break;
    case optimization::method::LN_NEWUOA:
        method = NLOPT_LN_NEWUOA;
        break;
    case optimization::method::LN_NEWUOA_BOUND:
        method = NLOPT_LN_NEWUOA_BOUND;
        break;
    case optimization::method::LN_NELDERMEAD:
        method = NLOPT_LN_NELDERMEAD;
        break;
    case optimization::method::LN_SBPLX:
        method = NLOPT_LN_SBPLX;
        break;
    case optimization::method::LN_AUGLAG:
        method = NLOPT_LN_AUGLAG;
        break;
    case optimization::method::LN_AUGLAG_EQ:
        method = NLOPT_LN_AUGLAG_EQ;
        break;
    case optimization::method::LN_BOBYQA:
        method = NLOPT_LN_BOBYQA;
        break;
    case optimization::method::LN_PRAXIS:
        method = NLOPT_LN_PRAXIS;
        break;
    case optimization::method::LD_MMA:
        method = NLOPT_LD_MMA;
        break;
    case optimization::method::LD_SLSQP:
        method = NLOPT_LD_SLSQP;
        break;
    case optimization::method::LD_LBFGS:
        method = NLOPT_LD_LBFGS;
        break;
    case optimization::method::GN_DIRECT:
        method = NLOPT_GN_DIRECT;
        break;
    case optimization::method::GN_DIRECT_L:
        method = NLOPT_GN_DIRECT_L;
        break;
    case optimization::method::GN_DIRECT_L_RAND:
        method = NLOPT_GN_DIRECT_L_RAND;
        break;
    case optimization::method::GN_DIRECT_NOSCAL:
        method = NLOPT_GN_DIRECT_NOSCAL;
        break;
    case optimization::method::GN_DIRECT_L_NOSCAL:
        method = NLOPT_GN_DIRECT_L_NOSCAL;
        break;
    case optimization::method::GN_DIRECT_L_RAND_NOSCAL:
        method = NLOPT_GN_DIRECT_L_RAND_NOSCAL;
        break;
    case optimization::method::GN_ORIG_DIRECT:
        method = NLOPT_GN_ORIG_DIRECT;
        break;
    case optimization::method::GN_ORIG_DIRECT_L:
        method = NLOPT_GN_ORIG_DIRECT_L;
        break;
    case optimization::method::GN_MLSL:
        method = NLOPT_GN_MLSL;
        break;
    case optimization::method::GN_MLSL_LDS:
        method = NLOPT_GN_MLSL_LDS;
        break;
    case optimization::method::GN_ISRES:
        method = NLOPT_GN_ISRES;
        break;
    case optimization::method::AUGLAG:
        method = NLOPT_AUGLAG;
        break;
    case optimization::method::AUGLAG_EQ:
        method = NLOPT_AUGLAG_EQ;
        break;
    case optimization::method::G_MLSL:
        method = NLOPT_G_MLSL;
        break;
    case optimization::method::G_MLSL_LDS:
        method = NLOPT_G_MLSL_LDS;
        break;
    case optimization::method::GN_ESCH:
        method = NLOPT_GN_ESCH;
        break;
    case optimization::method::GN_CRS2_LM:
        method = NLOPT_GN_CRS2_LM;
        break;
    case optimization::method::GN_AGS:
        method = NLOPT_GN_AGS;
        break;
    default:
        method = NLOPT_LN_COBYLA;
        break;
    }
    return method;
}

const real_block& optimization::nlopt_solve(optimization::method m, double eps, size_t max_iter)
{
    size_t dim = this->point.rows();
    nlopt_algorithm loc_method = (nlopt_algorithm)this->select_nlopt_method(optimization::default_local_method), method = (nlopt_algorithm)this->select_nlopt_method(m);
    nlopt_opt opt_loc = nlopt_create(loc_method, dim);
    nlopt_opt opt = nlopt_create(method, dim);
    nlopt_set_local_optimizer(opt, opt_loc);
    help_t obj;
    obj.filter = &this->filter_cb;
    obj.grad = &this->grad_cb;
    funcation_t auto_obj = [&](const real_block& x) {
        double sum = 0.0;
        if (optimization::enable_auto_eq_objector) {
            for (auto& i : this->eq_fun) {
                sum += pow(i(x), 2);
            }
        }
        if (optimization::enable_auto_ineq_objector) {
            for (auto& i : this->ineq_fun) {
                sum += pow(std::max(0.0, i(x)), 2);
            }
        }

        return this->cb(x) + sum / eps;
    };
    obj.fun = &auto_obj;
    nlopt_set_min_objective(opt, optimization::instance_fun, &obj);
    nlopt_set_xtol_rel(opt, eps);
    nlopt_set_ftol_abs(opt, eps);
    nlopt_set_force_stop(opt, eps);
    nlopt_set_maxeval(opt, max_iter);
    nlopt_set_population(opt, optimization::default_population);
    double lb[dim], ub[dim];
    if (!this->range.empty()) {
        for (size_t i = 0; i < dim; ++i) {
            range_t& cur_range = this->range[i];
            lb[i] = cur_range.first;
            ub[i] = cur_range.second;
        }
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);
    }

    std::vector<help_t> eq_help, ineq_help;
    for (size_t i = 0; i < this->eq_fun.size(); ++i) {
        help_t h;
        h.filter = &this->filter_cb;
        h.fun = &this->eq_fun[i];
        eq_help.emplace_back(h);
    }
    for (size_t i = 0; i < eq_help.size(); ++i) {
        nlopt_add_equality_constraint(opt, instance_eq_fun, &eq_help[i], eps);
    }

    for (size_t i = 0; i < this->ineq_fun.size(); ++i) {
        help_t h;
        h.filter = &this->filter_cb;
        h.fun = &this->ineq_fun[i];
        ineq_help.emplace_back(h);
    }
    for (size_t i = 0; i < ineq_help.size(); ++i) {
        nlopt_add_inequality_constraint(opt, instance_ineq_fun, &ineq_help[i], eps);
    }

    double ret[dim];
    for (size_t i = 0; i < dim; ++i) {
        ret[i] = this->point(i, 0);
    }

    if (nlopt_optimize(opt, ret, &this->fval) >= 0) {
        this->ok = true;
        for (size_t i = 0; i < dim; ++i) {
            this->point(i, 0) = ret[i];
        }
    }
    if (this->filter_cb) {
        this->filter_cb(this->point);
    }
    this->ok = this->ok && this->check(this->point, eps);
    nlopt_destroy(opt);
    nlopt_destroy(opt_loc);
    return this->point;
}

const real_block& optimization::solve(optimization::method m, double eps, size_t max_iter)
{
    if (this->solver == optimization::solver_t::NLOPT) {
        return this->nlopt_solve(m, eps, max_iter);
    }
    return this->nlopt_solve(m, eps, max_iter);
}

const optimization::history_t& optimization::get_history() const
{
    return this->history;
}

optimization& optimization::set_enable_integer_filter()
{
    this->filter_cb = [&](real_block& x) {
        size_t m = x.rows();
        for (size_t i = 0; i < m; ++i) {
            x(i, 0) = round(x(i, 0));
        }
    };
    return *this;
}
optimization& optimization::set_enable_binary_filter()
{
    if (this->range.empty()) {
        for (size_t i = 0; i < this->point.rows(); ++i) {
            this->range.push_back({ 0, 1 });
        }
    } else {
        for (auto& i : this->range) {
            i.first = 0;
            i.second = 1;
        }
    }

    this->filter_cb = [&](real_block& x) {
        size_t m = x.rows();
        for (size_t i = 0; i < m; ++i) {
            x(i, 0) = round(x(i, 0));
        }
    };
    return *this;
}
}