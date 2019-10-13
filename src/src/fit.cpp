#include "fit.hpp"
#include <iostream>
namespace anyprog {

fit::fit(const real_block& x, size_t m)
    : solver(optimization::solver_t::NLOPT)
    , dat(x)
    , X(x.rows(), m + 1)
    , point(m + 1, 1)
    , cb()
    , filter_cb()
    , eq_fun()
    , ineq_fun()
{
    for (size_t i = 0; i < X.rows(); ++i) {
        for (size_t j = 0; j < X.cols(); ++j) {
            X(i, j) = pow(x(i, 0), m - j);
        }
    }
    this->cb = [&, m](const real_block& row, const real_block& ret) {
        double sum = 0;
        for (size_t i = 0; i < ret.rows(); ++i) {
            sum += ret(i, 0) * pow(row(0, 0), m - i);
        }
        return sum;
    };
}
fit::fit(const real_block& x, const std::vector<funcation_t>& fun, const real_block& param)
    : solver(optimization::solver_t::NLOPT)
    , dat(x)
    , X(x.rows(), fun.size())
    , point(param)
    , cb()
    , filter_cb()
    , eq_fun()
    , ineq_fun()
{
    for (size_t i = 0; i < X.rows(); ++i) {
        for (size_t j = 0; j < X.cols(); ++j) {
            X(i, j) = fun[j](x.row(i), param);
        }
    }
    this->cb = [&](const real_block& row, const real_block& ret) {
        double sum = 0.0;
        for (size_t i = 0; i < fun.size(); ++i) {
            sum += ret(i, 0) * fun[i](row, ret);
        }
        return sum;
    };
}
fit& fit::set_equation_condition(const std::vector<optimization::equation_condition_funcation_t>& eq_cond)
{
    this->eq_fun = eq_cond;
    return *this;
}
fit& fit::set_inequation_condition(const std::vector<optimization::inequation_condition_funcation_t>& ineq_cond)
{
    this->ineq_fun = ineq_cond;
    return *this;
}
fit& fit::set_filter_function(const fit::filter_funcation_t& cb)
{
    this->filter_cb = cb;
    return *this;
}

const real_block& fit::solve(const real_block& y)
{
    if (X.rows() == X.cols()) {
        this->point = X.ldlt().solve(y);
    }
    this->point = X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
    return this->point;
}

real_block fit::fitting(const real_block& ret) const
{
    real_block Y(this->dat.rows(), 1);
    for (size_t i = 0; i < Y.rows(); ++i) {
        Y(i, 0) = this->cb(this->dat.row(i), ret);
    }
    return Y;
}

fit& fit::set_solver(optimization::solver_t s)
{
    this->solver = s;
    return *this;
}

const real_block& fit::lssolve(const real_block& y, optimization::method m, double eps, size_t max_iter)
{
    optimization::funcation_t obj_fun = [&](const real_block& ret) {
        real_block Y(this->dat.rows(), 1);
        for (size_t i = 0; i < Y.rows(); ++i) {
            Y(i, 0) = this->cb(this->dat.row(i), ret);
        }
        return (Y - y).norm() / (2.0 * Y.rows());
    };
    optimization opt(obj_fun, this->point);
    opt.set_solver(this->solver);
    if (!this->eq_fun.empty()) {
        opt.set_equation_condition(this->eq_fun);
    }
    if (!this->ineq_fun.empty()) {
        opt.set_inequation_condition(this->ineq_fun);
    }

    if (this->filter_cb) {
        opt.set_filter_function(this->filter_cb);
    }

    this->point = opt.solve(m, eps, max_iter);
    return this->point;
}
const real_block& fit::lssolve(const real_block& y, const std::vector<optimization::range_t>& range, optimization::method m, double eps, size_t max_iter)
{
    optimization::funcation_t obj_fun = [&](const real_block& ret) {
        real_block Y(this->dat.rows(), 1);
        for (size_t i = 0; i < Y.rows(); ++i) {
            Y(i, 0) = this->cb(this->dat.row(i), ret);
        }
        return (Y - y).norm() / (2.0 * Y.rows());
    };
    optimization opt(obj_fun, this->point, range);
    opt.set_solver(this->solver);
    if (!this->eq_fun.empty()) {
        opt.set_equation_condition(this->eq_fun);
    }
    if (!this->ineq_fun.empty()) {
        opt.set_inequation_condition(this->ineq_fun);
    }

    if (this->filter_cb) {
        opt.set_filter_function(this->filter_cb);
    }

    this->point = opt.solve(m, eps, max_iter);
    return this->point;
}

const real_block& fit::lssearch(const real_block& y, const std::vector<optimization::range_t>& range, size_t max_random_iter, size_t max_not_changed, double s, optimization::method m, double eps, size_t max_iter)
{
    optimization::funcation_t obj_fun = [&](const real_block& ret) {
        real_block Y(this->dat.rows(), 1);
        for (size_t i = 0; i < Y.rows(); ++i) {
            Y(i, 0) = this->cb(this->dat.row(i), ret);
        }
        return (Y - y).norm() / (2.0 * Y.rows());
    };
    optimization opt(obj_fun, this->point, range);
    opt.set_solver(this->solver);
    if (!this->eq_fun.empty()) {
        opt.set_equation_condition(this->eq_fun);
    }
    if (!this->ineq_fun.empty()) {
        opt.set_inequation_condition(this->ineq_fun);
    }

    if (this->filter_cb) {
        opt.set_filter_function(this->filter_cb);
    }

    this->point = opt.search(max_random_iter, max_not_changed, s, m, eps, max_iter);
    return this->point;
}
}