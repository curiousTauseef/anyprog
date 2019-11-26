#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page506.htm
//The global minima: x* = (1,1,…,1,3,3,3,1), f(x*) = -15.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double a = 0, b = 0, c = 0;
        for (size_t i = 0; i < 4; ++i) {
            a += x(i);
        }
        for (size_t i = 0; i < 4; ++i) {
            b += x(i) * x(i);
        }
        for (size_t i = 4; i < x.rows(); ++i) {
            c += x(i);
        }
        return 5 * a - 5 * b - c;
    };
    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return 2 * x(0) + 2 * x(1) + x(9) + x(10) - 10;
        },
        [](const anyprog::real_block& x) {
            return 2 * x(0) + 2 * x(2) + x(9) + x(11) - 10;
        },
        [](const anyprog::real_block& x) {
            return 2 * x(1) + 2 * x(2) + x(10) + x(11) - 10;
        },
        [](const anyprog::real_block& x) {
            return -8 * x(0) + x(9);
        },
        [](const anyprog::real_block& x) {
            return -8 * x(1) + x(10);
        },
        [](const anyprog::real_block& x) {
            return -8 * x(2) + x(11);
        },
        [](const anyprog::real_block& x) {
            return -2 * x(3) - x(4) + x(9);
        },
        [](const anyprog::real_block& x) {
            return -2 * x(5) - x(6) + x(10);
        },
        [](const anyprog::real_block& x) {
            return -2 * x(7) - x(8) + x(11);
        }
    };

    size_t dim = 13;
    anyprog::real_block u(dim, 1);
    u.fill(1);
    u(9) = 100;
    u(10) = 100;
    u(11) = 100;

    std::vector<anyprog::optimization::range_t> range;
    for (size_t i = 0; i < dim; ++i) {
        range.push_back({ 0, u(i) });
    }

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}