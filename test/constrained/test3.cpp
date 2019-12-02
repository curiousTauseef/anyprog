#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1582.htm
//The global minima: x* =  (78,33,29.995,45,36.7758),f(x*) = -30665.539.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return 5.3578547 * x(2) * x(2) + 0.8356891 * x(0) * x(4) + 37.293239 * x(0) - 40792.141;
    };
    auto u = [](const anyprog::real_block& x) {
        return 85.334407 + 0.0056858 * x(1) * x(4) + 0.0006262 * x(0) * x(3) - 0.0022053 * x(2) * x(4);
    };
    auto v = [](const anyprog::real_block& x) {
        return 80.51249 + 0.0071317 * x(1) * x(4) + 0.0029955 * x(0) * x(1) + 0.0021813 * x(2) * x(2);
    };
    auto w = [](const anyprog::real_block& x) {
        return 9.300961 + 0.0047026 * x(2) * x(4) + 0.0012547 * x(0) * x(2) + 0.0019085 * x(2) * x(3);
    };
    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [&](const anyprog::real_block& x) {
            return -u(x);
        },
        [&](const anyprog::real_block& x) {
            return u(x) - 92;
        },
        [&](const anyprog::real_block& x) {
            return -v(x) + 90;
        },
        [&](const anyprog::real_block& x) {
            return v(x) - 110;
        },
        [&](const anyprog::real_block& x) {
            return -w(x) + 20;
        },
        [&](const anyprog::real_block& x) {
            return w(x) - 25;
        }
    };

    size_t dim = 5;
    anyprog::real_block L(dim, 1), U(dim, 1);
    L << 78, 33, 27, 27, 27;
    U << 102, 45, 45, 45, 45;

    std::vector<anyprog::optimization::range_t> range;
    for (size_t i = 0; i < dim; ++i) {
        range.push_back({ L(i), U(i) });
    }

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}