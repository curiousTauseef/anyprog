#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page3235.htm
//The global minima: x∗ = (579.3167, 1359.943, 5110.071, 182.0174, 295.5985, 217.9799, 286.4162,395.5979), f (x∗) = 7049.3307.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return x(0) + x(1) + x(2);
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return -1 + 0.0025 * (x(3) + x(5));
        },
        [](const anyprog::real_block& x) {
            return -1 + 0.0025 * (-x(3) + x(4) + x(6));
        },
        [](const anyprog::real_block& x) {
            return -1 + 0.01 * (-x(4) + x(7));
        },
        [](const anyprog::real_block& x) {
            return 100 * x(0) - x(0) * x(5) + 833.33252 * x(3) - 83333.333;
        },
        [](const anyprog::real_block& x) {
            return x(1) * x(3) - x(1) * x(6) - 1250 * x(3) + 1250 * x(4);
        },
        [](const anyprog::real_block& x) {
            return x(2) * x(4) - x(2) * x(7) - 2500 * x(4) + 1250000;
        }
    };

    size_t dim = 8;
    anyprog::real_block L(dim, 1), U(dim, 1);
    L << 10, 100, 100, 1, 1, 1, 1, 1;
    U << 10, 10, 10, 1, 1, 1, 1, 1;

    std::vector<anyprog::optimization::range_t> range;
    for (size_t i = 0; i < dim; ++i) {
        range.push_back({ 10 * L(i), 1000 * U(i) });
    }

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}