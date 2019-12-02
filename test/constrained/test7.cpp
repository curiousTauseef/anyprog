#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page898.htm
//The global minima: x∗ = (2.330499, 1.951372,−0.4775414, 4.365726,−0.6244870, 1.038131,1.594227), f (x∗) = 680.6300573.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(x(0) - 10, 2) + 5 * pow(x(1) - 12, 2) + pow(x(2), 4) + 3 * pow(x(3) - 11, 2) + 10 * pow(x(4), 6) + 7 * pow(x(5), 2) + pow(x(6), 4) - 4 * x(5) * x(6) - 10 * x(5) - 8 * x(6);
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            double v1 = 2 * x(0) * x(0), v2 = x(1) * x(1);
            return v1 + 3 * v2 * v2 + x(2) + 4 * x(3) * x(3) + 5 * x(4) - 127;
        },
        [](const anyprog::real_block& x) {
            return 7 * x(0) + 3 * x(1) + 10 * x(2) * x(2) + x(3) - x(4) - 282;
        },
        [](const anyprog::real_block& x) {
            double v1 = 2 * x(0) * x(0), v2 = x(1) * x(1);
            return 23 * x(0) + v2 + 6 * x(5) * x(5) - 8 * x(6) - 196;
        },
        [](const anyprog::real_block& x) {
            double v1 = 2 * x(0) * x(0), v2 = x(1) * x(1);
            return 2 * v1 + v2 - 3 * x(0) * x(1) + 2 * x(2) * x(2) + 5 * x(5) - 11 * x(6);
        }
    };

    size_t dim = 7;
    anyprog::optimization::range_t range = { -10, 10 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}