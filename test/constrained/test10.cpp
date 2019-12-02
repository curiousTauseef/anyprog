#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page4447.htm
//The global minima: x∗ = (−1.717143, 1.595709, 1.827247, −0.7636413, −0.763645) , f (x∗) = 0.0539498.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return exp(x.prod());
    };

    std::vector<anyprog::optimization::equation_condition_function_t> eq = {
        [](const anyprog::real_block& x) {
            double sum = 0;
            for (size_t i = 0; i < x.rows(); ++i) {
                sum += x(i) * x(i);
            }
            return sum - 10;
        },
        [](const anyprog::real_block& x) {
            return x(1) * x(2) - 5 * x(3) * x(4);
        },
        [](const anyprog::real_block& x) {
            return pow(x(0), 3) + pow(x(1), 3) + 1;
        },
    };

    std::vector<anyprog::optimization::range_t> range = {
        { -2.3, 2.3 }, { -2.3, 2.3 }, { -3.2, 3.2 }, { -3.2, 3.2 }, { -3.2, 3.2 }
    };

    anyprog::optimization opt(obj, range);
    opt.set_equation_condition(eq);
    auto ret = opt.search(10, 5);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}