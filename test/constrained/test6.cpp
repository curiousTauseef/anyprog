#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page3389.htm
//The global minima: x∗ = (1.2279713, 4.2453733) , f(x∗) = 0.095825.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return -(pow(sin(2 * M_PI * x(0)), 3) * sin(2 * M_PI * x(1))) / (pow(x(0), 3) * (x(0) + x(1)));
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return x(0) * x(0) - x(1) + 1;
        },
        [](const anyprog::real_block& x) {
            return 1 - x(0) + pow(x(1) - 4, 2);
        }
    };

    size_t dim = 2;
    anyprog::optimization::range_t range = { 0, 10 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 5);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}