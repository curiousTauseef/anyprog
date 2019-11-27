#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page4001.htm
// The global minima: x* =  ±(1/20.5, 1/2), f(x*) = 0.75.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return x(0) * x(0) + (x(1) - 1) * (x(1) - 1);
    };

    std::vector<anyprog::optimization::equation_condition_function_t> eq = {
        [](const anyprog::real_block& x) {
            return x(1) - x(0) * x(0);
        }
    };

    size_t dim = 2;

    anyprog::optimization::range_t range = { -1, 1 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_equation_condition(eq);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}