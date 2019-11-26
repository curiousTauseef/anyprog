#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1760.htm
//The global minima: x* =  (0, -1), f(x*) = 3.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double a = 1 + pow(x(0) + x(1) + 1, 2) * (19 - 14 * x(0) + 3 * x(0) * x(0) - 14 * x(1) + 6 * x(0) * x(1) + 3 * x(1) * x(1));
        double b = 30 + pow(2 * x(0) - 3 * x(1), 2) * (18 - 32 * x(0) + 12 * x(0) * x(0) + 48 * x(1) - 36 * x(0) * x(1) + 27 * x(1) * x(1));
        return a * b;
    };

    anyprog::optimization::range_t range = { -2, 2 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}