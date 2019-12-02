#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1621.htm
//The global minima: x* =  (0.0898, -0.7126), (-0.0898, 0.7126), f(x*) = 0

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return 1.0316285 + 4 * x(0) * x(0) - 2.1 * pow(x(0), 4) + pow(x(0), 6) / 3 + x(0) * x(1) - 4 * x(1) * x(1) + 4 * pow(x(1), 4);
    };

    anyprog::optimization::range_t range = { -5, 5 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}