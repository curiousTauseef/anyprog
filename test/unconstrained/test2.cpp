#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page288.htm
//The global minimum: x* =  (3, 0.5), f(x*) = 0.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(1.5 - x(0) * (1 - x(1)), 2) + pow(2.25 - x(0) * (1 - pow(x(1), 2)), 2) + pow(2.625 - x(0) * (1 - pow(x(1), 3)), 2);
    };

    anyprog::optimization::range_t range = { -4.5, 4.5 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}