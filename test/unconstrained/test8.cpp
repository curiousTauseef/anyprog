#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1361.htm
//The global minima: x* =  (π, π), f(x*) = - 1.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return -cos(x(0)) * cos(x(1)) * exp(-pow(x(0) - M_PI, 2) - pow(x(1) - M_PI, 2));
    };

    anyprog::optimization::range_t range = { -100, 100 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search(100, 80);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}