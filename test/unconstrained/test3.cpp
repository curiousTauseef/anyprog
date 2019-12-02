#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page595.htm
//The global minima: x* =  (0, 0), f(x*) = 0

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(x(0), 2) + 2 * pow(x(1), 2) - 0.3 * cos(3 * M_PI * x(0)) - 0.4 * cos(4 * M_PI * x(1)) + 0.7;
    };

    anyprog::optimization::range_t range = { -100, 100 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}