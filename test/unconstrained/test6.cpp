#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1016.htm
//The global minimum: x* =  (1, …, 1), f(x*) = 0

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return 100 * pow(pow(x(0), 2) - x(1), 2)
            + pow(x(0) - 1, 2) + pow(x(2) - 1, 2) + 90 * pow(pow(x(2), 2) - x(3), 2)
            + 10.1 * (pow(x(1) - 1, 2) + pow(x(3) - 1, 2)) + 19.8 * (x(1) - 1) * (x(3) - 1);
    };

    anyprog::optimization::range_t range = { -10, 10 };
    size_t dim = 4;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}