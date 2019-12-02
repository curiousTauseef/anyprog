#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2533.htm
//The global minima: x* =  (14.095,0.84296), f(x*) = -6961.81388.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(x(0) - 10, 3) + pow(x(1) - 20, 3);
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return -pow(x(0) - 5, 2) - pow(x(1) - 5, 2) + 100;
        },
        [](const anyprog::real_block& x) {
            return pow(x(0) - 6, 2) + pow(x(1) - 5, 2) - 82.81;
        }
    };

    std::vector<anyprog::optimization::range_t> range = { { 13, 100 }, { 0, 100 } };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}