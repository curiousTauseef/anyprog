#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page913.htm
//The global minima: x* =  (-π , 12.275), (π , 2.275), (9.42478, 2.475),f(x*) = 0.397887

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(x(1) - (5.1 / (4 * M_PI * M_PI)) * pow(x(0), 2) + 5 * x(0) / M_PI - 6, 2) + 10 * (1 - 1 / (8 * M_PI)) * cos(x(0)) + 10;
    };

    std::vector<anyprog::optimization::range_t> range = { { -5, 10 }, { 0, 15 } };

    anyprog::optimization opt(obj, range);
    auto ret = opt.search(10 ,3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}