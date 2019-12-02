#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1882.htm
//The global minima: f(x*) = -186.7309.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double s1 = 0, s2 = 0;
        for (size_t i = 1; i <= 5; ++i) {
            s1 += i * cos((i + 1) * x(0) + i);
            s2 += i * cos((i + 1) * x(1) + i);
        }
        return s1 * s2;
    };

    anyprog::optimization::range_t range = { -10, 10 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search(10, 8);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}