#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1113.htm
// The global minima: x* =  (0, …, 0), f(x*) = 0.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double s = 0;
        for (size_t i = 0; i < x.rows(); ++i) {
            s += x(i) * x(i);
        }
        return s;
    };

    anyprog::optimization::range_t range = { -5.12, 5.12 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}