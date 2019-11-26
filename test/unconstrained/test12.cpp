#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2537.htm
//The global minima: x* =  (1, …, 1), f(x*) = 0.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double s = 0;
        for (size_t i = 0; i < x.rows() - 1; ++i) {
            s += 100 * pow(x(i) * x(i) - x(i + 1), 2) + pow(x(i) - 1, 2);
        }
        return s;
    };

    anyprog::optimization::range_t range = { -5, 10 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search();
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}