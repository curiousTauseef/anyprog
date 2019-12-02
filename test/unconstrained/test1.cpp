#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page295.htm
//The global minimum: x* =  (0, …, 0), f(x*) = 0.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double a = 20, b = 0.2, c = 2 * M_PI, s1 = 0, s2 = 0;
        double dim = x.rows();
        for (size_t i = 0; i < dim; ++i) {
            s1 += x(i) * x(i);
            s2 += cos(c * x(i));
        }
        return -a * exp(-b * sqrt(1.0 / dim * s1)) - exp(1.0 / dim * s2) + a + exp(1.0);
    };

    anyprog::optimization::range_t range = { -15, 30 };
    size_t dim = 2;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search(10,3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}