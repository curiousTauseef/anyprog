#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page1240.htm
//The global minima: f(x*) = 0.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double sum = 0;
        for (size_t i = 1; i < x.rows(); ++i) {
            sum += i * pow(2 * x(i) * x(i) - x(i - 1), 2);
        }
        return pow(x(0) - 1, 2) + sum;
    };

    anyprog::optimization::range_t range = { -10, 10 };
    size_t dim = 4;
    anyprog::optimization opt(obj, range, dim);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}