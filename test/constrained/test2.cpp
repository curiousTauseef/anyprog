#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2613.htm
//The global minima: f(x*) = 1.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return -pow(sqrt(x.rows()), x.rows()) * x.prod();
    };
    std::vector<anyprog::optimization::inequation_condition_function_t> eq = {
        [](const anyprog::real_block& x) {
            double sum = 0;
            for (size_t i = 0; i < x.rows(); ++i) {
                sum += x(i) * x(i);
            }
            return sum - 1;
        }
    };

    size_t dim = 20;
    anyprog::optimization::range_t range = { 0, 1 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_equation_condition(eq);
    auto ret = opt.search(10, 5);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}