#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 3;
    anyprog::optimization::function_t obj = [&](const anyprog::real_block& x) {
        return pow(x(0) - 1, 2) + pow(x(1) - 2, 2) + pow(x(2) - 4, 2);
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return -x(0) + x(1) + x(2) - 3.5;
        },
        [](const anyprog::real_block& x) {
            return x(0) + x(1) - x(2) - 6;
        }
    };

    anyprog::optimization::range_t range = { 0, 100 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(ineq);
    opt.set_enable_integer_filter();

    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);
    return 0;
}