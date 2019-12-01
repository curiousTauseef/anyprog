#include "../help.hpp"



int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [&](const anyprog::real_block& x) {
        return -0.7 * x(2) + 5 * pow(x(0) - 0.5, 2) + 0.8;
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return -exp(x(0) - 0.2) - x(1);
        },
        [](const anyprog::real_block& x) {
            return x(1) + 1.1 * x(2) + 1;
        },
        [](const anyprog::real_block& x) {
            return x(0) - 1.2 * x(2) - 0.2;
        }

    };

    std::vector<anyprog::optimization::range_t> range = { { 0.2, 1 }, { -2.22554, -1 }, { 0, 1 } };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    opt.set_filter_function([](anyprog::real_block& x) {
        x(2) = round(x(2));
    });
    auto ret = opt.search(5, 3);
    anyprog::print(opt.is_ok(), ret, obj);
    return 0;
}