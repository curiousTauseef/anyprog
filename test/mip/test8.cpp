#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 5;
    anyprog::optimization::function_t obj = [&](const anyprog::real_block& x) {
        return -(pow(x(0), 2) + pow(x(1), 2) + 3 * pow(x(2), 2) + 4 * pow(x(3), 2) + 2 * pow(x(4), 2) - 8 * x(0) - 2 * x(1) - 3 * x(2) - x(3) - 2 * x(4));
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return x.sum() - 400;
        },
        [](const anyprog::real_block& x) {
            return x(0) + 2 * x(1) + 2 * x(2) + x(3) + 6 * x(4) - 800;
        },
        [](const anyprog::real_block& x) {
            return 2 * x(0) + x(1) + 6 * x(2) - 200;
        },
        [](const anyprog::real_block& x) {
            return x(2) + x(3) + 5 * x(4) - 200;
        }

    };

    anyprog::optimization::range_t range = { 0, 99 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(ineq);
    opt.set_enable_integer_filter();
    auto ret = opt.search(5, 3);
    anyprog::print(opt.is_ok(), ret, obj);
    return 0;
}