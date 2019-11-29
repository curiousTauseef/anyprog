#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 2;
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return -log(x(0)) - log(x(1));
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [&](const anyprog::real_block& x) {
            return x(0) - x(1);
        }
    };

    std::vector<anyprog::optimization::equation_condition_function_t> eq = {
        [&](const anyprog::real_block& x) {
            return x(0) + 2 * x(1) - 5;
        }
    };

    anyprog::optimization::range_t range = { 0, 10 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(ineq);
    opt.set_equation_condition(eq);
    auto ret = opt.search(10,3);
    anyprog::print(opt.is_ok(), ret, obj);
    return 0;
}