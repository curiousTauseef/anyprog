#include "../help.hpp"

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(x(0) - 1, 2) + pow(x(1) - 1, 2) + pow(x(2) - 1, 2) - log(1 + x(3)) + pow(x(4) - 1, 2) + pow(x(5) - 2, 2) + pow(x(6) - 3, 2);
    };
    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return x.sum() - x(3) - 5;
        },
        [](const anyprog::real_block& x) {
            return pow(x(2), 2) + pow(x(4), 2) + pow(x(5), 2) + pow(x(6), 2) - 5.5;
        },
        [](const anyprog::real_block& x) {
            return x(0) + x(4) - 1.2;
        },
        [](const anyprog::real_block& x) {
            return x(1) + x(5) - 1.8;
        },
        [](const anyprog::real_block& x) {
            return x(2) + x(6) - 2.5;
        },
        [](const anyprog::real_block& x) {
            return x(3) + x(4) - 1.2;
        },
        [](const anyprog::real_block& x) {
            return pow(x(1), 2) + pow(x(5), 2) - 1.64;
        },
        [](const anyprog::real_block& x) {
            return pow(x(2), 2) + pow(x(6), 2) - 4.25;
        },
        [](const anyprog::real_block& x) {
            return pow(x(1), 2) + pow(x(6), 2) - 4.64;
        }
    };

    std::vector<anyprog::optimization::range_t> range = { { 0, 1 }, { 0, 1 }, { 0, 1 }, { 0, 1 }, { 0, 10 }, { 0, 10 }, { 0, 10 } };
    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    opt.set_filter_function([](anyprog::real_block& x) {
        x(0, 0) = round(x(0, 0));
        x(1, 0) = round(x(1, 0));
        x(2, 0) = round(x(2, 0));
        x(3, 0) = round(x(3, 0));
    });
    auto ret = opt.search(10, 8);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}