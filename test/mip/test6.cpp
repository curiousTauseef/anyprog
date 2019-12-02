#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 2;
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return pow(x(0), 4) + pow(x(1), 4) + 16 * (x(0) * x(1) + pow(4 + x(1), 2));
    };
    anyprog::optimization::range_t range = { -10, 10 };
    anyprog::optimization opt(obj, range, dim);
    opt.set_enable_integer_filter();
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);
}