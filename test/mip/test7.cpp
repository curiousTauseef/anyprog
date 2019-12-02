#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 3;
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        double sum = 0;
        for (size_t i = 0; i < x.rows(); ++i) {
            sum += pow(x(i), 4) - 4.9 * pow(x(i), 2);
        }
        return sum;
    };
    anyprog::optimization::range_t range = { -5, 5 };
    anyprog::optimization opt(obj, range, dim);
    opt.set_enable_integer_filter();
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);
}