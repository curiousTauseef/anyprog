#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 2;
    anyprog::real_block H(dim, dim), c(dim, 1);
    H << 1, -1,
        -1, 2;
    c << -2, -6;

    anyprog::optimization::function_t obj = [&](const anyprog::real_block& x) {
        anyprog::real_block ret = 0.5 * (x.transpose() * H * x) + c.transpose() * x;
        return ret(0, 0);
    };

    anyprog::real_block A(3, dim), b(3, 1);
    A << 1, 1,
        -1, 2,
        2, 1;
    b << 2, 2, 3;

    anyprog::optimization::range_t range = { 0, 100 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(A, b);
    auto ret = opt.search(10,3);
    anyprog::print(opt.is_ok(), ret, obj);
    return 0;
}