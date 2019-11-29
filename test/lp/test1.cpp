#include "../help.hpp"

int main(int argc, char** argv)
{
    size_t dim = 2;
    anyprog::real_block obj(dim, 1), A(3, dim), b(3, 1);
    obj << 8, 1;
    A << -1, -2,
        -4, -1,
        2, 1;
    b << 14, -33, 20;

    anyprog::optimization::range_t range = { 0, 100 };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(A, b);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);
    return 0;
}