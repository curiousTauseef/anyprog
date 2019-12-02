#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2424.htm
//The global minima: x*=(2.171996, 2.363683, 8.773926, 5.095984, 0.9906548, 1.430574,1.321644, 9.828726, 8.280092, 8.375927),f(x*) = 24.3062091.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return x(0) * x(0) + x(1) * x(1) + x(0) * x(1) - 14 * x(0) - 16 * x(1) + (x(2) - 10) * (x(2) - 10) + 4 * (x(3) - 5) * (x(3) - 5) + (x(4) - 3) * (x(4) - 3) + 2 * (x(5) - 1) * (x(5) - 1) + 5 * x(6) * x(6) + 7 * (x(7) - 11) * (x(7) - 11) + 2 * (x(8) - 10) * (x(8) - 10) + (x(9) - 7) * (x(9) - 7) + 45;
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return 4 * x(0) + 5 * x(1) - 3 * x(6) + 9 * x(7) - 105;
        },
        [](const anyprog::real_block& x) {
            return 10 * x(0) - 8 * x(1) - 17 * x(6) + 2 * x(7);
        },
        [](const anyprog::real_block& x) {
            return -8 * x(0) + 2 * x(1) + 5 * x(8) - 2 * x(9) - 12;
        },
        [](const anyprog::real_block& x) {
            return 3 * (x(0) - 2) * (x(0) - 2) + 4 * (x(1) - 3) * (x(1) - 3) + 2 * x(2) * x(2) - 7 * x(3) - 120;
        },
        [](const anyprog::real_block& x) {
            return 5 * x(0) * x(0) + 8 * x(1) + (x(2) - 6) * (x(2) - 6) - 2 * x(3) - 40;
        },
        [](const anyprog::real_block& x) {
            return 0.5 * (x(0) - 8) * (x(0) - 8) + 2 * (x(1) - 4) * (x(1) - 4) + 3 * x(4) * x(4) - x(5) - 30;
        },
        [](const anyprog::real_block& x) {
            return x(0) * x(0) + 2 * (x(1) - 2) * (x(1) - 2) - 2 * x(0) * x(1) + 14 * x(4) - 6 * x(5);
        },
        [](const anyprog::real_block& x) {
            return -3 * x(0) + 6 * x(1) + 12 * (x(8) - 8) * (x(8) - 8) - 7 * x(9);
        }
    };

    size_t dim = 10;
    anyprog::optimization::range_t range = { -10, 10 };

    anyprog::optimization opt(obj, range, dim);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}