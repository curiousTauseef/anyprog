#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page4917.htm
//http://ijssst.info/Vol-17/No-43/paper5.pdf
//The best minima as we all known: x∗ = (0.73416509,0.36346997,38.03065892,234.73387898) , f(x∗) = 5821.19232.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return 0.6224 * x(0) * x(2) * x(3) + 1.7781 * x(1) * pow(x(2), 2) + 3.1661 * pow(x(0), 2) * x(3) + 19.84 * pow(x(0), 2) * x(2);
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return -x(0) + 0.0193 * x(2);
        },
        [](const anyprog::real_block& x) {
            return -x(1) + 0.00954 * x(2);
        },
        [](const anyprog::real_block& x) {
            return -M_PI * pow(x(2), 2) * x(3) - (4.0 / 3.0) * M_PI * pow(x(2), 3) + 1296000;
        },
        [](const anyprog::real_block& x) {
            return x(3) - 240;
        },
    };

    std::vector<anyprog::optimization::range_t> range = {
        { 0.0625, 200 }, { 0, 99 * 0.0625 }, { 10, 200 }, { 10, 250 }
    };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}