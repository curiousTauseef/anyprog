#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page5161.htm
//http://ceur-ws.org/Vol-2255/paper2.pdf
//The global minima: x∗ = (0.051689156131,0.356720026419, 11.288831695483) , f(x∗) = 0.0126652327883.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return x(1) * x(0) * x(0) * (x(2) + 2);
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [](const anyprog::real_block& x) {
            return 1 - (pow(x(1), 3) * x(2)) / (71785 * pow(x(0), 4));
        },
        [](const anyprog::real_block& x) {
            return (4 * x(1) * x(1) - x(0) * x(1)) / (12566 * (x(1) * pow(x(0), 3) - pow(x(0), 4))) + 1.0 / (5108 * x(0) * x(0)) - 1;
        },
        [](const anyprog::real_block& x) {
            return 1 - (140.45 * x(0)) / (pow(x(1), 2) * x(2));
        },
        [](const anyprog::real_block& x) {
            return (x(1) + x(0)) / 1.5 - 1;
        },
    };

    std::vector<anyprog::optimization::range_t> range = {
        { 0.05, 2 }, { 0.25, 1.3 }, { 2, 15 }
    };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10,3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}