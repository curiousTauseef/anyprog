#include "../help.hpp"

//http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page4679.htm
//The best minima as we all known: x∗ = (0.20573,7.0924,9.03663,0.20573) , f(x∗) = 2.21813.

int main(int argc, char** argv)
{
    anyprog::optimization::function_t obj = [](const anyprog::real_block& x) {
        return 1.10471 * x(0) * x(0) * x(1) + 0.04811 * x(2) * x(3) * (14.0 + x(1));
    };

    double P = 6000, L = 14, E = 30e+6, G = 12e+6, t_max = 13600, s_max = 30000, d_max = 0.25;

    auto M = [&](const anyprog::real_block& x) {
        return P * (L + x(1) / 2.0);
    };
    auto R = [&](const anyprog::real_block& x) {
        return sqrt(0.25 * (pow(x(1), 2) + pow(x(0) + x(2), 2)));
    };
    auto J = [&](const anyprog::real_block& x) {
        return 2.0 / sqrt(2.0) * x(0) * x(1) * (x(1) * x(1) / 12.0 + 0.25 * pow(x(0) + x(2), 2));
    };
    auto P_c = [&](const anyprog::real_block& x) {
        return (4.013 * E / (6 * L * L)) * x(2) * pow(x(3), 3) * (1 - 0.25 * x(2) * sqrt(E / G) / L);
    };
    auto t1 = [&](const anyprog::real_block& x) {
        return P / (sqrt(2.0) * x(0) * x(1));
    };
    auto t2 = [&](const anyprog::real_block& x) {
        return M(x) * R(x) / J(x);
    };
    auto t = [&](const anyprog::real_block& x) {
        return sqrt(pow(t1(x), 2) + t1(x) * t2(x) * x(1) / R(x) + pow(t2(x), 2));
    };
    auto s = [&](const anyprog::real_block& x) {
        return 6 * P * L / (x(3) * x(2) * x(2));
    };
    auto d = [&](const anyprog::real_block& x) {
        return 4 * P * pow(L, 3) / (E * x(3) * pow(x(2), 3));
    };

    std::vector<anyprog::optimization::inequation_condition_function_t> ineq = {
        [&](const anyprog::real_block& x) {
            return t(x) - t_max;
        },
        [&](const anyprog::real_block& x) {
            return s(x) - s_max;
        },
        [&](const anyprog::real_block& x) {
            return x(0) - x(3);
        },
        [&](const anyprog::real_block& x) {
            return 0.10471 * x(0) * x(0) + 0.04811 * x(2) * x(3) * (14.0 + x(1)) - 5.0;
        },
        [&](const anyprog::real_block& x) {
            return d(x) - d_max;
        },
        [&](const anyprog::real_block& x) {
            return P - P_c(x);
        }
    };

    std::vector<anyprog::optimization::range_t> range = {
        { 0.125, 2 }, { 0.1, 10 }, { 0.1, 10 }, { 0.1, 10 }
    };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(10, 3);
    anyprog::print(opt.is_ok(), ret, obj);

    return 0;
}