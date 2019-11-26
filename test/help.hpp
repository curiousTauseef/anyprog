#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

namespace anyprog {
void print(bool ok, const real_block& ret, const optimization::function_t& obj)
{
    if (ok) {
        std::cout << "object=\t" << obj(ret) << "\n";
        std::cout << "solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }

    } else {
        std::cout << "Not Found.\n";
    }
}

void print(bool ok, const real_block& ret, const real_block& obj)
{
    if (ok) {
        std::cout << "object=\t" << obj.transpose() * ret << "\n";
        std::cout << "solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }

    } else {
        std::cout << "Not Found.\n";
    }
}
}