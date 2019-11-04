#ifndef ANYPROG_BLOCK_HPP
#define ANYPROG_BLOCK_HPP

#include "Eigen/Eigen"
#include <complex>
#include <vector>
#include <functional>

namespace anyprog {

typedef Eigen::MatrixXd real_block;
typedef Eigen::MatrixXcd complex_block;

namespace block{
    real_block reshape(const real_block& ret, size_t rows, size_t cols);
    complex_block reshape(const complex_block& ret, size_t rows, size_t cols);
}


}
#endif