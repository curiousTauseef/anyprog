#include "block.hpp"

namespace anyprog {
namespace block {

    real_block reshape(const real_block& ret, size_t rows, size_t cols)
    {
        return Eigen::Map<const Eigen::MatrixXd>(ret.data(), rows, cols);
    }
    complex_block reshape(const complex_block& ret, size_t rows, size_t cols)
    {
        return Eigen::Map<const Eigen::MatrixXcd>(ret.data(), rows, cols);
    }
}
}