#include "block.hpp"
#include "random.hpp"

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
    real_block linspace(double l, double u, size_t len)
    {
        real_block ret(len, 1);
        double s = double(len - 1);
        for (size_t i = 0; i < len; ++i) {
            ret(i, 0) = l + (u - l) * i / s;
        }
        return ret;
    }
    complex_block linspace(const std::complex<double>& l, const std::complex<double>& u, size_t len)
    {
        complex_block ret(len, 1);
        double s = double(len - 1), re_l = l.real(), re_u = u.real(), im_l = l.imag(), im_u = u.imag();
        for (size_t i = 0; i < len; ++i) {
            ret(i, 0) = std::complex<double>(re_l + (re_u - re_l) * i / s, im_l + (im_u - im_l) * i / s);
        }
        return ret;
    }
}
}