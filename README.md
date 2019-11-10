# anyprog
A C++ scientific library for mathematical programming,data fitting and solving nonlinear equations

# feature
- [anyprog](#anyprog)
- [feature](#feature)
- [usage](#usage)
  - [mathematical programming](#mathematical-programming)
    - [nonlinear-unconstrained-optimization](#nonlinear-unconstrained-optimization)
      - [example-1](#example-1)
      - [example-2](#example-2)
      - [example-3](#example-3)
    - [nonlinear-constrained-optimization](#nonlinear-constrained-optimization)
      - [example-1](#example-1-1)
      - [example-2](#example-2-1)
      - [example-3](#example-3-1)
      - [example-4](#example-4)
    - [linear-optimization](#linear-optimization)
    - [quadratic-optimazition](#quadratic-optimazition)
    - [mixed-integer-optimazition](#mixed-integer-optimazition)
      - [example-1](#example-1-2)
      - [example-2](#example-2-2)
      - [example-3](#example-3-2)
    - [any-optimazition](#any-optimazition)
      - [example-1](#example-1-3)
      - [example-2](#example-2-3)
      - [example-3](#example-3-3)
      - [example-4](#example-4-1)
      - [example-5](#example-5)
  - [data fitting](#data-fitting)
    - [polynomial-fitting](#polynomial-fitting)
    - [nonlinear-fitting](#nonlinear-fitting)
  - [solving nonlinear equations](#solving-nonlinear-equations)
    - [example-1](#example-1-4)
    - [example-2](#example-2-4)


# usage
`pkg-config --libs --cflags anyprog`

## mathematical programming

### nonlinear-unconstrained-optimization
#### example-1
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    auto output = [](const anyprog::real_block& ret, const anyprog::optimization::funcation_t& obj) {
        std::cout << "solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }
        std::cout << "object=\t" << obj(ret) << "\n\n";
    };

    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return 100 * pow(x(1, 0) - pow(x(0, 0), 2), 2) + pow(1 - x(0, 0), 2);
    };

    std::vector<anyprog::optimization::range_t> range = { { -10, 10 }, { -10, 10 } };
    output(anyprog::optimization::fminbnd(obj, range, 1e-10), obj);

    anyprog::real_block param(2, 1);
    param(0, 0) = 0;
    param(1, 0) = 0;
    output(anyprog::optimization::fminunc(obj, param, 1e-10), obj);

    auto grad = [&](const anyprog::real_block& x) {
        anyprog::real_block ret(x.rows(), 1);
        ret(0, 0) = 2 * x(0) - 400 * x(0) * (x(1) - pow(x(0), 2)) - 2;
        ret(1, 0) = 200 * x(1) - 200 * pow(x(0), 2);
        return ret;
    };
    output(anyprog::optimization::fminunc(obj, grad, param), obj);

    return 0;
}

```
```txt
solution:
x(0)=	0.999986
x(1)=	0.999971
object=	2.16042e-10

solution:
x(0)=	1
x(1)=	1
object=	3.05615e-14

solution:
x(0)=	1
x(1)=	1
object=	3.56293e-17

```

#### example-2
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return 100 * pow(x(1, 0) - pow(x(0, 0), 2), 2) + pow(1 - x(0, 0), 2);
    };
    std::vector<anyprog::optimization::range_t> range = { { -10, 10 }, { -10, 10 } };
    anyprog::optimization opt(obj, range);
    auto ret = opt.search();

    std::cout << "global solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n\n";

    auto history = opt.get_history();
    std::cout << "search history:\n";
    for (auto& iter : history) {
        for (size_t i = 0; i < iter.second.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << iter.second(i, 0) << "\n";
        }
        std::cout << "object=\t" << iter.first << "\n\n";
    }

    return 0;
}
```
```txt
global solution:
x(0)=	1.0025
x(1)=	1.00506
object=	6.45722e-06

search history:
x(0)=	0.845278
x(1)=	0.713799
object=	0.0239873

x(0)=	1.09488
x(1)=	1.19905
object=	0.00901008

x(0)=	1.05482
x(1)=	1.11285
object=	0.00300943

x(0)=	0.963414
x(1)=	0.928027
object=	0.00134049

x(0)=	1.02463
x(1)=	1.05004
object=	0.000609679

x(0)=	1.02144
x(1)=	1.04345
object=	0.000460868

x(0)=	1.0025
x(1)=	1.00506
object=	6.45722e-06

```
#### example-3
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return x(0, 0) * sin(x(0, 0)) * cos(2.0 * x(0, 0)) - 2.0 * x(0, 0) * sin(3.0 * x(0, 0));
    };

    std::vector<anyprog::optimization::range_t> range = { { 0, 20 } };

    anyprog::optimization opt(obj, range);
    auto ret = opt.search(100, 10);

    std::cout << "global solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n\n";

    auto history = opt.get_history();
    std::cout << "search history:\n";
    for (auto& iter : history) {
        for (size_t i = 0; i < iter.second.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << iter.second(i, 0) << "\n";
        }
        std::cout << "object=\t" << iter.first << "\n\n";
    }

    return 0;
}
```

```txt
global solution:
x(0)=	19.4114
object=	-34.0963

search history:
x(0)=	15.1613
object=	-26.6283

x(0)=	19.4114
object=	-34.0963
```

### nonlinear-constrained-optimization
#### example-1
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return -log(x(0)) - log(x(1));
    };

    std::vector<anyprog::optimization::inequation_condition_funcation_t> ineq;
    ineq.emplace_back([&](const anyprog::real_block& x) {
        return x(0) - x(1);
    });

    std::vector<anyprog::optimization::equation_condition_funcation_t> eq;
    eq.emplace_back([&](const anyprog::real_block& x) {
        return x(0) + 2 * x(1) - 5;
    });

    std::vector<anyprog::optimization::range_t> range = { { 0, 10 }, { 0, 10 } };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    opt.set_equation_condition(eq);
    auto ret = opt.solve();

    std::cout << "solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n";

    return 0;
}
```
```txt
solution:
x(0)=	1.66667
x(1)=	1.66667
object=	-1.02165

```
#### example-2
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return pow(1 - x(0), 2) + 100 * pow(x(1) - pow(x(0), 2), 2);
    };

    std::vector<anyprog::optimization::inequation_condition_funcation_t> ineq;
    ineq.emplace_back([&](const anyprog::real_block& x) {
        return pow(x(0) - 1. / 3., 2) + pow(x(1) - 1. / 3., 2) - pow(1. / 3., 2);
    });

    ineq.emplace_back([&](const anyprog::real_block& x) {
        return x(0) + x(1) - 2;
    });

    ineq.emplace_back([&](const anyprog::real_block& x) {
        return x(0) - 2 * x(1) - 3;
    });

    std::vector<anyprog::optimization::range_t> range = { { 0, 0.5 }, { 0.2, 0.8 } };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.solve();

    std::cout << "solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n";
    return 0;
}
```
```txt
solution:
x(0)=	0.5
x(1)=	0.250026
object=	0.25

```
#### example-3
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    size_t dim = 6;
    anyprog::optimization::funcation_t obj = [&](const anyprog::real_block& x) {
        return sqrt(pow(x(0) - x(1), 2) + pow(x(3) - x(4), 2)) + sqrt(pow(x(0) - x(2), 2) + pow(x(3) - x(5), 2)) + sqrt(pow(x(2) - x(1), 2) + pow(x(5) - x(4), 2));
    };
    std::vector<anyprog::optimization::inequation_condition_funcation_t> ineq;
    ineq.emplace_back([&](const anyprog::real_block& x) {
        return pow(x(0) - 5, 2) + pow(x(3) - 4, 2) - 4;
    });
    ineq.emplace_back([&](const anyprog::real_block& x) {
        return pow(x(1) + 5, 2) + pow(x(4) + 3, 2) - 1;
    });
    ineq.emplace_back([&](const anyprog::real_block& x) {
        return pow(x(2) + 1, 2) + pow(x(5) - 1, 2) - 1;
    });

    std::vector<anyprog::optimization::range_t> range;
    for (size_t i = 0; i < dim; ++i) {
        range.push_back({ -100, 100 });
    }
    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(ineq);
    auto ret = opt.search(100, 10);

    std::cout << "global solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n\n";

    auto history = opt.get_history();
    std::cout << "search history:\n";
    for (auto& iter : history) {
        for (size_t i = 0; i < iter.second.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << iter.second(i, 0) << "\n";
        }
        std::cout << "object=\t" << iter.first << "\n\n";
    }

    return 0;
}
```
```txt
global solution:
x(0)=	3.35966
x(1)=	-4.18107
x(2)=	-0.310501
x(3)=	2.85576
x(4)=	-2.4261
x(5)=	0.286241
object=	18.4131

search history:
x(0)=	3.36641
x(1)=	-4.18061
x(2)=	-0.441214
x(3)=	2.84614
x(4)=	-2.42676
x(5)=	0.184949
object=	18.4132

x(0)=	3.36015
x(1)=	-4.18141
x(2)=	-0.244818
x(3)=	2.85507
x(4)=	-2.42562
x(5)=	0.344485
object=	18.4131

x(0)=	3.35966
x(1)=	-4.18107
x(2)=	-0.310501
x(3)=	2.85576
x(4)=	-2.4261
x(5)=	0.286241
object=	18.4131

```
#### example-4
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    //https://www.coin-or.org/Ipopt/documentation/node20.html
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return x(0) * x(3) * (x(0) + x(1) + x(2)) + x(2);
    };
    std::vector<anyprog::optimization::inequation_condition_funcation_t> ineq;
    ineq.emplace_back([](const anyprog::real_block& x) {
        return 25 - x.prod();
    });
    std::vector<anyprog::optimization::equation_condition_funcation_t> eq;
    eq.emplace_back([](const anyprog::real_block& x) {
        size_t dim = x.rows();
        double sum = 0;
        for (size_t i = 0; i < dim; ++i) {
            sum += pow(x(i), 2);
        }
        return sum - 40;
    });
    std::vector<anyprog::optimization::range_t> range = { { 1, 5 }, { 1, 5 }, { 1, 5 }, { 1, 5 } };
    anyprog::real_block param(4, 1);
    param << 1, 5, 5, 1;
    anyprog::optimization opt(obj, param, range);
    opt.set_equation_condition(eq).set_inequation_condition(ineq);
    auto ret = opt.solve();
    if (opt.is_ok()) {
        std::cout << "solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }
        std::cout << "object:\t" << obj(ret) << "\n";
    } else {
        std::cout << "Not found.\n";
    }

    return 0;
}
```
```txt
solution:
x(0)=	1
x(1)=	4.74303
x(2)=	3.82112
x(3)=	1.37941
object:	17.014
```

### linear-optimization
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::real_block obj(3, 1);
    obj << -5, -4, -6;

    anyprog::real_block A(3, 3), b(3, 1);
    A << 1, -1, 1, 3, 2, 4, 3, 2, 0;
    b << 20, 42, 30;

    std::vector<anyprog::optimization::range_t> range = { { 0, 20 }, { 0, 20 }, { 0, 20 } };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(A, b);
    anyprog::real_block ret = opt.solve();

    std::cout << "solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n";

    return 0;
}
```
```txt
solution:
x(0)=	3.9968e-15
x(1)=	15
x(2)=	3
object=	-78
```

### quadratic-optimazition
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::real_block h(2, 2);
    h << 2, -1, -1, 4;
    anyprog::real_block c(1, 2);
    c << -1, -10;

    anyprog::optimization::funcation_t obj = [&](const anyprog::real_block& x) {
        anyprog::real_block ret = 0.5 * (x.transpose() * h * x) + c * x;
        return ret(0, 0);
    };

    anyprog::real_block A(1, 2);
    A << 3, 2;
    anyprog::real_block b(1, 1);
    b << 6;

    std::vector<anyprog::optimization::range_t> range = { { 0, 10 }, { 0, 10 } };

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(A, b);

    auto ret = opt.solve();

    std::cout << "solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n";
    return 0;
}
```
```txt
solution:
x(0)=	0.499887
x(1)=	2.25017
object=	-13.75

```

### mixed-integer-optimazition
#### example-1
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    size_t dim = 4;
    anyprog::real_block obj(dim, 1);
    obj << 3, 7, -1, 1;

    anyprog::real_block A(3, dim);
    A << -2, 1, -1, 1, -1, 1, -6, -4, -5, -3, 0, -1;

    anyprog::real_block b(3, 1);
    b << -1, -6, -5;

    std::vector<anyprog::optimization::range_t> range;
    for (size_t i = 0; i < dim; ++i) {
        range.push_back({ 0, 1 });
    }

    anyprog::optimization opt(obj, range);
    opt.set_inequation_condition(A, b)
        .set_enable_binary_filter();

    auto ret = opt.search(100, 10);

    std::cout << "global solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n\n";

    auto history = opt.get_history();
    std::cout << "search history:\n";
    for (auto& iter : history) {
        for (size_t i = 0; i < iter.second.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << iter.second(i, 0) << "\n";
        }
        std::cout << "object=\t" << iter.first << "\n\n";
    }

    return 0;
}
```
```txt
global solution:
x(0)=	1
x(1)=	0
x(2)=	1
x(3)=	0
object=	2

search history:
x(0)=	1
x(1)=	1
x(2)=	1
x(3)=	0
object=	9

x(0)=	1
x(1)=	0
x(2)=	1
x(3)=	1
object=	3

x(0)=	1
x(1)=	0
x(2)=	1
x(3)=	0
object=	2

```
#### example-2
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return pow(x(0) - 1, 2) + pow(x(1) - 1, 2) + pow(x(2) - 1, 2) - log(1 + x(3)) + pow(x(4) - 1, 2) + pow(x(5) - 2, 2) + pow(x(6) - 3, 2);
    };
    std::vector<anyprog::optimization::inequation_condition_funcation_t> ineq;
    ineq.emplace_back([](const anyprog::real_block& x) {
        return x.sum() - x(3) - 5;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return pow(x(2), 2) + pow(x(4), 2) + pow(x(5), 2) + pow(x(6), 2) - 5.5;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return x(0) + x(4) - 1.2;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return x(1) + x(5) - 1.8;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return x(2) + x(6) - 2.5;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return x(3) + x(4) - 1.2;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return pow(x(1), 2) + pow(x(5), 2) - 1.64;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return pow(x(2), 2) + pow(x(6), 2) - 4.25;
    });
    ineq.emplace_back([](const anyprog::real_block& x) {
        return pow(x(1), 2) + pow(x(6), 2) - 4.64;
    });

    anyprog::real_block param(7, 1);
    param << 1, 1, 1, 1, 1, 1, 1;

    anyprog::optimization opt(obj, param);
    opt.set_inequation_condition(ineq);
    opt.set_filter_function([](anyprog::real_block& x) {
        x(0, 0) = round(x(0, 0));
        x(1, 0) = round(x(1, 0));
        x(2, 0) = round(x(2, 0));
        x(3, 0) = round(x(3, 0));
    });
    auto ret = opt.search();
    if (opt.is_ok()) {
        std::cout << "golbal solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }
        std::cout << "object:\t" << obj(ret) << "\n";
    } else {
        std::cout << "Not found.\n";
    }

    return 0;
}
```
```txt
golbal solution:
x(0)=	1
x(1)=	0
x(2)=	0
x(3)=	1
x(4)=	0.2
x(5)=	1.28062
x(6)=	1.95448
object:	3.55746
```
#### example-3
![mixed_integer_programming.png](doc/mixed_integer_programming.png)
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    size_t dim = 4;
    anyprog::real_block obj(dim, 1);
    obj << -1, -1, -2, 2;
    anyprog::real_block Aeq(1, dim), beq(1, 1), A(3, dim), b(3, 1);
    Aeq << 1, 1, 1, 1;
    beq << 10;
    A << 1, 0, 2, 0, 0, 2, -8, 0, 0, -1, 2, -1;
    b << 700, 0, -1;

    std::vector<anyprog::optimization::range_t> range;
    for (size_t i = 0; i < dim; ++i) {
        range.push_back({ 0, 10 });
    }
    anyprog::optimization opt(obj, range);
    opt.set_equation_condition(Aeq, beq).set_inequation_condition(A, b);
    opt.set_enable_integer_filter();
    anyprog::real_block ret = opt.search();
    if (opt.is_ok()) {
        std::cout << "global solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }
        std::cout << "object:\t" << ret.transpose() * obj << "\n\n";

        auto history = opt.get_history();
        std::cout << "search history:\n";
        for (auto& iter : history) {
            for (size_t i = 0; i < iter.second.rows(); ++i) {
                std::cout << "x(" << i << ")=\t" << iter.second(i, 0) << "\n";
            }
            std::cout << "object=\t" << iter.first << "\n\n";
        }
    } else {
        std::cout << "Not found.\n";
    }

    return 0;
}
```
```txt
global solution:
x(0)=	0
x(1)=	7
x(2)=	3
x(3)=	0
object:	-13

search history:
x(0)=	4
x(1)=	3
x(2)=	1
x(3)=	2
object=	-5

x(0)=	0
x(1)=	6
x(2)=	3
x(3)=	1
object=	-10

x(0)=	5
x(1)=	4
x(2)=	1
x(3)=	0
object=	-11

x(0)=	3
x(1)=	5
x(2)=	2
x(3)=	0
object=	-12

x(0)=	0
x(1)=	7
x(2)=	3
x(3)=	0
object=	-13
```


### any-optimazition
#### example-1
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::optimization::funcation_t obj = [](const anyprog::real_block& x) {
        return x(0) * x(3) * (x(0) + x(1) + x(2)) + x(2);
    };

    std::vector<anyprog::optimization::equation_condition_funcation_t> eq;
    eq.emplace_back([](const anyprog::real_block& x) {
        double sum = 0;
        for (size_t i = 0; i < 4; ++i) {
            sum += pow(x(i), 2);
        }
        return sum - 40;
    });

    std::vector<anyprog::optimization::inequation_condition_funcation_t> ineq;
    ineq.emplace_back([](const anyprog::real_block& x) {
        return 25 - x(0) * x(1) * x(2) * x(3);
    });

    std::vector<anyprog::optimization::range_t> range = { { 1, 5 }, { 1, 5 }, { 1, 5 }, { 1, 5 } };

    anyprog::optimization opt(obj, range);

    opt.set_equation_condition(eq)
        .set_inequation_condition(ineq)
        .set_filter_function([](anyprog::real_block& x) {
            x(0, 0) = round(x(0, 0));
        });

    auto ret = opt.search();
    std::cout << "global solution:\n";
    for (size_t i = 0; i < ret.rows(); ++i) {
        std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
    }
    std::cout << "object=\t" << opt.obj(ret) << "\n\n";

    auto history = opt.get_history();
    std::cout << "search history:\n";
    for (auto& iter : history) {
        for (size_t i = 0; i < iter.second.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << iter.second(i, 0) << "\n";
        }
        std::cout << "object=\t" << iter.first << "\n\n";
    }

    return 0;
}
```
```txt
global solution:
x(0)=	1
x(1)=	4.74309
x(2)=	3.82103
x(3)=	1.37942
object=	17.014

search history:
x(0)=	3
x(1)=	5
x(2)=	2.23607
x(3)=	1
object=	32.9443

x(0)=	2
x(1)=	3.16228
x(2)=	5
x(3)=	1
object=	25.3246

x(0)=	1
x(1)=	4.74309
x(2)=	3.82103
x(3)=	1.37942
object=	17.014
```
#### example-2
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    //https://uk.mathworks.com/help/optim/ug/mixed-integer-linear-programming-basics.html
    anyprog::real_block obj(8, 1);
    obj << 350 * 5, 330 * 3, 310 * 4, 280 * 6, 500, 450, 400, 100;
    anyprog::real_block Aeq(3, 8), beq(3, 1);
    Aeq << 5, 3, 4, 6, 1, 1, 1, 1,
        5 * 0.05, 3 * 0.04, 4 * 0.05, 6 * 0.03, 0.08, 0.07, 0.06, 0.03,
        5 * 0.03, 3 * 0.03, 4 * 0.04, 6 * 0.04, 0.06, 0.07, 0.08, 0.09;
    beq << 25, 1.25, 1.25;
    std::vector<anyprog::optimization::range_t> range = { { 0, 1 }, { 0, 1 }, { 0, 1 }, { 0, 1 }, { 0, 10 }, { 0, 10 }, { 0, 10 }, { 0, 10 } };
    anyprog::optimization opt(obj, range);
    opt.set_equation_condition(Aeq, beq);
    opt.set_filter_function([](anyprog::real_block& x) {
        x(0, 0) = round(x(0, 0));
        x(1, 0) = round(x(1, 0));
        x(2, 0) = round(x(2, 0));
        x(3, 0) = round(x(3, 0));
    });
    auto ret = opt.search();
    if (opt.is_ok()) {
        std::cout << "solution:\n"
                  << ret << "\n";
        std::cout << "object:\t" << ret.transpose() * obj << "\n";
    } else {
        std::cout << "Not found.\n";
    }

    return 0;
}
```
```txt
solution:
       1
       1
       0
       1
 7.13349
0.233014
0.133493
     3.5
object:	8495

```
#### example-3
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    //https://www.geeksforgeeks.org/job-assignment-problem-using-branch-and-bound/
    anyprog::real_block c(4, 4);
    c << 9, 2, 7, 8,
        6, 4, 3, 7,
        5, 8, 1, 8,
        7, 6, 9, 4;
    bool ok = false;
    double fval = 0;
    auto ret = anyprog::optimization::assignment(c, ok, fval);
    if (ok) {
        std::cout << fval << "\n";
        std::cout << ret << "\n";
    } else {
        std::cout << "Not Found\n";
    }
    return 0;
}
```
```txt
13
0 1 0 0
1 0 0 0
0 0 1 0
0 0 0 1
```
#### example-4
![tsp.png](doc/tsp.png)
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>


int main(int argc, char** argv)
{
    //Travelling Salesman Problem
    anyprog::real_block c(4, 4);
    double inf = 10000;
    c << inf, 500, 600, 100,
        100, inf, 800, 500,
        1000, 200, inf, 2000,
        400, 400, 100, inf;
    std::cout << c << "\n\n";
    anyprog::optimization::tsp tsp(c);
    auto path = tsp.solve();
    std::cout << "path=|";
    for (auto& i : path) {
        std::cout << i << "|";
    }
    std::cout << "\ndistance="
              << tsp.obj() << "\n";
    return 0;
}
```
```txt
10000   500   600   100
  100 10000   800   500
 1000   200 10000  2000
  400   400   100 10000

path=|0|3|2|1|0|
distance=500
```

#### example-5
![tsp2.png](doc/tsp2.png)
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    //https://blog.csdn.net/haohaoliu_/article/details/83024664
    std::vector<std::string> city = {
        "沈阳市",
        "长春市",
        "哈尔滨市",
        "北京市",
        "天津市",
        "呼和浩特市",
        "银川市",
        "太原市",
        "石家庄市",
        "济南市",
        "郑州市",
        "西安市",
        "武汉市",
        "南京市",
        "合肥市",
        "上海市",
        "长沙市",
        "南昌市",
        "杭州市",
        "福州市",
        "广州市",
        "台北市",
        "海口市",
        "南宁市",
        "重庆市",
        "昆明市",
        "贵阳市",
        "成都市",
        "兰州市",
        "西宁市",
        "拉萨市",
        "乌鲁木齐市",
        "香港",
        "澳门"
    };

    std::vector<std::pair<double, double>> gps = {
        { 123.429092, 41.796768 },
        { 125.324501, 43.886841 },
        { 126.642464, 45.756966 },
        { 116.405289, 39.904987 },
        { 117.190186, 39.125595 },
        { 111.751990, 40.841490 },
        { 106.232480, 38.486440 },
        { 112.549248, 37.857014 },
        { 114.502464, 38.045475 },
        { 117.000923, 36.675808 },
        { 113.665413, 34.757977 },
        { 108.948021, 34.263161 },
        { 114.298569, 30.584354 },
        { 118.76741, 32.041546 },
        { 117.283043, 31.861191 },
        { 121.472641, 31.231707 },
        { 112.982277, 28.19409 },
        { 115.892151, 28.676493 },
        { 120.15358, 30.287458 },
        { 119.306236, 26.075302 },
        { 113.28064, 23.125177 },
        { 121.5200760, 25.0307240 },
        { 110.199890, 20.044220 },
        { 108.320007, 22.82402 },
        { 106.504959, 29.533155 },
        { 102.71225, 25.040609 },
        { 106.713478, 26.578342 },
        { 104.065735, 30.659462 },
        { 103.834170, 36.061380 },
        { 101.777820, 36.617290 },
        { 91.11450, 29.644150 },
        { 87.616880, 43.826630 },
        { 114.165460, 22.275340 },
        { 113.549130, 22.198750 }
    };
    // https://github.com/janantala/GPS-distance
    // http://www.movable-type.co.uk/scripts/latlong-vincenty.html
    auto dis = anyprog::optimization::tsp::gps_distance(gps);
    anyprog::optimization::tsp tsp(dis);
    auto path = tsp.solve();
    std::cout << "path=|";
    for (auto& i : path) {
        std::cout << city[i] << "|";
    }
    std::cout << "\ndistance="
              << tsp.obj() << " km\n";
    return 0;
}
```
```txt
path=|沈阳市|长春市|哈尔滨市|上海市|杭州市|南京市|合肥市|南昌市|武汉市|长沙市|广州市|澳门|香港|海口市|南宁市|贵阳市|重庆市|成都市|兰州市|西宁市|昆明市|银川市|西安市|呼和浩特市|太原市|郑州市|石家庄市|北京市|天津市|济南市|福州市|台北市|拉萨市|乌鲁木齐市|沈阳市|
distance=47037 km
```

## data fitting
plotting data with python3:
```python

import matplotlib.pyplot as plt
import numpy as np

x, y, fit = np.loadtxt('dat.csv',delimiter=',', unpack=True)
plt.plot(x, y, 'ro', label='dat')
plt.plot(x, fit, 'b-', label='fit')
plt.xlabel('x')
plt.ylabel('y')
plt.title('anyprog data fitting')
plt.legend()
plt.show()
```
### polynomial-fitting

```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::real_block x(11, 1), y(11, 1);
    x << 0, 0.3000, 0.6000, 0.9000, 1.2000, 1.5000, 1.8000, 2.1000, 2.4000, 2.7000, 3.0000;
    y << 2.0000, 2.3780, 3.9440, 7.3460, 13.2320, 22.2500, 35.0480, 52.2740, 74.5760, 102.6020, 137.0000;
    anyprog::fit fit(x, 3);
    auto ret = fit.solve(y), fitting = fit.fitting(ret);
    std::cout << ret << std::endl;
    std::ofstream out("dat.csv", std::ios::trunc | std::ios::ate);
    for (size_t i = 0; i < x.rows(); ++i) {
        out << x(i) << "," << y(i) << "," << fitting(i) << "\n";
    }
    return 0;
}
```
```txt
           4
           3
-1.06581e-14
           2
```
![polyfitting](doc/polyfitting.png)


### nonlinear-fitting
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    anyprog::real_block x(18, 1), y(18, 1);
    x << 0, 0.4, 1.2, 2, 2.8, 3.6, 4.4, 5.2, 6, 7.2, 8, 9.2, 10.4, 11.6, 12.4, 13.6, 14.4, 15;
    y << 1, 0.85, 0.29, -0.27, -0.53, -0.4, -0.12, 0.17, 0.28, 0.15, -0.03, -0.15, -0.071, 0.059, 0.08, 0.032, -0.015, -0.02;
    anyprog::real_block param(3, 1);
    param(0, 0) = 1;
    param(1, 0) = 1;
    param(2, 0) = -0.1;
    std::vector<anyprog::fit::funcation_t> fun;
    fun.emplace_back([](const anyprog::real_block& x, const anyprog::real_block& p) {
        return (cos(p(1, 0) * x(0, 0)) * pow(M_E, p(2, 0) * x(0, 0))) / p(0, 0);
    });
    anyprog::fit fit(x, fun, param);
    auto ret = fit.lssolve(y), fitting = fit.fitting(ret);
    std::cout << ret << std::endl;
    std::ofstream out("dat.csv", std::ios::trunc | std::ios::ate);
    for (size_t i = 0; i < x.rows(); ++i) {
        out << x(i) << "," << y(i) << "," << fitting(i) << "\n";
    }
    return 0;
}
```
```txt
   0.9375
  1.00111
-0.207208
```
![nlfitting](doc/nlfitting.png)

## solving nonlinear equations
### example-1
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    std::vector<anyprog::equation::funcation_t> eq;
    eq.emplace_back([](const anyprog::real_block& x){
        return exp(-exp(x(0)+x(1)))-x(1)*(1+x(0)*x(0));
    });
    eq.emplace_back([](const anyprog::real_block&x){
        return x(0)*cos(x(1))+x(1)*sin(x(0))-0.5;
    });
    anyprog::real_block param(2,1);
    param<<0,0;
    anyprog::equation solver(eq,param);
    auto ret = solver.solve();
    if (solver.is_ok()) {
        std::cout << "solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }
    } else {
        std::cout << "Not found.\n";
    }

    return 0;
}
```
```txt
solution:
x(0)=	0.444489
x(1)=	0.139079
```
### example-2
```cpp
#include <anyprog/anyprog.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    std::vector<anyprog::equation::funcation_t> eq;
    eq.emplace_back([](const anyprog::real_block& x) {
        return x(0) - 0.6 * sin(x(0)) - 0.3 * cos(x(1));
    });
    eq.emplace_back([](const anyprog::real_block& x) {
        return x(1) - 0.6 * cos(x(0)) + 0.3 * sin(x(1));
    });
    anyprog::real_block param(2, 1);
    param << 0.5, 0.5;
    anyprog::equation solver(eq, param);
    auto ret = solver.solve();
    if (solver.is_ok()) {
        std::cout << "solution:\n";
        for (size_t i = 0; i < ret.rows(); ++i) {
            std::cout << "x(" << i << ")=\t" << ret(i, 0) << "\n";
        }
    } else {
        std::cout << "Not found.\n";
    }

    return 0;
}
```
```txt
solution:
x(0)=	0.635445
x(1)=	0.373439

```