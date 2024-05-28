#pragma once
#include <vector>
#include "constants.hpp"

class GridCells2D
{
public:
    GridCells2D();
    ~GridCells2D();

    // 定义类型别名
    using Matrix3D = std::vector<std::vector<std::vector<double>>>;

    float u[SIZE], v[SIZE];
    float u0[SIZE], v0[SIZE];
    float fx[SIZE], fy[SIZE];
    float dens[SIZE];

    // LBM相关的分布函数
    Matrix3D f, feq, ftemp;

    void initializeLBM(int nx, int ny);
};
