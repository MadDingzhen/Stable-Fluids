#pragma once
#include <vector>
#include "constants.hpp"

class GridCells2D
{
public:
    GridCells2D();
    ~GridCells2D();

    using Matrix3D = std::vector<std::vector<std::vector<double>>>;

    float u[SIZE], v[SIZE];
    float u0[SIZE], v0[SIZE];
    float fx[SIZE], fy[SIZE];
    float dens[SIZE];
    float initialColor[SIZE][3];  // 新增：记录每个网格点的初始颜色

    // LBM相关的分布函数
    Matrix3D f, feq, ftemp;

    void initializeLBM(int nx, int ny);
};
