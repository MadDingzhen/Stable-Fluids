#include "GridCells2D.hpp"

// 构造函数，初始化网格单元和 LBM 所需的数据结构
GridCells2D::GridCells2D() 
    : u(), v(), u0(), v0(), fx(), fy(), dens(),
      // 初始化三维向量，用于存储 LBM 的分布函数
      f(N, std::vector<std::vector<double>>(N, std::vector<double>(9))),
      feq(N, std::vector<std::vector<double>>(N, std::vector<double>(9))),
      ftemp(N, std::vector<std::vector<double>>(N, std::vector<double>(9)))
{
    // 初始化 LBM 的分布函数
    initializeLBM(N, N);
    for (int i = 0; i < SIZE; ++i)
    {
        // 初始颜色设为蓝色
        initialColor[i][0] = 0.0f;
        initialColor[i][1] = 0.0f;
        initialColor[i][2] = 1.0f;
    }
}

// 初始化 LBM 的分布函数
void GridCells2D::initializeLBM(int nx, int ny)
{
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int i = 0; i < 9; i++) {
                f[x][y][i] = 1.0 / 9.0;  // 将每个格点的分布函数初始化为 1/9
            }
        }
    }
}

// 析构函数
GridCells2D::~GridCells2D()
{
}
