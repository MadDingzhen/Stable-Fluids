#include "GridCells2D.hpp"

GridCells2D::GridCells2D() 
    : u(), v(), u0(), v0(), fx(), fy(), dens(),
      f(N, std::vector<std::vector<double>>(N, std::vector<double>(9))),
      feq(N, std::vector<std::vector<double>>(N, std::vector<double>(9))),
      ftemp(N, std::vector<std::vector<double>>(N, std::vector<double>(9)))
{
    initializeLBM(N, N);
}

void GridCells2D::initializeLBM(int nx, int ny)
{
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int i = 0; i < 9; i++) {
                f[x][y][i] = 1.0 / 9.0;
            }
        }
    }
}

GridCells2D::~GridCells2D()
{
}
