#include <vector>
#include <cmath>

// 二维 D2Q9 模型的 Lattice Boltzmann 类
class LatticeBoltzmann {
private:
    int nx, ny;  // 网格尺寸
    std::vector<std::vector<std::vector<double>>> f, feq, ftemp;  // 分布函数和临时存储
    const double tau = 0.6;  // 松弛时间
    // D2Q9 模型的权重
    const double w[9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};
    // D2Q9 模型的 X 和 Y 方向速度
    const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

public:
    // 构造函数，初始化网格尺寸和分布函数
    LatticeBoltzmann(int nx, int ny) : nx(nx), ny(ny), 
        f(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9))),
        feq(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9))),
        ftemp(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9))) {}

    // 初始化分布函数
    void initialize() {
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int i = 0; i < 9; i++) {
                    f[x][y][i] = w[i];  // 初始化分布函数
                }
            }
        }
    }

    // 碰撞和传播步骤
    void collideAndStream() {
        // 碰撞步骤
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                double rho = 0, ux = 0, uy = 0;
                // 计算局部密度和速度
                for (int i = 0; i < 9; i++) {
                    rho += f[x][y][i];
                    ux += f[x][y][i] * cx[i];
                    uy += f[x][y][i] * cy[i];
                }
                ux /= rho; uy /= rho;
                // 计算平衡态分布函数
                for (int i = 0; i < 9; i++) {
                    double cu = 3 * (cx[i] * ux + cy[i] * uy);
                    feq[x][y][i] = rho * w[i] * (1 + cu + 0.5 * cu * cu - 1.5 * (ux * ux + uy * uy));
                    // 执行碰撞操作
                    ftemp[x][y][i] = f[x][y][i] - (f[x][y][i] - feq[x][y][i]) / tau;
                }
            }
        }

        // 传播步骤
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int i = 0; i < 9; i++) {
                    int nx = (x + cx[i] + nx) % nx;
                    int ny = (y + cy[i] + ny) % ny;
                    // 执行传播操作
                    f[nx][ny][i] = ftemp[x][y][i];
                }
            }
        }
    }
};
