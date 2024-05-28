#include <vector>
#include <cmath>

class LatticeBoltzmann {
private:
    int nx, ny;
    std::vector<std::vector<std::vector<double>>> f, feq, ftemp;
    const double tau = 0.6;  // Relaxation time
    const double w[9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};  // Weights for D2Q9
    const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};  // X velocities for D2Q9
    const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};  // Y velocities for D2Q9

public:
    LatticeBoltzmann(int nx, int ny) : nx(nx), ny(ny), f(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9))),
                                       feq(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9))),
                                       ftemp(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9))) {}

    void initialize() {
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int i = 0; i < 9; i++) {
                    f[x][y][i] = w[i];  // Initialize distribution function
                }
            }
        }
    }

    void collideAndStream() {
        // Collide
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                double rho = 0, ux = 0, uy = 0;
                for (int i = 0; i < 9; i++) {
                    rho += f[x][y][i];
                    ux += f[x][y][i] * cx[i];
                    uy += f[x][y][i] * cy[i];
                }
                ux /= rho; uy /= rho;
                for (int i = 0; i < 9; i++) {
                    double cu = 3 * (cx[i] * ux + cy[i] * uy);
                    feq[x][y][i] = rho * w[i] * (1 + cu + 0.5 * cu * cu - 1.5 * (ux * ux + uy * uy));
                    ftemp[x][y][i] = f[x][y][i] - (f[x][y][i] - feq[x][y][i]) / tau;
                }
            }
        }

        // Stream
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int i = 0; i < 9; i++) {
                    int nx = (x + cx[i] + nx) % nx;
                    int ny = (y + cy[i] + ny) % ny;
                    f[nx][ny][i] = ftemp[x][y][i];
                }
            }
        }
    }
};
