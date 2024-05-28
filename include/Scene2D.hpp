#pragma once
#include <iostream>
#include <string>
#include <GLFW/glfw3.h>
#include "GridCells2D.hpp"

class Scene2D
{
public:
    Scene2D(GridCells2D *grid_cells);
    ~Scene2D();

    void update();
    void draw();
    void writeData();
    float getTime() const { return m_time; }  // 新增：获取时间的函数

private:
    void drawGrid();
    void drawVelocity();
    void drawDensity();
    void writeData_inVtiFormat();

    int m_file_num;
    GridCells2D *m_grid_cells;
    float m_time;  // 时间变量
};
