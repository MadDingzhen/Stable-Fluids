#include <iostream>
#include <GLFW/glfw3.h>
#include "GridCells2D.hpp"
#include "Scene2D.hpp"
#include "Simulator2D.hpp"
#include "time_utils.hpp"  // 引入全局时间工具

// 设置密度模式的函数声明
void setDensityMode(int argc, char *argv[], EMode *mode);

int main(int argc, char *argv[])
{
    // 设置默认密度模式
    EMode mode = E_Continuous;
    setDensityMode(argc, argv, &mode);

    // 创建网格单元、场景和模拟器对象
    GridCells2D *grid_cell = new GridCells2D();
    Scene2D scene(grid_cell);
    Simulator2D *simulator = new Simulator2D(grid_cell, mode);

    // 初始化 OpenGL
    if (!glfwInit())
    {
        fprintf(stderr, "Initialization failed!\n");
        exit(EXIT_FAILURE);
    }

    // 创建窗口
    GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, TITLE, nullptr, nullptr);

    if (!window)
    {
        fprintf(stderr, "Window creation failed!");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    // 注册事件回调函数
    glfwSetMouseButtonCallback(window, simulator->mouseEvent);
    glfwSetCursorPosCallback(window, simulator->mouseMoveEvent);

    // 初始化场景
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    std::cout << "\n*** START SIMULATION ***\n";

    while (!glfwWindowShouldClose(window))
    {
        TimeUtils::currentTime += DT;  // 更新全局时间
        // std::cout << "Main loop time: " << TimeUtils::getCurrentTime() << std::endl;

        simulator->update();
        scene.update();  // 确保 Scene2D 的 update 方法被调用
        scene.draw();

        // 交换绘图缓冲区
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    std::cout << "*** END ***\n\n";

    // 清理资源
    if (simulator)
    {
        delete simulator;
    }
    if (grid_cell)
    {
        delete grid_cell;
    }

    glfwTerminate();
    return 0;
}

// 设置密度模式的函数实现
void setDensityMode(int argc, char *argv[], EMode *mode)
{
    if (argc > 2)
    {
        fprintf(stderr, "too many arguments\n");
        exit(EXIT_FAILURE);
    }
    else if (argc == 2)
    {
        char *p = argv[1];
        if (*p == '-')
        {
            p++;
            // 设置密度模式
            switch (*p)
            {
            case 'o':
                *mode = E_Once;
                break;
            case 'c':
                *mode = E_Continuous;
                break;
            default:
                fprintf(stderr, "invalid argument\n");
                exit(EXIT_FAILURE);
                break;
            }
        }
        else
        {
            fprintf(stderr, "invalid argument\n");
            exit(EXIT_FAILURE);
        }
    }
}
