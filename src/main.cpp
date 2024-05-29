#include <iostream>
#include <GLFW/glfw3.h>
#include "GridCells2D.hpp"
#include "Scene2D.hpp"
#include "Simulator2D.hpp"
#include "time_utils.hpp"
#include "performance_monitor.hpp"

void setDensityMode(int argc, char *argv[], EMode *mode);

int main(int argc, char *argv[])
{
    EMode mode = E_Continuous;
    setDensityMode(argc, argv, &mode);

    GridCells2D *grid_cell = new GridCells2D();
    Scene2D scene(grid_cell);
    Simulator2D *simulator = new Simulator2D(grid_cell, mode);

    if (!glfwInit())
    {
        fprintf(stderr, "Initialization failed!\n");
        exit(EXIT_FAILURE);
    }

    GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, TITLE, nullptr, nullptr);

    if (!window)
    {
        fprintf(stderr, "Window creation failed!");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glfwSetMouseButtonCallback(window, simulator->mouseEvent);
    glfwSetCursorPosCallback(window, simulator->mouseMoveEvent);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    std::cout << "\n*** START SIMULATION ***\n";

    perf_monitor.start();
    while (!glfwWindowShouldClose(window))
    {
        TimeUtils::currentTime += DT;
        simulator->update();
        scene.update();
        scene.draw();

        perf_monitor.logFrame();  // 记录每一帧的性能数据

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    perf_monitor.stop();

    std::cout << "*** END ***\n\n";

    if (simulator) delete simulator;
    if (grid_cell) delete grid_cell;

    glfwTerminate();
    return 0;
}

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
