#pragma once
#include <chrono>
#include <fstream>
#include <vector>

class PerformanceMonitor {
public:
    PerformanceMonitor();
    ~PerformanceMonitor();

    void start();
    void stop();
    void logFrame();
    void printStats() const;

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> last_frame_time;
    std::vector<double> frame_times;

    void logCpuUsage() const;
    void logMemoryUsage() const;
    void logBandwidthUsage() const;
};

extern PerformanceMonitor perf_monitor;
