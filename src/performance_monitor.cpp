#include "performance_monitor.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <sstream>
#include <sys/sysinfo.h>
#include <unistd.h>

PerformanceMonitor perf_monitor;

PerformanceMonitor::PerformanceMonitor() {}

PerformanceMonitor::~PerformanceMonitor() {
    printStats();
}

void PerformanceMonitor::start() {
    start_time = std::chrono::high_resolution_clock::now();
    last_frame_time = start_time;
}

void PerformanceMonitor::stop() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Total elapsed time: " << elapsed.count() << " seconds" << std::endl;
}

void PerformanceMonitor::logFrame() {
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> frame_time = now - last_frame_time;
    frame_times.push_back(frame_time.count());
    last_frame_time = now;

    logCpuUsage();
    logMemoryUsage();
    logBandwidthUsage();
}

void PerformanceMonitor::printStats() const {
    double avg_frame_time = 0;
    for (const auto &time : frame_times) {
        avg_frame_time += time;
    }
    avg_frame_time /= frame_times.size();

    std::cout << "Average frame time: " << avg_frame_time << " seconds" << std::endl;
    std::cout << "FPS: " << 1.0 / avg_frame_time << std::endl;
}

void PerformanceMonitor::logCpuUsage() const {
    std::ifstream stat_file("/proc/stat");
    std::string line;
    if (std::getline(stat_file, line)) {
        std::istringstream ss(line);
        std::string cpu;
        long user, nice, system, idle;
        ss >> cpu >> user >> nice >> system >> idle;
        long total = user + nice + system + idle;
        long usage = total - idle;
        std::cout << "CPU Usage: " << (usage * 100 / total) << "%" << std::endl;
    }
}

void PerformanceMonitor::logMemoryUsage() const {
    struct sysinfo mem_info;
    sysinfo(&mem_info);
    long total_memory = mem_info.totalram * mem_info.mem_unit;
    long free_memory = mem_info.freeram * mem_info.mem_unit;
    long used_memory = total_memory - free_memory;
    std::cout << "Memory Usage: " << (used_memory * 100 / total_memory) << "%" << std::endl;
}

void PerformanceMonitor::logBandwidthUsage() const {
    // Implement bandwidth usage logging
}
