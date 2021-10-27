/*
 * @file Timer.cpp
 *
 * @date Oct 22, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Timer.hpp"

#include <chrono>
#include <string>
#include <sstream>
#include <map>
#include <ctime>

namespace Nextsim {
// Static main clock
Timer Timer::main("main");

Timer::Timer()
    : Timer("")
{}

Timer::Timer(const Key& baseTimerName)
    : root()
    , current(&root)
{
    root.name = baseTimerName;
}

void Timer::tick(const Timer::Key& timerName)
{
    // Descend to the next child node, creating it if necessary
    TimerNode* parent = current;
    current = &current->childNodes[timerName];
    // redundant assignment, except when a new node has been created above.
    current->name  = timerName;
    current->parent = parent;

    // Mark the time, ticks and run status
    current->wallHack = std::chrono::high_resolution_clock::now();
    current->cpuHack = std::clock();
    ++current->ticks;
    current->running = true;
}

void Timer::tock(const std::string& timerName)
{
    // Calculate the durations and stop running
    current->wallTime += std::chrono::duration_cast<WallTimeDuration>(
            std::chrono::high_resolution_clock::now() - current->wallHack);
    current->cpuTime += static_cast<double>(std::clock() - current->cpuHack) / CLOCKS_PER_SEC;
    current->running = false;

    // Ascend to the parent
    current = current->parent;
}

//TODO: implement lap and elapsed
double Timer::lap(const Key& timerName) const { return 0; }
double Timer::elapsed(const Key& timerName) const { return 0; }


void Timer::additionalTime(
    const TimerPath& path, WallTimeDuration wallAdd, CpuTimeDuration cpuAdd, int ticksAdd)
{
    // Descend the given path
    TimerNode& cursor = root;
    for (auto& nodeName : path) {
        cursor = cursor.childNodes[nodeName];
    }
    cursor.wallTime += wallAdd;
    cursor.cpuTime += cpuAdd;
    cursor.ticks += ticksAdd;
}

Timer::TimerPath Timer::currentTimerNodePath() const
{
    TimerPath path;
    TimerNode* cursor = current;
    while (cursor != &root) {
        path.push_front(cursor->name);
        cursor = cursor->parent;
    }
    return path;
}

Timer::TimerPath Timer::pathToFirstMatch(const Key& timerName) const
{
    return root.searchDescendants(timerName);
}

std::ostream& Timer::report(const Key& timerName, std::ostream& os) const
{
    return report(pathToFirstMatch(timerName), os);
}

std::ostream& Timer::report(std::ostream& os) const
{
    return root.reportAll(os, "");
}

std::ostream& Timer::report(const TimerPath& path, std::ostream& os) const
{
    const TimerNode* cursor = &root;
    for (auto& element: path) {
        cursor = &cursor->childNodes.at(element);
    }
    return cursor->report(os, "");
}

Timer::TimerNode::TimerNode()
    : wallTime(WallTimeDuration::zero())
    , cpuTime(0)
    , ticks(0)
    , running(false)
    , parent(nullptr)
{ }

Timer::TimerPath Timer::TimerNode::searchDescendants(const Key& timerName) const
{
    TimerPath path;
    for (auto& children: childNodes) {
        if (children.first == timerName) {
            path.push_front(children.first);
        } else if (!children.second.searchDescendants(timerName).empty()) {
            path.push_front(name);
        }
    }
    return path;
}

std::ostream& Timer::TimerNode::report(std::ostream& os, const std::string& prefix) const
{
    return os << name << ": ticks = " << ticks << " wall time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(wallTime).count() * 1e-6
              << " s cpu time = " << cpuTime << " s" << std::endl;
}

std::ostream& Timer::TimerNode::reportAll(std::ostream& os, const std::string& prefix) const
{
    report(os, prefix);
    for (auto& child: childNodes) {
        child.second.reportAll(os, " " + prefix);
    }
    return os;
}
/*    std::stringstream builder;
    builder << timerName << ": is " << (running[timerName] ? "" : "not ")
            << wallTime[timerName] << " s wall time, "
            << cpuTime[timerName] << " s CPU time,"
            << calls[timerName] << " timed calls"
            << std::endl;

    return builder.str();
}*/
}

std::ostream& operator<<(std::ostream& os, const Nextsim::Timer& tim)
{
    return tim.report(os);
}
