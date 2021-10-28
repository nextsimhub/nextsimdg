/*
 * @file Timer.cpp
 *
 * @date Oct 22, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Timer.hpp"

#include "include/Chrono.hpp"
#include <chrono>
#include <ctime>
#include <map>
#include <regex>
#include <sstream>
#include <string>

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
    root.tick();
}

void Timer::tick(const Timer::Key& timerName)
{
    // Descend to the next child node, creating it if necessary
    TimerNode* parent = current;
    current = &current->childNodes[timerName];
    // redundant assignment, except when a new node has been created above.
    current->name  = timerName;
    current->parent = parent;

    // Mark the time, ticks and run status of the node
    current->tick();
}

void Timer::tock(const std::string& timerName)
{
    tock();
}

void Timer::tock()
{
    // Calculate the durations and stop running
    current->tock();
    // Ascend to the parent
    current = current->parent;
}

void Timer::TimerNode::tick()
{
    timeKeeper.start();
}

void Timer::TimerNode::tock()
{
    timeKeeper.stop();
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
    cursor.timeKeeper.extraWallTime(wallAdd);
    cursor.timeKeeper.extraCpuTime(cpuAdd);
    cursor.timeKeeper.extraTicks(ticksAdd);
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
    : parent(nullptr)
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
    os << prefix;
    // Get the wall time in seconds
    WallTimeDuration wallTimeNow = timeKeeper.wallTime();
    CpuTimeDuration cpuTimeNow = timeKeeper.cpuTime();
    os << name << ": ticks = " << timeKeeper.ticks();
    os << " wall time " << std::chrono::duration_cast<std::chrono::microseconds>(wallTimeNow).count() * 1e-6 << " s";
    os << " cpu time " << cpuTimeNow << " s";
    return os;
}

static std::string branch = "├";
static std::string spc = " ";
static std::string cont = "│";
static std::string last = "└";

std::string extendPrefix(const std::string& prefix)
{
    std::string newPrefix = std::regex_replace(prefix, std::regex(branch), cont);
    return std::regex_replace(newPrefix, std::regex(last), spc);
}

std::ostream& Timer::TimerNode::reportAll(std::ostream& os, const std::string& prefix) const
{
    report(os, prefix);
    os << std::endl;

    int nNodes = childNodes.size();
    int iNode = 0;

    for (auto& child: childNodes) {
        std::string lastBranch = (++iNode == nNodes) ? last : branch;
        child.second.reportAll(os, extendPrefix(prefix) + lastBranch);
    }
    return os;
}
}

std::ostream& operator<<(std::ostream& os, const Nextsim::Timer& tim)
{
    return tim.report(os);
}
