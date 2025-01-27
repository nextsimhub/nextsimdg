/*!
 * @file stackTraceException.cpp
 * @date 27 Jan 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef STACKTRACEEXCEPTION_HPP
#define STACKTRACEEXCEPTION_HPP

#include <boost/exception/all.hpp>
#include <boost/stacktrace.hpp>

typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;

template <class E> void throw_with_trace(const E& e)
{
    throw boost::enable_error_info(e) << traced(boost::stacktrace::stacktrace());
}

#endif // STACKTRACEEXCEPTION_HPP
