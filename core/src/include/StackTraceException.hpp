/*!
 * @file stackTraceException.cpp
 * @date 30 Jan 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef STACKTRACEEXCEPTION_HPP
#define STACKTRACEEXCEPTION_HPP

// This is taken directly from the boost source: boost/stacktrace/detail/collect_unwind.ipp:33:2:
// version 1.76
#if defined(_GNU_SOURCE) || defined(BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED)                      \
    || defined(BOOST_WINDOWS)

#define STACKTRACEEXCEPTION_ACTIVE

#include <boost/exception/all.hpp>
#include <boost/stacktrace.hpp>

typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;

template <class E> void throw_with_trace(const E& e)
{
    throw boost::enable_error_info(e) << traced(boost::stacktrace::stacktrace());
}
#else
template <class E> void throw_with_trace(const E& e) { throw e; }
#endif

#endif // STACKTRACEEXCEPTION_HPP
