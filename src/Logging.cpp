/*
 * Logging.cpp
 *
 *  Created on: 20 jan. 2022
 *      Author: impor
 */

#include "Logging.h"

#define BOOST_LOG_DYN_LINK 1

#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
#include <iomanip>

namespace logging  = boost::log;
namespace attrs    = boost::log::attributes;
namespace expr     = boost::log::expressions;
namespace keywords = boost::log::keywords;

namespace dude_log {

// Convert file path to only the filename
std::string path_to_filename(const std::string& path) {
	return path.substr(path.find_last_of("/\\")+1);
}

void init() {
	// Attributes that hold filename and line number
	logging::core::get()->add_thread_attribute("File", attrs::mutable_constant<std::string>(""));
	logging::core::get()->add_thread_attribute("Line", attrs::mutable_constant<int>(0));
	logging::core::get()->add_thread_attribute("Tijd", attrs::mutable_constant<double>(0));

#if 1
	logging::add_console_log(
			std::clog,
			keywords::filter = expr::attr< logging::trivial::severity_level >("Severity") >= logging::trivial::warning,
			keywords::format = (
					expr::stream
					<< "t=" << std::left << std::setw(5) << expr::attr<double>("Tijd")
					<< " <" << std::left << std::setw(7) << logging::trivial::severity << "> "
					<< expr::attr<std::string>("File")
					<< ':' << expr::attr<int>("Line") << " "
					<< expr::smessage
			)
	);
#endif
	logging::add_file_log(
			keywords::file_name = "sample.log",
			keywords::format = (
					expr::stream
					<< "t=" << std::left << std::setw(5) << expr::attr<double>("Tijd")
					<< " <" << std::left << std::setw(7) << logging::trivial::severity << "> "
					<< expr::attr<std::string>("File")
					<< ':' << expr::attr<int>("Line") << " "
					<< expr::smessage
			)
	);
}

logging::sources::severity_logger<logging::trivial::severity_level> lg;

} // dude_log
