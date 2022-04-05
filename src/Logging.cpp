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
#include <boost/algorithm/string.hpp>
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

logging::trivial::severity_level string_to_severity(const std::string& str) {
	const auto lc = boost::algorithm::to_lower_copy(str);
	if (lc == "fatal") return logging::trivial::fatal;
	if (lc == "error") return logging::trivial::error;
	if (lc == "warning") return logging::trivial::warning;
	if (lc == "info") return logging::trivial::info;
	if (lc == "debug") return logging::trivial::debug;

	// some warning here...
	return logging::trivial::debug;
}
void init(const std::string& file_name, const std::string& file_sev, const std::string& console_sev) {
	// Attributes that hold filename and line number
	logging::core::get()->add_thread_attribute("File", attrs::mutable_constant<std::string>(""));
	logging::core::get()->add_thread_attribute("Line", attrs::mutable_constant<int>(0));
	logging::core::get()->add_thread_attribute("Tijd", attrs::mutable_constant<double>(0));

#if 1
	logging::add_console_log(
			std::clog,
			keywords::filter = expr::attr< logging::trivial::severity_level >("Severity") >= string_to_severity(console_sev),
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
	if (! file_name.empty()) {
		logging::add_file_log(
				keywords::file_name = file_name,
				keywords::filter = expr::attr< logging::trivial::severity_level >("Severity") >= string_to_severity(file_sev),
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
}

logging::sources::severity_logger<logging::trivial::severity_level> lg;

} // dude_log
