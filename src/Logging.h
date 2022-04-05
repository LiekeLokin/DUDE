/*
 * Logging.h
 *
 *  Created on: 20 jan. 2022
 *      Author: impor
 */

#ifndef LOGGING_H_
#define LOGGING_H_

#define BOOST_LOG_DYN_LINK 1

#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>

namespace logging  = boost::log;
namespace attrs    = boost::log::attributes;

namespace dude_log {
// Set attribute and return the new value
template<typename ValueType>
ValueType set_get_attrib(const char* name, ValueType value) {
	auto attr = logging::attribute_cast<attrs::mutable_constant<ValueType>>(logging::core::get()->get_thread_attributes()[name]);
	attr.set(value);
	return attr.get();
}

// Convert file path to only the filename
std::string path_to_filename(const std::string& path);

void init(const std::string& file_name, const std::string& file_sev, const std::string& console_sev);

extern logging::sources::severity_logger<logging::trivial::severity_level> lg;

} // dude_log

#define SHOW_VAR(var) #var " = " << var << "; "
#define SHOW_2VARS(var1, var2) SHOW_VAR(var1) << SHOW_VAR(var2)
#define SHOW_3VARS(var1, var2, var3) SHOW_2VARS(var1, var2) << SHOW_VAR(var3)
#define SHOW_4VARS(var1, var2, var3, var4) SHOW_3VARS(var1, var2, var3) << SHOW_VAR(var4)

// Our macro that includes severity, filename and line number
#define DUDE_LOG(sev) \
		BOOST_LOG_WITH_PARAMS( \
				(dude_log::lg), \
				(dude_log::set_get_attrib("File", dude_log::path_to_filename(__FILE__))) \
				(dude_log::set_get_attrib("Line", __LINE__)) \
				(dude_log::set_get_attrib("Tijd", tijd)) \
				(::boost::log::keywords::severity = (boost::log::trivial::sev)) \
		)

#endif /* LOGGING_H_ */
