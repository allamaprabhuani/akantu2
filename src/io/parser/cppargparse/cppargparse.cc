/**
 * @file   cppargparse.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Fri Mar 28 17:08:38 2014
 *
 * @brief  implementation of the ArgumentParser
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "cppargparse.hh"

#include <cstdlib>
#include <cstring>
#include <libgen.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <queue>
#include <sstream>
#include <algorithm>

#include <exception>
#include <stdexcept>

namespace cppargparse {

/* -------------------------------------------------------------------------- */
static inline std::string to_upper(const std::string & str) {
  std::string lstr = str;
  std::transform(lstr.begin(),
		 lstr.end(),
		 lstr.begin(),
		 (int(*)(int))std::toupper);
  return lstr;
}


/* -------------------------------------------------------------------------- */
/* ArgumentParser                                                             */
/* -------------------------------------------------------------------------- */
  ArgumentParser::ArgumentParser() : external_exit(NULL), prank(0), psize(1) {
  this->addArgument("-h;--help", "show this help message and exit", 0, _boolean, false, true);
}

ArgumentParser::~ArgumentParser() {
  for(_Arguments::iterator it = arguments.begin(); it != arguments.end(); ++it) {
    delete it->second;
  }
}

void ArgumentParser::setParallelContext(int prank, int psize) {
  this->prank = prank;
  this->psize = psize;
}

void ArgumentParser::_exit(const std::string & msg, int status) {
  if(prank == 0) {
    if(msg != "") {
      std::cerr << msg << std::endl;
      std::cerr << std::endl;
    }

    this->print_help(std::cerr);
  }

  if(external_exit)
    (*external_exit)(status);
  else {
    exit(status);
  }
}

const ArgumentParser::Argument & ArgumentParser::operator[](const std::string & name) const {
  Arguments::const_iterator it = success_parsed.find(name);

  if(it != success_parsed.end()) {
    return *(it->second);
  } else {
    throw std::range_error("No argument named \'" + name
			   + "\' was found in the parsed argument,"
			   + " make sur to specify it \'required\'"
			   + " or to give it a default value");
  }
}

/* -------------------------------------------------------------------------- */
void ArgumentParser::addArgument(const std::string & name_or_flag,
				 const std::string & help,
				 int nargs,
				 ArgumentType type) {
  _addArgument(name_or_flag, help, nargs, type);
}

ArgumentParser::_Argument & ArgumentParser::_addArgument(const std::string & name,
							 const std::string & help,
							 int nargs,
							 ArgumentType type) {
  _Argument * arg = NULL;

  switch(type) {
  case _string:  { arg = new ArgumentStorage<std::string>(); break; }
  case _float:   { arg = new ArgumentStorage<double>();      break; }
  case _integer: { arg = new ArgumentStorage<int>();         break; }
  case _boolean: { arg = new ArgumentStorage<bool>();        break; }
  }

  arg->help        = help;
  arg->nargs       = nargs;
  arg->type        = type;

  std::stringstream sstr(name);
  std::string item;
  std::vector<std::string> tmp_keys;
  while (std::getline(sstr, item, ';')) {
    tmp_keys.push_back(item);
  }

  int long_key  = -1;
  int short_key = -1;
  bool problem = (tmp_keys.size() > 2) || (name == "");
  for (std::vector<std::string>::iterator it = tmp_keys.begin(); it != tmp_keys.end(); ++it) {
    if(it->find("--") == 0) {
      problem |= (long_key != -1);
      long_key = it - tmp_keys.begin();
    } else if(it->find("-") == 0) {
      problem |= (long_key != -1);
      short_key = it - tmp_keys.begin();
    }
  }

  problem |= ((tmp_keys.size() == 2) && (long_key == -1 || short_key == -1));

  if(problem) {
    throw std::invalid_argument("Synthax of name or flags is not correct. Possible synthax are \'-f\', \'-f;--foo\', \'--foo\', \'bar\'");
  }

  if(long_key != -1) {
    arg->name = tmp_keys[long_key];
    arg->name.erase(0, 2);
  } else if(short_key != -1) {
    arg->name = tmp_keys[short_key];
    arg->name.erase(0, 1);
  } else {
    arg->name = tmp_keys[0];
    pos_args.push_back(arg);
    arg->required = true;
    arg->is_positional = true;
  }

  arguments[arg->name] = arg;

  if(!arg->is_positional) {
    if(short_key != -1) {
      std::string key = tmp_keys[short_key];
      key_args[key] = arg;
      arg->keys.push_back(key);
    }

    if(long_key != -1) {
      std::string key = tmp_keys[long_key];
      key_args[key] = arg;
      arg->keys.push_back(key);
    }
  }

  return *arg;
}


/* -------------------------------------------------------------------------- */
void ArgumentParser::parse(int& argc, char**& argv, int flags) {
  bool stop_in_not_parsed = flags & _stop_on_not_parsed;
  bool remove_parsed = flags & _remove_parsed;

  std::vector<std::string> argvs;

  for (int i = 0; i < argc; ++i) {
    argvs.push_back(std::string(argv[i]));
  }

  unsigned int current_position = 0;
  if(this->program_name == "") {
    std::string prog = argvs[current_position];
    const char * c_prog = prog.c_str();
    char * c_prog_tmp = strdup(c_prog);
    std::string base_prog(basename(c_prog_tmp));
    this->program_name = base_prog;
    std::free(c_prog_tmp);
  }

  std::queue<_Argument *> positional_queue;
  for(PositionalArgument::iterator it = pos_args.begin();
      it != pos_args.end(); ++it)
    positional_queue.push(*it);

  std::vector<int> argvs_to_remove;
  ++current_position; // consume argv[0]
  while(current_position < argvs.size()) {
    std::string arg = argvs[current_position];
    ++current_position;

    ArgumentKeyMap::iterator key_it = key_args.find(arg);

    bool is_positional = false;
    _Argument * argument_ptr = NULL;
    if(key_it == key_args.end()) {
      if(positional_queue.empty()) {
	if(stop_in_not_parsed) this->_exit("Argument " + arg + " not recognized", EXIT_FAILURE);
	continue;
      } else {
	argument_ptr = positional_queue.front();
	is_positional = true;
	--current_position;
      }
    } else {
      argument_ptr = key_it->second;
    }

    if(remove_parsed && ! is_positional) {
      argvs_to_remove.push_back(current_position - 1);
    }

    _Argument & argument = *argument_ptr;

    unsigned int min_nb_val = 0, max_nb_val = 0;
    switch(argument.nargs) {
    case _one_if_possible: max_nb_val = 1; break;                       // "?"
    case _at_least_one:    min_nb_val = 1;                              // "+"
    case _any:	           max_nb_val = argc - current_position; break; // "*"
    default: min_nb_val = max_nb_val = argument.nargs;                  // "N"
    }

    std::vector<std::string> values;
    unsigned int arg_consumed = 0;
    if(max_nb_val <= (argc - current_position)) {
      for(; arg_consumed < max_nb_val; ++arg_consumed) {
	std::string v = argvs[current_position];
	++current_position;
	bool is_key = key_args.find(v) != key_args.end();
	bool is_good_type = checkType(argument.type, v);

	if(!is_key && is_good_type) {
	  values.push_back(v);
	  if(remove_parsed) argvs_to_remove.push_back(current_position - 1);
	} else {
	  // unconsume not parsed argument for optional
	  if(!is_positional || (is_positional && is_key)) --current_position; 
	  break;
	}
      }
    }

    if(arg_consumed < min_nb_val) {
      if(!is_positional) {
	this->_exit("Not enought values for the argument "
		    + argument.name + " where provided", EXIT_FAILURE);
      } else {
	if(stop_in_not_parsed) this->_exit("Argument " + arg + " not recognized", EXIT_FAILURE);
      }
    } else {
      if (is_positional) positional_queue.pop();

      if(!argument.parsed) {
	success_parsed[argument.name] = &argument;
	argument.parsed = true;
	if((argument.nargs == _one_if_possible || argument.nargs == 0) && arg_consumed == 0) {
	  if(argument.has_const)
	    argument.setToConst();
	  else if(argument.has_default)
	    argument.setToDefault();
	} else {
	  argument.setValues(values);
	}
      } else {
	this->_exit("Argument " + argument.name
		    + " already present in the list of argument", EXIT_FAILURE);
      }
    }
  }

  for (_Arguments::iterator ait = arguments.begin();
       ait != arguments.end(); ++ait) {
    _Argument & argument = *(ait->second);
    if(!argument.parsed) {
      if(argument.has_default) {
	argument.setToDefault();
	success_parsed[argument.name] = &argument;
      }

      if(argument.required) {
	this->_exit("Argument " + argument.name + " required but not given!", EXIT_FAILURE);
      }
    }
  }


  // removing the parsed argument if remove_parsed is true
  if(argvs_to_remove.size()) {
    std::vector<int>::const_iterator next_to_remove = argvs_to_remove.begin();
    for (int i = 0, c = 0; i < argc; ++i) {
      if(next_to_remove == argvs_to_remove.end() || i != *next_to_remove) {
	argv[c] = argv[i];
	++c;
      } else {
	if (next_to_remove != argvs_to_remove.end()) ++next_to_remove;
      }
    }

    argc -= argvs_to_remove.size();
  }

  if(this->arguments["help"]->parsed) {
    this->_exit();
  }
}

/* -------------------------------------------------------------------------- */
bool ArgumentParser::checkType(ArgumentType type, const std::string & value) const {
  std::stringstream sstr(value);
  switch(type) {
  case _string:  { std::string s; sstr >> s; break; }
  case _float:   { double d;      sstr >> d; break; }
  case _integer: { int i;         sstr >> i; break; }
  case _boolean: { bool b;        sstr >> b; break; }
  }

  return (sstr.fail() == false);
}


/* -------------------------------------------------------------------------- */
void ArgumentParser::printself(std::ostream & stream) const {
  for(Arguments::const_iterator it = success_parsed.begin();
      it != success_parsed.end(); ++it) {
    const Argument & argument = *(it->second);
    argument.printself(stream);
    std::cout << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
void ArgumentParser::print_usage(std::ostream & stream) const {
  stream << "Usage: " << this->program_name;

  std::stringstream sstr_opt;
  // print shorten usage
  for(_Arguments::const_iterator it = arguments.begin();
      it != arguments.end(); ++it) {
    const _Argument & argument = *(it->second);
    if(!argument.is_positional) {
      if(!argument.required) stream << " [";
      stream << argument.keys[0];
      this->print_usage_nargs(stream, argument);
      if(!argument.required) stream << "]";
    }
  }

  for(PositionalArgument::const_iterator it = pos_args.begin();
      it != pos_args.end(); ++it) {
    const _Argument & argument = **it;
    this->print_usage_nargs(stream, argument);
  }

  stream << std::endl;
}

/* -------------------------------------------------------------------------- */
void ArgumentParser::print_usage_nargs(std::ostream & stream,
				      const _Argument & argument) const {
  std::string u_name = to_upper(argument.name);
  switch(argument.nargs) {
  case _one_if_possible: stream << " [" << u_name << "]"; break;
  case _at_least_one: stream << " " << u_name;
  case _any:  stream << " [" << u_name << " ...]"; break;
  default:
    for(int i = 0; i < argument.nargs; ++i) {
      stream << " " << u_name;
    }
  }
}

void ArgumentParser::print_help(std::ostream & stream) const {
  this->print_usage(stream);
  if(!pos_args.empty()) {
    stream << std::endl;
    stream << "positional arguments:" << std::endl;
    for(PositionalArgument::const_iterator it = pos_args.begin();
	it != pos_args.end(); ++it) {
      const _Argument & argument = **it;
      this->print_help_argument(stream, argument);
    }
  }

  if(!key_args.empty()) {
    stream << std::endl;
    stream << "optional arguments:" << std::endl;
    for(_Arguments::const_iterator it = arguments.begin();
	it != arguments.end(); ++it) {
      const _Argument & argument = *(it->second);
      if(!argument.is_positional) {
	this->print_help_argument(stream, argument);
      }
    }
  }
}

void ArgumentParser::print_help_argument(std::ostream & stream, const _Argument & argument) const {
  std::string key("");
  if(argument.is_positional)
    key =  argument.name;
  else {
    std::stringstream sstr;
    for(unsigned int i = 0; i < argument.keys.size(); ++i) {
      if(i != 0) sstr << ", ";
      sstr << argument.keys[i];
      this->print_usage_nargs(sstr, argument);
    }
    key = sstr.str();
  }

  stream << "  " << std::left << std::setw(15) << key << "  " << argument.help;
  argument.printDefault(stream);

  stream << std::endl;
}


}
