/** @file excompiler.cpp
 *
 *  Functions to facilitate the conversion of a ex to a function pointer suited for
 *  fast numerical integration.
 *
 */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "excompiler.h"

#include <stdexcept>
#include <ios>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBDL
#include <dlfcn.h>
#endif // def HAVE_LIBDL

#include "ex.h"
#include "lst.h"
#include "operators.h"
#include "relational.h"
#include "symbol.h"

namespace GiNaC {

#ifdef HAVE_LIBDL
	
/**
 * Small class that manages modules opened by libdl. It is used by compile_ex
 * and link_ex in order to have a clean-up of opened modules and their
 * associated source and so-files at the time of program termination. It is
 * supposed to be statically instantiated once (see definition of
 * global_excompiler below). The dtor of that object then performs the
 * clean-up. On top of that it provides some functionality shared between
 * different compile_ex and link_ex specializations.
 */
class excompiler
{
	/**
	 * Holds all necessary information about opened modules.
	 */
	struct filedesc
	{
		void* module;
		std::string name; /**< filename with .so suffix */
		bool clean_up; /**< if true, source and so-file will be deleted */
	};
	std::vector<filedesc> filelist; /**< List of all opened modules */
public:
	/**
	 * Complete clean-up of opend modules is done on destruction.
	 */
	~excompiler()
	{
		for (std::vector<filedesc>::const_iterator it = filelist.begin(); it != filelist.end(); ++it) {
			clean_up(it);
		}
	}
	/**
	 * Adds a new module to the list.
	 */
	void add_opened_module(void* module, const std::string& name, bool clean_up)
	{
		filedesc fd;
		fd.module = module;
		fd.name = name;
		fd.clean_up = clean_up;
		filelist.push_back(fd);
	}
	/**
	 * Closes a module and deletes the so-file if the associated clean-up flag is true.
	 */
	void clean_up(const std::vector<filedesc>::const_iterator it)
	{
		dlclose(it->module);
		if (it->clean_up) {
			remove(it->name.c_str());
		}
	}
	/**
	 * Creates a new C source file and adds a standard header. If filename is
	 * empty, a unique random name is produced and used.
	 */
	void create_src_file(std::string& filename, std::ofstream& ofs)
	{
		if (filename.empty()) {
			// fill filename with unique random word
			const char* filename_pattern = "./GiNaCXXXXXX";
			char* new_filename = new char[strlen(filename_pattern)+1];
			strcpy(new_filename, filename_pattern);
			if (!mktemp(new_filename)) {
				delete new_filename;
				throw std::runtime_error("mktemp failed");
			}
			filename = std::string(new_filename);
			ofs.open(new_filename, std::ios::out);
			delete new_filename;
		} else {
			// use parameter as filename
			ofs.open(filename.c_str(), std::ios::out);
		}
		
		if (!ofs) {
			throw std::runtime_error("could not create source code file for compilation");
		}

		ofs << "#include <stddef.h> " << std::endl;
		ofs << "#include <stdlib.h> " << std::endl;
		ofs << "#include <math.h> " << std::endl;
		ofs << std::endl;
	}
	/**
	 * Calls the shell script 'ginac-excompiler' to compile the produced C
	 * source file into an linkable so-file.  On demand the C source file is
	 * deleted.
	 */
	void compile_src_file(const std::string filename, bool clean_up)
	{
		std::string strcompile = "ginac-excompiler " + filename;
		if (system(strcompile.c_str())) {
			throw std::runtime_error("excompiler::compile_src_file: error compiling source file!");
		}
		if (clean_up) {
			remove(filename.c_str());
		}
	}
	/**
	 * Links a so-file whose filename is given.
	 */
	void* link_so_file(const std::string filename, bool clean_up)
	{
		void* module = NULL;
		module = dlopen(filename.c_str(), RTLD_NOW);
		if (module == NULL)	{
			throw std::runtime_error("excompiler::link_so_file: could not open compiled module!");
		}

		add_opened_module(module, filename, clean_up);

		return dlsym(module, "compiled_ex");
	}
	/**
	 * Removes a modules from the module list. Performs a clean-up before that.
	 * Every module with the given name will be affected.
	 */
	void unlink(const std::string filename)
	{
		for (std::vector<filedesc>::iterator it = filelist.begin(); it != filelist.end();) {
			if (it->name == filename) {
				clean_up(it);
				filelist.erase(it);
			} else {
				++it;
			}
		}
	}
};

/**
 * This static object manages the modules opened by the complile_ex and link_ex
 * functions. On program termination its dtor is called and all open modules
 * are closed. The associated source and so-files are eventually deleted then
 * as well.
 * In principle this could lead to a static deconstruction order fiasco, if
 * other code from this library uses the compile_ex and link_ex functions
 * (which it doesn't at the moment and won't in the likely future, so therefore
 * we ignore this issue).
 */
static excompiler global_excompiler;

void compile_ex(const ex& expr, const symbol& sym, FUNCP_1P& fp, const std::string filename)
{
	symbol x("x");
	ex expr_with_x = expr.subs(lst(sym==x));

	std::ofstream ofs;
	std::string unique_filename = filename;
	global_excompiler.create_src_file(unique_filename, ofs);

	ofs << "double compiled_ex(double x)" << std::endl;
	ofs << "{" << std::endl;
	ofs << "double res = ";
	expr_with_x.print(GiNaC::print_csrc_double(ofs));
	ofs << ";" << std::endl;
	ofs << "return(res); " << std::endl;
	ofs << "}" << std::endl;

	ofs.close();

	global_excompiler.compile_src_file(unique_filename, filename.empty());
	// This is not standard compliant! ... no conversion between
	// pointer-to-functions and pointer-to-objects ...
	fp = (FUNCP_1P) global_excompiler.link_so_file(unique_filename+".so", filename.empty());
}

void compile_ex(const ex& expr, const symbol& sym1, const symbol& sym2, FUNCP_2P& fp, const std::string filename)
{
	symbol x("x"), y("y");
	ex expr_with_xy = expr.subs(lst(sym1==x, sym2==y));

	std::ofstream ofs;
	std::string unique_filename = filename;
	global_excompiler.create_src_file(unique_filename, ofs);

	ofs << "double compiled_ex(double x, double y)" << std::endl;
	ofs << "{" << std::endl;
	ofs << "double res = ";
	expr_with_xy.print(GiNaC::print_csrc_double(ofs));
	ofs << ";" << std::endl;
	ofs << "return(res); " << std::endl;
	ofs << "}" << std::endl;

	ofs.close();

	global_excompiler.compile_src_file(unique_filename, filename.empty());
	// This is not standard compliant! ... no conversion between
	// pointer-to-functions and pointer-to-objects ...
	fp = (FUNCP_2P) global_excompiler.link_so_file(unique_filename+".so", filename.empty());
}

void compile_ex(const lst& exprs, const lst& syms, FUNCP_CUBA& fp, const std::string filename)
{
	lst replacements;
	for (int count=0; count<syms.nops(); ++count) {
		std::ostringstream s;
		s << "a[" << count << "]";
		replacements.append(syms.op(count) == symbol(s.str()));
	}

	std::vector<ex> expr_with_cname;
	for (int count=0; count<exprs.nops(); ++count) {
		expr_with_cname.push_back(exprs.op(count).subs(replacements));
	}

	std::ofstream ofs;
	std::string unique_filename = filename;
	global_excompiler.create_src_file(unique_filename, ofs);

	ofs << "void compiled_ex(const int* an, const double a[], const int* fn, double f[])" << std::endl;
	ofs << "{" << std::endl;
	for (int count=0; count<exprs.nops(); ++count) {
		ofs << "f[" << count << "] = ";
		expr_with_cname[count].print(GiNaC::print_csrc_double(ofs));
		ofs << ";" << std::endl;
	}
	ofs << "}" << std::endl;

	ofs.close();

	global_excompiler.compile_src_file(unique_filename, filename.empty());
	// This is not standard compliant! ... no conversion between
	// pointer-to-functions and pointer-to-objects ...
	fp = (FUNCP_CUBA) global_excompiler.link_so_file(unique_filename+".so", filename.empty());
}

void link_ex(const std::string filename, FUNCP_1P& fp)
{
	// This is not standard compliant! ... no conversion between
	// pointer-to-functions and pointer-to-objects ...
	fp = (FUNCP_1P) global_excompiler.link_so_file(filename, false);
}

void link_ex(const std::string filename, FUNCP_2P& fp)
{
	// This is not standard compliant! ... no conversion between
	// pointer-to-functions and pointer-to-objects ...
	fp = (FUNCP_2P) global_excompiler.link_so_file(filename, false);
}

void link_ex(const std::string filename, FUNCP_CUBA& fp)
{
	// This is not standard compliant! ... no conversion between
	// pointer-to-functions and pointer-to-objects ...
	fp = (FUNCP_CUBA) global_excompiler.link_so_file(filename, false);
}

void unlink_ex(const std::string filename)
{
	global_excompiler.unlink(filename);
}

#else // def HAVE_LIBDL

/*
 * In case no working libdl has been found by configure, the following function
 * stubs preserve the interface. Every function just raises an exception.
 */

void compile_ex(const ex& expr, const symbol& sym, FUNCP_1P& fp, const std::string filename)
{
	throw std::runtime_error("compile_ex has been disabled because of missing libdl!");
}

void compile_ex(const ex& expr, const symbol& sym1, const symbol& sym2, FUNCP_2P& fp, const std::string filename)
{
	throw std::runtime_error("compile_ex has been disabled because of missing libdl!");
}

void compile_ex(const lst& exprs, const lst& syms, FUNCP_CUBA& fp, const std::string filename)
{
	throw std::runtime_error("compile_ex has been disabled because of missing libdl!");
}

void link_ex(const std::string filename, FUNCP_1P& fp)
{
	throw std::runtime_error("link_ex has been disabled because of missing libdl!");
}

void link_ex(const std::string filename, FUNCP_2P& fp)
{
	throw std::runtime_error("link_ex has been disabled because of missing libdl!");
}

void link_ex(const std::string filename, FUNCP_CUBA& fp)
{
	throw std::runtime_error("link_ex has been disabled because of missing libdl!");
}

void unlink_ex(const std::string filename)
{
	throw std::runtime_error("unlink_ex has been disabled because of missing libdl!");
}

#endif // def HAVE_LIBDL

} // namespace GiNaC
