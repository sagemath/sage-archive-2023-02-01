/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/commentator.C
 * Copyright (C) 1999 B. David Saunders,
 *                    Jean-Guillaume Dumas
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * This file implements the C++ interface to commentators (for
 * providing runtime commentary to the user)
 */

#include "linbox/linbox-config.h"

#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "linbox/util/commentator.h"
#include "linbox/util/debug.h"
#include "linbox/util/timer.h"

namespace LinBox
{
	// -----------------------------------------------------
	// Mathematical routines
	// -----------------------------------------------------
	double nroot (double a, long r, double precision)
	{
		long rm = r - 1 ;
		double c1 = double (rm) / double (r), c2 = a / double (r);
		double g = 1, pgr = 1, err = a - 1;

		while (err > precision) {
			g = g * c1 + c2 / pgr;
			pgr = pow (g, (double) rm);
			err = a - pgr * g;
			if (err < 0)
				err = -err;
		}

		return g;
	}

	long isnpower (long& l, long a)
	{
		long r = 2;
		double g;

		while ((g = nroot (a, r, 0.1)) >= 2) {
			l = (long) floor (g);
			if (g-double (l) > 0.1)
				++l;
			if (pow ((double) l, (double) r) == a)
				return r;
			++r;
		}

		return 0;
	}

	Commentator::Commentator ()
		: cnull (new nullstreambuf), _estimationMethod (BEST_ESTIMATE), _format (OUTPUT_CONSOLE),
		  _show_timing (true), _show_progress (true), _show_est_time (true)
	{
		//registerMessageClass (BRIEF_REPORT,         std::clog, 1, LEVEL_IMPORTANT);
		registerMessageClass (BRIEF_REPORT,         _report, 1, LEVEL_IMPORTANT);
		registerMessageClass (PROGRESS_REPORT,      _report);
		registerMessageClass (TIMING_MEASURE,       _report);
		registerMessageClass (TIMING_ESTIMATE,      _report);
		registerMessageClass (PARTIAL_RESULT,       _report);
		registerMessageClass (INTERNAL_WARNING,     _report, 10, LEVEL_NORMAL);
		registerMessageClass (INTERNAL_ERROR,       _report, 10, LEVEL_NORMAL);
		registerMessageClass (INTERNAL_DESCRIPTION, _report);
	}
	Commentator::Commentator (std::ostream& out)
		: cnull (new nullstreambuf), _estimationMethod (BEST_ESTIMATE), _format (OUTPUT_CONSOLE),
		  _show_timing (true), _show_progress (true), _show_est_time (true)
	{
		registerMessageClass (BRIEF_REPORT,         out, 1, LEVEL_IMPORTANT);
		registerMessageClass (BRIEF_REPORT,         out, 1, LEVEL_IMPORTANT);
		registerMessageClass (PROGRESS_REPORT,      _report);
		registerMessageClass (TIMING_MEASURE,       _report);
		registerMessageClass (TIMING_ESTIMATE,      _report);
		registerMessageClass (PARTIAL_RESULT,       _report);
		registerMessageClass (INTERNAL_WARNING,     _report, 10, LEVEL_NORMAL);
		registerMessageClass (INTERNAL_ERROR,       _report, 10, LEVEL_NORMAL);
		registerMessageClass (INTERNAL_DESCRIPTION, _report);
	}

	Commentator::~Commentator() {
		std::map <const char *, MessageClass *, LessThanString>::iterator i;
		for (i = _messageClasses.begin (); i != _messageClasses.end (); ++i)
			delete i->second;
	}

	void Commentator::start (const char *description, const char *fn, unsigned long len)
	{
		if (fn == (const char *) 0 && _activities.size () > 0)
			fn = _activities.top ()->_fn;

		if (isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, INTERNAL_DESCRIPTION, fn))
			report (LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) //<< "Starting activity: "
			<< description << std::endl;

		Activity *new_act = new Activity (description, fn, len);

		if (isPrinted (_activities.size (), LEVEL_IMPORTANT, BRIEF_REPORT, fn))
			printActivityReport (*new_act);

		_activities.push (new_act);

		new_act->_timer.start ();
	}

	void Commentator::startIteration (unsigned int iter, unsigned long len)
	{
		std::ostringstream str;

		str << "Iteration " << iter << std::ends;

		_iteration_str = str.str ();
		start (_iteration_str.c_str (), (const char *) 0, len);
	}

	void Commentator::stop (const char *msg, const char *long_msg, const char *fn)
	{
		float realtime, usertime, systime;
		Activity *top_act;

		linbox_check (_activities.top () != (Activity *) 0);
		linbox_check (msg != (const char *) 0);

		if (long_msg == (const char *) 0)
			long_msg = msg;

		top_act = _activities.top ();

		top_act->_timer.stop ();

		realtime = top_act->_timer.realtime ();
		usertime = top_act->_timer.usertime ();
		systime = top_act->_timer.systime ();

		if (realtime < 0) realtime = 0;
		if (usertime < 0) usertime = 0;
		if (systime < 0) systime = 0;

		if (fn != (const char *) 0 &&
		    _activities.size () > 0 &&
		    top_act->_fn != (const char *) 0 &&
		    strcmp (fn, top_act->_fn) != 0)
		{
			report (LEVEL_IMPORTANT, INTERNAL_WARNING)
				<< "Activity report mismatch. Check that start () and stop () calls are paired correctly." << std::endl;
		}

		fn = top_act->_fn;

		_activities.pop ();

		if (isPrinted (_activities.size (), LEVEL_IMPORTANT, BRIEF_REPORT, fn))
		{
		//std::cout << "calling finishA ";
			finishActivityReport (*top_act, msg);
			//std::cout << std::endl;
		}

		if (isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, INTERNAL_DESCRIPTION, fn)) {
			std::ostream &output = report (LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
			output.precision (4);
			output << "Finished activity (rea: " << realtime << "s, cpu: ";
			output.precision (4);
			output << usertime << "s, sys: ";
			output.precision (4);
			output << systime << "s): " << long_msg << std::endl;
		}
		else if (isPrinted (_activities.size (), LEVEL_IMPORTANT, INTERNAL_DESCRIPTION, fn)) {
			std::ostream &output = report (LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
			output.precision (4);
			output << "Completed activity: " << top_act->_desc << " (r: " << realtime << "s, u: ";
			output.precision (4);
			output << usertime << "s, s: ";
			output.precision (4);
			output << systime << "s) " << long_msg << std::endl;
		}

		delete top_act;
	}

	void Commentator::progress (long k, long len)
	{
		linbox_check (_activities.top () != (Activity *) 0);

		Activity *act = _activities.top ();
                Timer tmp = act->_timer;
		act->_timer.stop ();

		if (k == -1)
			act->_progress++;
		else
			act->_progress = k;

		if (len != -1)
			act->_len = len;

		if (act->_progress > act->_len)
			act->_len = act->_progress;

		std::ostream &rep = report (LEVEL_IMPORTANT, PROGRESS_REPORT);
		rep.precision (3);
		rep.setf (std::ios::fixed);
		rep << "Progress: " << act->_progress << " out of " << act->_len
		    << " (" << act->_timer.realtime () << "s elapsed)" << std::endl;

		if (_show_progress && isPrinted (_activities.size () - 1, LEVEL_IMPORTANT, BRIEF_REPORT, act->_fn))
			updateActivityReport (*act);
                act->_timer = tmp;
	}

	std::ostream &Commentator::report (long level, const char *msg_class)
	{
		linbox_check (msg_class != (const char *) 0);

		if (!isPrinted (_activities.size (), level, msg_class,
				(_activities.size () > 0) ? _activities.top ()->_fn : (const char *) 0))
			return cnull;

		MessageClass &messageClass = getMessageClass (msg_class);

		return messageClass._stream;
	}

	void Commentator::indent (std::ostream &stream) const
	{
		unsigned int i;

		for (i = 0; i < _activities.size (); i++)
			stream << "  ";
	}

	void Commentator::restoreActivityState (ActivityState state)
	{
		std::stack<Activity *> backup;

		while (!_activities.empty () && _activities.top () != state._act) {
			backup.push (_activities.top ());
			_activities.pop ();
		}

		if (_activities.empty ()) {
			// Uh oh -- the state didn't give a valid activity

			while (!backup.empty ()) {
				_activities.push (backup.top ());
				backup.pop ();
			}
		}
	}

	void Commentator::setMaxDepth (long depth)
	{
		MessageClass &briefReportClass = getMessageClass (BRIEF_REPORT);
		std::map <const char *, MessageClass *, LessThanString>::iterator i;

		for (i = _messageClasses.begin (); i != _messageClasses.end (); ++i)
			if (i->second != &briefReportClass)
				i->second->setMaxDepth (depth);
	}

	void Commentator::setMaxDetailLevel (long level)
	{
		MessageClass &briefReportClass = getMessageClass (BRIEF_REPORT);
		std::map <const char *, MessageClass *, LessThanString>::iterator i;

		for (i = _messageClasses.begin (); i != _messageClasses.end (); ++i)
			if (i->second != &briefReportClass)
				i->second->setMaxDetailLevel (level);
	}

	MessageClass &Commentator::registerMessageClass (const char *msg_class, std::ostream &stream, unsigned long max_depth, unsigned long max_level)
	{
		linbox_check (msg_class != (const char *) 0);

		MessageClass *new_obj = new MessageClass (*this, msg_class, stream, max_depth, max_level);
		_messageClasses[msg_class] = new_obj;
		return *new_obj;
	}

	MessageClass &Commentator::cloneMessageClass (const char *new_msg_class, const char *msg_class)
	{
		linbox_check (new_msg_class != (const char *) 0);
		linbox_check (msg_class != (const char *) 0);

		MessageClass *new_obj = new MessageClass (getMessageClass (msg_class));
		new_obj->_msg_class = new_msg_class;
		_messageClasses[new_msg_class] = new_obj;
		return *new_obj;
	}

	MessageClass &Commentator::cloneMessageClass (const char *new_msg_class, const char *msg_class, std::ostream &stream)
	{
		linbox_check (new_msg_class != (const char *) 0);
		linbox_check (msg_class != (const char *) 0);

		MessageClass &old_obj = getMessageClass (msg_class);
		MessageClass *new_obj = new MessageClass (*this, new_msg_class, stream, old_obj._configuration);
		_messageClasses[new_msg_class] = new_obj;
		return *new_obj;
	}

	MessageClass &Commentator::getMessageClass (const char *msg_class)
		{ return *_messageClasses[msg_class]; }

	void Commentator::setPrintParameters (unsigned long depth, unsigned long level, const char *fn)
	{
		MessageClass &briefReportClass = getMessageClass (BRIEF_REPORT);
		std::map <const char *, MessageClass *, LessThanString>::iterator i;

		for (i = _messageClasses.begin (); i != _messageClasses.end (); ++i)
			if (i->second != &briefReportClass)
				i->second->setPrintParameters (depth, level, fn);
	}

	void Commentator::setBriefReportParameters (OutputFormat format, bool show_timing, bool show_progress, bool show_est_time)
	{
		_format        = format;
		_show_timing   = show_timing;
		_show_progress = show_progress;
		_show_est_time = show_est_time;
	}

	bool Commentator::isPrinted (unsigned long depth, unsigned long level, const char *msg_class, const char *fn)
	{
		if (_messageClasses.find (msg_class) == _messageClasses.end ())
			return false;

		MessageClass &messageClass = getMessageClass (msg_class);

		return messageClass.isPrinted (depth, level, fn);
	}

	void Commentator::setBriefReportStream (std::ostream &stream)
		{ setMessageClassStream (BRIEF_REPORT, stream); }

	void Commentator::setReportStream (std::ostream &stream)
	{
		setMessageClassStream (PROGRESS_REPORT,      stream);
		setMessageClassStream (TIMING_MEASURE,       stream);
		setMessageClassStream (TIMING_ESTIMATE,      stream);
		setMessageClassStream (PARTIAL_RESULT,       stream);
		setMessageClassStream (INTERNAL_ERROR,       stream);
		setMessageClassStream (INTERNAL_WARNING,     stream);
		setMessageClassStream (INTERNAL_DESCRIPTION, stream);

		if (stream == getMessageClass (BRIEF_REPORT)._stream)
			getMessageClass (BRIEF_REPORT).setMaxDepth (0);
	}

	void Commentator::setMessageClassStream (const char *msg_class, std::ostream &stream)
	{
                //temporarily fixed the bug in test-commentator, left memory leaking.
                MessageClass *old_msg_class = _messageClasses[msg_class];
                cloneMessageClass (msg_class, msg_class, stream);
                delete old_msg_class;

	}

	void Commentator::setDefaultReportFile (const char *filename)
	{
		_report.open (filename);
	}

	void Commentator::printActivityReport (Activity &activity)
	{
		MessageClass &messageClass = getMessageClass (BRIEF_REPORT);

		if (_format == OUTPUT_CONSOLE) {
			messageClass._stream << activity._desc << "...";

			if (messageClass.isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, activity._fn))
				messageClass._stream << std::endl;
			else if (_show_progress && activity._len > 0) {
				messageClass._stream << "  0%";
				_last_line_len = strlen ("  0%");
			}
			else
				_last_line_len = 0;

			messageClass._smart_streambuf.stream ().flush ();
		}
		else if (_format == OUTPUT_PIPE &&
			 (((_show_progress || _show_est_time) && activity._len > 0) ||
			  messageClass.isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, activity._fn)))
		{
			messageClass._stream << activity._desc << "...";

			if (_show_progress)
				messageClass._stream << std::endl;
		}
	}

	void Commentator::updateActivityReport (Activity &activity)
	{
		MessageClass &messageClass = getMessageClass (BRIEF_REPORT);
		unsigned int i, old_len;
		std::ostringstream str;
		double percent = (double) activity._progress / (double) activity._len * 100.0;

		if (_format == OUTPUT_CONSOLE) {
			if (!messageClass.isPrinted (_activities.size (), LEVEL_IMPORTANT, activity._fn)) {
				if (_show_progress) {
					for (i = 0; i < _last_line_len; i++)
						messageClass._stream << '\b';
					str.width (3);
					str << floor (percent + 0.5) << '%' << std::ends;
					old_len = _last_line_len;
					_last_line_len = strlen (str.str ().c_str ());
					messageClass._stream << str.str ();
					for (int i = 0; i < (int) (old_len - _last_line_len); i++)
						messageClass._stream << ' ';
				}
			}
			else if (messageClass.isPrinted (_activities.size () - 1, LEVEL_UNIMPORTANT, activity._fn)) {
#if 0
				if (_show_est_time)
					messageClass._stream << activity._estimate.front ()._time
							     << " remaining" << std::endl;
#endif
			}

			messageClass._smart_streambuf.stream ().flush ();
		}
		else if (_format == OUTPUT_PIPE) {
			if (_show_progress) {
				messageClass._stream << floor (percent + 0.5) << "% done";
#if 0
				if (_show_est_time)
					messageClass._stream << " (" << activity._estimate.front ()._time
							     << " remaining)";
#endif
				messageClass._stream << std::endl;
			}
#if 0
			else if (_show_est_time)
				messageClass._stream << activity._estimate.front ()._time
						     << " remaining" << std::endl;
#endif
		}
	}

	void Commentator::finishActivityReport (Activity &activity, const char *msg)
	{
	//std::cout << "finishA " << _show_progress << _show_timing << std::endl;
		MessageClass &messageClass = getMessageClass (BRIEF_REPORT);
		unsigned int i;

		if (_format == OUTPUT_CONSOLE) {
			if (!messageClass.isPrinted (_activities.size () + 1, LEVEL_UNIMPORTANT, activity._fn)) {
				if (_show_progress)
					for (i = 0; i < _last_line_len; i++)
						messageClass._stream << '\b';

				messageClass._stream << msg;

				if (_show_timing)
					messageClass._stream << " (" << activity._timer.usertime () << " s)" << std::endl;
				else
					messageClass._stream << std::endl;
			}
			else if (messageClass.isPrinted (_activities.size (), LEVEL_UNIMPORTANT, activity._fn)) {
				for (i = 0; i < _activities.size (); i++)
					messageClass._stream << "  ";

				messageClass._stream << msg;
				//messageClass._stream << "Done: " << msg;

				if (_show_timing)
					messageClass._stream << " (" << activity._timer.usertime () << " s)" << std::endl;
				else
					messageClass._stream << std::endl;
			}

			messageClass._smart_streambuf.stream ().flush ();
		}
		else if (_format == OUTPUT_PIPE) {
			for (i = 0; i < _activities.size (); i++)
				messageClass._stream << "  ";

			if (((_show_progress || _show_est_time) && activity._len > 0) ||
			    messageClass.isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, activity._fn))
				messageClass._stream << "Done: " << msg << std::endl;
			else
				messageClass._stream << activity._desc << ": " << msg << std::endl;
		}
	}

	MessageClass::MessageClass (const Commentator &comm,
				    const char *msg_class,
				    std::ostream &stream,
				    unsigned long max_depth,
				    unsigned long max_level)
		: _msg_class (msg_class),
		  _smart_streambuf (comm, stream),
		  _stream (&_smart_streambuf),
		  _max_level (max_level),
		  _max_depth (max_depth)
	{
		fixDefaultConfig ();
	}

	void MessageClass::setMaxDepth (long depth)
	{
		_max_depth = (unsigned long) depth;
		fixDefaultConfig ();
	}

	void MessageClass::setMaxDetailLevel (long level)
	{
		_max_level = (unsigned long) level;
		fixDefaultConfig ();
	}

	void MessageClass::setPrintParameters (unsigned long depth, unsigned long level, const char *fn)
	{
		if (fn == (const char *) 0)
			fn = "";

		std::list <std::pair <unsigned long, unsigned long> > &config = _configuration[fn];
		std::list <std::pair <unsigned long, unsigned long> >::iterator i, j;

		i = config.begin ();

		// Iterate through preceeding elements in the std::list and remove
		// any that specify a lower level than we are using
		while (i != config.end () && i->first <= depth) {
			if (i->second <= level) {
				j = i++;
				config.erase (j);
			} else {
				++i;
			}
		}

		// Insert our new directive into the std::list
		if (i == config.end () || i->second != level)
			config.insert (i, std::pair <unsigned long, unsigned long> (depth, level));

		// Iterate through following elements in the std::list and remove any
		// that specify a higher level than we are using
		while (i != config.end ()) {
			if (i->second > level) {
				j = i++;
				config.erase (j);
			} else {
				++i;
			}
		}

		// End result: The std::list should be monotonically increasing in
		// the first parameter and decreasing in the second
	}

	bool MessageClass::isPrinted (unsigned long depth, unsigned long level, const char *fn)
	{
 		if (checkConfig (_configuration[""], depth, level))
		{	//std::cout << " fn=\"\", d " << depth << ", l " << level << " false" << std::endl;
 			return true;
		}
		else if (fn != (const char *) 0)
			//return checkConfig (_configuration[fn], depth, level);
			{ bool ans = checkConfig (_configuration[fn], depth, level);
			  if (ans)
		 		{	//std::cout << " fn=" << fn << ", d " << depth << ", l " << level << " true" << std::endl;
			 		return true;
		 		}
			  else
		 		{	//std::cout << " fn=" << fn << ", d " << depth << ", l " << level << " false" << std::endl;
			 		return false;
		 		}
		 	}

		else
		{	//std::cout << " fn=0, d " << depth << ", l " << level << " false" << std::endl;
			return false;
		}
	}

	MessageClass::MessageClass (const Commentator &comm,
				    const char *msg_class,
				    std::ostream &stream,
				    Configuration configuration)
		: _msg_class (msg_class),
		  _smart_streambuf (comm, stream),
		  _stream (&_smart_streambuf),
		  _configuration (configuration)
	{}

	void MessageClass::fixDefaultConfig ()
	{
		std::list <std::pair <unsigned long, unsigned long> > &config = _configuration[""];

		config.clear ();
		config.push_back (std::pair <unsigned long, unsigned long> (_max_depth, _max_level));
		config.push_back (std::pair <unsigned long, unsigned long> ((unsigned long) -1, Commentator::LEVEL_ALWAYS));
	}

	bool MessageClass::checkConfig (std::list <std::pair <unsigned long, unsigned long> > &config, unsigned long depth, unsigned long level)
	{
		std::list <std::pair <unsigned long, unsigned long> >::iterator i;

		i = config.begin ();
		while (i != config.end ()) {
			if (depth < i->first) {
				if (level <= i->second)
					return true;
				else
					return false;
			}

			i++;
		}

		return false;
	}

	void MessageClass::dumpConfig () const
	{
		Configuration::const_iterator i;
		std::list <std::pair <unsigned long, unsigned long> >::const_iterator j;

		for (i = _configuration.begin (); i != _configuration.end (); i++) {
			std::cerr << "Configuration (" << (*i).first << "):" << std::endl;

			for (j = (*i).second.begin (); j != (*i).second.end (); j++)
				std::cerr << "  Depth: " << (*j).first << ", Level: " << (*j).second << std::endl;

			std::cerr << std::endl;
		}
	}

	int MessageClass::smartStreambuf::sync ()
	{
		std::streamsize n = pptr () - pbase ();
		return (n && writeData (pbase (), n) != n) ? EOF : 0;
	}

	int MessageClass::smartStreambuf::overflow (int ch)
	{
		std::streamsize n = pptr () - pbase ();

		if (n && sync ())
			return EOF;

		if (ch != EOF) {
			char cbuf[1];
			cbuf[0] = ch;
			if (writeData (cbuf, 1) != 1)
				return EOF;
		}

		pbump (-n);
		return 0;
	}

	std::streamsize MessageClass::smartStreambuf::xsputn (const char *text, std::streamsize n)
	{
		return (sync () == EOF) ? 0 : writeData (text, n);
	}

	int MessageClass::smartStreambuf::writeData (const char *text, std::streamsize n)
	{
		std::streamsize idx;
		std::streamsize m = n;

		if (_indent_next) {
			_comm.indent (_stream);
			_indent_next = false;
		}

		for (idx = 0; (idx < m) &&(text[idx] != '\n') ; ++idx);

		while (idx < m) {
			_stream.write (text, idx + 1);
			m -= idx + 1;

			if (m > 0)
				_comm.indent (_stream);
			else
				_indent_next = true;

			text += idx + 1;
			for (idx = 0; idx != '\n' && idx < m; ++idx);
		}

		_stream.write (text, m);

		_stream.flush ();

		return n;
	}

	// Default global commentator
	//Commentator commentator;
}
