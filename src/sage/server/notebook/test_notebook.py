r"""nodoctest

THIS IS NOT USED ANYMORE -- should be deleted -- refers to very old version of notebook.

Test the notebook to analyze how it is behaving
after making some changes to it.

Given a Notebook session that had already been
run with the below command:

    sage: notebook(log_server=True)   # not tested

we then take the resulting 'server_log' and:
1) Initiate (spawn) a new Notebook session.
2) Pass all the inputs that the original session had.
3) We again record this *new* session.
4) Compare the original and new server logs.

"""

###########################################################################
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

#WARNING/NOTE: The Notebook testing functionality is a
#              *Work in Progress*. -Alex

#standard libraries
import httplib, urllib
import re, os, time

#SAGE libraries
import notebook
from sage.structure.sage_object import SageObject, load
from sage.interfaces.sage0 import Sage
from sage.plot.plot import bar_chart, text, graphics_array, Graphics

class Playback_Notebook_Session(SageObject):
    def __init__(self, dir="sage_notebook", sobj="nb.sobj"):
        """
        Test the Notebook by taking an already existing
        Notebook directory that was produced by running
        a Notebook session like so:

        "notebook(log_server=True)"

        After that session is completed a log will have
        been formed that has every GET and POST that the
        Notebook server handled.

        If the Notebook dir is not specified then
        dir='sage_notebook' and  sobj='nb.sobj'

        """
        self._NB = None
        self._NBTEST = None
        self._nb_dir = dir
        self._nb_sobj = sobj
        self._nb = "%s/%s" % (self._nb_dir, self._nb_sobj)

    def load_nb_sobj(self, dir):
        try:
            nb = load(dir)
        except AttributeError:
            print "Notebook did not load correctly."
            return
        return nb

    def get_NB_server_log(self):
        if self._NB is None:
            self._NB = self.load_nb_sobj(self._nb)
        return self._NB.server_log()

    def get_NBTEST_server_log(self):
        if self._NBTEST is None:
            print "You must run 'playback' first to get the NBTEST server_log."
            return
        return self._NBTEST.server_log()

    def move_old_NBTEST(self):
        """
        Safely 'remove' an old testing notebook...could be done differently.
        """
        if os.path.exists('NBTEST'):
            os.system('mv NBTEST OLDNBTEST')

    def playback(self, port=8000):
        """
        Play back the server log of a previously run
        notebook that had been run in the following way:

        EXAMPLES:
           sage: notebook(log_server=True)    # not tested
        """
        #get the test notebooks log
        server_log = self.get_NB_server_log()
        #remove (safely) any old testing notebooks
        self.move_old_NBTEST()
        #start a new SAGE process to run notebook
        S = Sage()
        print "Trying to open test notebook (%s) at port %s ..."%('NBTEST', port)
        ss = "N = notebook(dir='%s', port=%s, log_server=True)"%('NBTEST', port)
        print "sending this string: " + ss
        S._send(ss)
        time.sleep(0.1)
        #assuming S._get() opens the notebook, find at which port it was opened at:
        s = S._get()[1] #the first element is the response string
        port = re.compile(r'[0-9]{4}').findall(s)[0] #find port, which is a 4 digit number
        print "... the test notebook was opened at port %s"%port
        conn = httplib.HTTPConnection("localhost:%s"%port)
        time.sleep(0.1)
        #send requests in server_log to the new notebook, which is recorded
        for log in server_log:
            print log
            if log[0] == 'POST': #log[1] = path, log[2] = postvars, log[3] = 'content-type
                params = urllib.urlencode(log[2])
                headers = {"content-type":log[3], "accept":"text/plain"}
                conn.request('POST', log[1], params, headers)
                time.sleep(0.1)
                resp = conn.getresponse()
                if log[1] == '/cell_update': #give some time for the javascript??
                    time.sleep(0.1)
            else:
                conn.request('GET', log[1])
                resp = conn.getresponse()
                time.sleep(0.1)
            print S._get()[1]
        time.sleep(2)
        #S._send("N.save()")
        #g = S._get()[1]
        #print g
        conn.close()
        S._keyboard_interrupt()
        #g2 = S._get()[1]
        #print g2

    def compare_log_lengths(self):
        """
        Compare the lengths of the nblog and the testnblog.

        """
        nblog = self.get_NB_server_log()
        self._NBTEST = load("NBTEST/nb.sobj")
        nbtestlog = get_NBTEST_server_log()
        print "Original log has length: %s, Test log has length: %s"%(len(nblog), len(nbtestlog))
        uniqorig = len(set([l[1] for l in nblog]))
        uniqtest = len(set([l[1] for l in nbtestlog]))
        print "Orig log has %s unique requests, Test log has %s unique requests"%(uniqorig, uniqtest)

    def _gets_posts(self, log, gets=True, posts=True):
        """
        Given a server log, this function goes through
        the log and extracts either all the GET requests,
        all the POST requests, or both.

        """
        L = [] #in case both gets and posts are False
        if gets and posts:
            L = [r[1] for r in log]
        if gets and not posts:
            L = [r[1] for r in log if r[0] != "POST"]
        if posts and not gets:
            L = [r[1] for r in log if r[0] != "GET"]
        return L


    def _unique_occurances(self, log, gets=True, posts=True):
        """
        Counts unqiue occurances of GET, POST or both.

        """
        L = self._gets_posts(log, gets=gets, posts=posts)
        OD = {} #the occurances dictionary
        for req in L:
            OD[req] = OD.get(req, 0)+1
        return OD.items()


    def _track_occurances(self, log, gets=True, posts=True):
        """
        Given a server log, this function counts,
        in temporal order, all requests, (gets, posts, both)
        and returns a list of tuples, where each tuple
        has ("the_request", n), where n is the number of
        times "the_request" occured in a row.

        """
        L = self._gets_posts(log, gets=gets, posts=posts)
        OL = [] #the occurances list
        cnt = 1 #counter
        for req, n in zip(L, range(len(L))):
            try: #we will catch IndexError when we are at the end of OL.
                if req == L[n+1]: #if the next requests are the same
                    cnt += 1      #count them up
                else:
                    OL.append((req, cnt)) #the above req was different the last
                    cnt = 1               #so add tuple to OL and restart counter
            except IndexError:       #we are at the end of the OL list
                OL.append((req, cnt)) #so append the last req,cnt and be done.
        return OL

    def server_log_bar_chart(self, gets=True, posts=True):
        """
        Make a bar chart from a server log.

        """
        from matplotlib import colors
        nblog = self.get_NB_server_log()
        #self._NBTEST = load("NBTEST/nb.sobj")
        #self._NBTEST = load("test/nb.sobj")
        #nblog = self._NBTEST.server_log()
        #nbtestlog = get_NBTEST_server_log()
        #print nblog
        OL = self._track_occurances(nblog, gets=gets, posts=posts)
        #TL = self._track_occurances(nbtestlog, gets=gets, posts=posts)
        UL = list(set([t[0] for t in OL])) #unique request
        CL = colors.cnames.keys()[:len(UL)] #get colors for unique requests
        G1 = Graphics()
        G2 = Graphics()
        cnt = 1
        ll = len(OL)
        for req, col in zip(UL, CL): #loop over requests and colors
            f = lambda r,n: n if (r == req) else 0 #python 2.5!
            tl = [f(r,n) for r,n in OL]
            G1 += bar_chart(tl, rgbcolor=col)
            G2 += text(req, (3*ll/2, cnt), rgbcolor=col, fontsize=6)
            cnt += 1
        G = G1 + G1
        #G.axes_label(l=["time", "requests"])
        #return graphics_array([G1, G2])
        return G







notebook_playback = Playback_Notebook_Session

