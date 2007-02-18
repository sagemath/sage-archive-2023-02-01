import os

from sage.dsage.database.job import Job
from sage.dsage.dist_functions.dist_function import DistributedFunction

class DistributedPOVRay(DistributedFunction):
    r"""
    DistributedPOVRay distributes rendering of a .pov file.

    Parameters:
        DSage -- a DSage object
        name -- the name of your job
        files -- a list of files required for your pov job, including
                the pov file itself
        splits -- number of jobs you want to split into
        **kwargs -- parameters you wish to pass to povray
            Note: You must pass in at least a width and height.

    OUTPUT:
        A number of .ppm files (depending on split) and a final .ppm file
        which is the combination of all the rendered parts.

    AUTHOR:
        Yi Qiang

    """

    def __init__(self, DSage, name, files, splits, **kwargs):
        DistributedFunction.__init__(self, DSage)
        self.name = name
        self.files = files
        for f in self.files:
            if f.endswith('.pov'):
                self.pov_fname = f

        self.splits = splits
        self.width = kwargs['W']
        self.height = kwargs['H']
        self.kwargs = kwargs

        # Figure out how to split up the image into jobs
        self.remainder = self.height % self.splits
        self.step = self.height // self.splits
        self.sr = 1
        self.er = self.sr + self.step
        self.n = 0
        self.n1 = 0
        for i in xrange(splits):
            job = self.next_job()
            if not job == None:
                self.outstanding_jobs.append(job)

    def next_job(self):
        if self.n + 1 == self.splits:
            self.er = self.height
        if self.n + 1 > self.splits:
            return

        job_file = "povray('%s', outfile='%s_%04d.ppm', " % (self.pov_fname,
                                                             self.name,
                                                             self.n)

        for k, v in self.kwargs.iteritems():
            job_file += '%s=%s, ' % (k, v)

        job_file += 'SR=%s, ER=%s, ' % (self.sr, self.er)
        job_file += ')'
        job_file += '\n'
        job_file += "tmp = open('%s_%04d.ppm', 'rb').read()\n" % (self.name,
                                                                  self.n)
        job_file += "save(tmp, '%s_%04d')\n" % (self.name, self.n)
        job_file += "DSAGE_RESULT = '%s_%04d' + '.sobj'\n" % (self.name,
                                                              self.n)
        job_file += '\n'

        self.job_files.append(job_file)

        job = Job(file=job_file, name='%s_%04d.ppm' % (self.name, self.n))

        for file in self.files:
            job.attach_file(file)
        self.n += 1

        self.sr = self.er + 1
        self.er += self.step

        return job

    def process_result(self, job):
        if self.done:
            return
        f = open(job.name, 'wb')
        f.write(job.result)
        f.close()
        self.n1 += 1
        if len(self.waiting_jobs) == 0:
            print "Got all the images, stiching them now.\n"
            cmd = "combineppm %s > %s.ppm" % (self.name, self.name)
            os.system(cmd)
            self.done = True
