import wave
from sage.plot.plot import list_plot
from sage.structure.sage_object import SageObject
from sage.misc.html import html

class Wave(SageObject):
    """
    A class wrapping a wave audio file.

    INPUT:
        You must call Wave() with either data = filename, where
        filename is the name of a wave file, or with each of the
        following options:

            channels  -- the number of channels in the wave file (1 for
                        mono, 2 for stereo, etc...
            width     -- the number of bytes per sample
            framerate -- the number of frames per second
            nframes   -- the number of frames in the data stream
            bytes     -- a string object containing the bytes of the
                         data stream

    Slicing:
        Slicing a Wave object returns a new wave object that has been
        trimmed to the bytes that you have given it.

    Indexing:
        Getting the $n$th item in a Wave object will give you the value
        of the $n$th frame.
    """
    def __init__(self, data=None, **kwds):
        if data is not None:
            self._filename = data
            wv = wave.open(data, "rb")
            self._nchannels = wv.getnchannels()
            self._width = wv.getsampwidth()
            self._framerate = wv.getframerate()
            self._nframes = wv.getnframes()
            self._bytes = wv.readframes(self._nframes)
            self._channel_data = self.separate_channels()
            wv.close()
        elif kwds:
            try:
                self._nchannels = kwds['nchannels']
                self._width = kwds['width']
                self._framerate = kwds['framerate']
                self._nframes = kwds['nframes']
                self._bytes = kwds['bytes']
                self._channel_data = kwds['channel_data']
            except KeyError:
                raise KeyError, "invalid input to Wave initializer"
        else:
            raise ValueError, "Must give a filename"

    def channel_data(self, n):
        return self._channel_data[n]

    def separate_channels(self):
        data = self._bytes
        l = len(data) / (self._width)
        channel_data = [[] for i in xrange(self._nchannels)]
        if self._width == 1:
            # handle the one byte case
            for n in xrange(l):
                channel_data[n % self._nchannels].append(ord(data[n])-127)

        elif self._width == 2:
            for n in xrange(l):
                # compute the value as an integer
                x = ord(data[2*n]) + 256 * ord(data[2*n + 1])
                if x > 32768:
                    x -= 65536
                channel_data[n % self._nchannels].append(x)
        else:
            raise NotImplementedError, "greater than 16-bit wavs not supported"

        return channel_data

    def getnchannels(self):
        return self._nchannels

    def getsampwidth(self):
        return self._width

    def getframerate(self):
        return self._framerate

    def getnframes(self):
        return self._nframes

    def readframes(self, nframes):
        return self._bytes[:nframes*self._width]

    def getlength(self):
        return float(self._nframes) / (self._nchannels * float(self._framerate))

    def _repr_(self):
        nc = self.getnchannels()
        return "Wave file with %s channel%s of length %s seconds%s" % \
        (nc, "" if nc == 1 else "s", self.getlength(), "" if nc == 1 else " each")

    def plot(self, channel=0, npoints=None, plotjoined=True, **kwds):
        """
        Plots the audio data.
        """
        if npoints == None:
            npoints = self._nframes
        # figure out on what intervals to sample the data
        seconds = float(self._nframes) / float(self._width)
        sample_step = seconds / float(npoints)

        domain = [float(n * sample_step) / float(self._framerate) for n in xrange(npoints)]
        # now, how many of the frames do we sample?
        frame_skip = self._nframes / npoints
        # the values of the function at each point in the domain

        values = [self.channel_data(channel)[frame_skip*i] for i in xrange(npoints)]
        # now scale the values
        scale = 1 << (8*self._width -1)
        values = [float(s) / float(scale) for s in values]
        points = zip(domain, values)

        return list_plot(points, plotjoined=plotjoined, **kwds)

    def plot_raw(self, channel=0, npoints=None, plotjoined=True, **kwds):
        if npoints == None:
            npoints = self._nframes
        seconds = float(self._nframes) / float(self._width)
        sample_step = seconds / float(npoints)
        domain = [float(n*sample_step) / float(self._framerate) for n in xrange(npoints)]
        frame_skip = self._nframes / npoints
        values = [self.channel_data(channel)[frame_skip*i] for i in xrange(npoints)]
        points = zip(domain, values)

        return list_plot(points, plotjoined=plotjoined, **kwds)

    # returns the ith frame of data in the wave, in the form of a string
    def __getitem__(self, i):
        n = i*self._width
        return self._bytes[n:n+self._width]

    def slice_seconds(self, start, stop):
        """
        Slices the wave from start to stop.

        INPUT:
            start -- the time index from which to begin the slice
            stop  -- the time index from which to end the slice

        OUTPUT:
            A Wave object whose data is this objects's data,
            sliced between the given time idices
        """
        start = int(start*self.getframerate())
        stop = int(stop*self.getframerate())
        return self[start:stop]

    def __getslice__(self, start, stop):
        return self.__copy__(start, stop)

    # start and stop are frame numbers
    def __copy__(self, start, stop):
        start = start * self._width
        stop = stop * self._width
        channels_sliced = [self._channel_data[i][start:stop] for i in range(self._nchannels)]

        return Wave(nchannels = self._nchannels,
                    width = self._width,
                    framerate = self._framerate,
                    bytes = self._bytes[start:stop],
                    nframes = stop - start,
                    channel_data = channels_sliced)
