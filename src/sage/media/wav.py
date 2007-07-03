r"""
Work with WAV files.

A WAV file is a header specifying format information, followed by a
sequence of bytes, representing the state of some audio signal over a
lenght of time.

A WAV file may have any number of channels. Typically, they have 1
(mono) or 2 (for stereo). The data of a WAV file is given as a
sequence of frames. A frame consists of samples. There is one sample
per channel, per frame. Every wav file has a sample width, or, the
number of bits per sample. Typically these are either 8 or 16 bits.

The wav module supplies more convenient access to this data. In
particular, see the docstring for \code{Wave.channel_data()}.

The header contains information necessary for playing the WAV file,
including the number of frames per second, the number of bytes per
sample, and the number of channels in the file.

AUTHORS:
    -- Bobby Moretti and Gonzolo Tornaria (2007-07-01): First version
    -- William Stein (2007-07-03): add more
    -- Bobby Moretti (2007-07-03): add doctests
"""

import math
import os
import wave

from sage.plot.plot import list_plot
from sage.structure.sage_object import SageObject
from sage.misc.all import srange
from sage.misc.html import html
from sage.rings.all import RDF

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
            self._name = os.path.split(data)[1]
            wv = wave.open(data, "rb")
            self._nchannels = wv.getnchannels()
            self._width = wv.getsampwidth()
            self._framerate = wv.getframerate()
            self._nframes = wv.getnframes()
            self._bytes = wv.readframes(self._nframes)
            self._channel_data = self._separate_channels()
            wv.close()
        elif kwds:
            try:
                self._nchannels = kwds['nchannels']
                self._width = kwds['width']
                self._framerate = kwds['framerate']
                self._nframes = kwds['nframes']
                self._bytes = kwds['bytes']
                self._channel_data = kwds['channel_data']
                self._name = kwds['name']
            except KeyError:
                raise KeyError, "invalid input to Wave initializer"
        else:
            raise ValueError, "Must give a filename"

    def save(self, filename='sage.wav'):
        r"""
        Save this wave file to disk, either as a SAGE sobj or as a .wav file.

        INPUT:
            filename -- the path of the file to save. If filename ends
                        with 'wav', then save as a wave file,
                        otherwise, save a SAGE object.

            If no input is given, save the file as 'sage.wav'.

        """
        if not filename.endswith('.wav'):
            SageObject.save(self, filename)
            return
        wv = wave.open(filename, 'wb')
        wv.setnchannels(self._nchannels)
        wv.setsampwidth(self._width)
        wv.setframerate(self._framerate)
        wv.setnframes(self._nframes)
        wv.writeframes(self._bytes)
        wv.close()

    def listen(self):
        """
        Listen to (or download) this wave file.

        Creates a link to this wave file in the notebook.
        """
        from sage.misc.html import html
        i = 0
        fname = 'sage%s.wav'%i
        while os.path.exists(fname):
            i += 1
            fname = 'sage%s.wav'%i

        self.save(fname)
        return html('<a href="cell://%s">Click to listen to %s</a>'%(fname, self._name))

    def channel_data(self, n):
        """
        Get the data from a given channel.

        INPUT:
            n -- the channel number to get

        OUTPUT:
            A list of signed ints, each containing the value of a frame.
        """
        return self._channel_data[n]

    def _separate_channels(self):
        """
        Separates the channels. This is an internal helper method for
        precomputing some data while initializing the wav object.
        """
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
        """
        Returns the number of channels in this wave object.

        OUTPUT:
            The number of channels in this wave file.
        """
        return self._nchannels

    def getsampwidth(self):
        """
        Returns the number of bytes per sample in this wave object.

        OUTPUT:
            The number of bytes in each sample.
        """
        return self._width

    def getframerate(self):
        """
        Returns the number of frames per second in this wave object.

        OUTPUT:
            The frame rate of this sound file.
        """
        return self._framerate

    def getnframes(self):
        """
        The total number of frames in this wave object.

        OUTPUT:
            The number of frames in this WAV.
        """
        return self._nframes

    def readframes(self, n):
        """
        Reads out the raw data for the first $n$ frames of this wave
        object.

        INPUT:
            n -- the number of frames to return

        OUTPUT:
            A list of bytes (in string form) representing the raw wav data.
        """
        return self._bytes[:nframes*self._width]

    def getlength(self):
        """
        Returns the length of this file (in seconds).

        OUTPUT:
            The running time of the entire WAV object.
        """
        return float(self._nframes) / (self._nchannels * float(self._framerate))

    def _repr_(self):
        nc = self.getnchannels()
        return "Wave file %s with %s channel%s of length %s seconds%s" % \
        (self._name, nc, "" if nc == 1 else "s", self.getlength(), "" if nc == 1 else " each")

    def _normalize_npoints(self, npoints):
        """
        Used internally while plotting to normalize the number of
        """
        return npoints if npoints else self._nframes

    def domain(self, npoints=None):
        """
        Used internally for plotting. Get the x-values for the various points to plot.
        """
        npoints = self._normalize_npoints(npoints)
        # figure out on what intervals to sample the data
        seconds = float(self._nframes) / float(self._width)
        sample_step = seconds / float(npoints)

        domain = [float(n * sample_step) / float(self._framerate) for n in xrange(npoints)]
        return domain

    def values(self, npoints=None, channel=0):
        """
        Used internally for plotting. Get the y-values for the various points to plot.
        """
        npoints = self._normalize_npoints(npoints)

        # now, how many of the frames do we sample?
        frame_skip = int(self._nframes / npoints)
        # the values of the function at each point in the domain

        values = [self.channel_data(channel)[frame_skip*i] for i in xrange(npoints)]
        # now scale the values
        scale = 1 << (8*self._width -1)
        values = [float(s) / float(scale) for s in values]
        return values

    def vector(self, npoints=None, channel=0):
        npoints = self._normalize_npoints(npoints)

        V = RDF**npoints
        return V(self.values(npoints=npoints, channel=channel))

    def plot(self, npoints=None, channel=0, plotjoined=True, **kwds):
        """
        Plots the audio data.

        INPUT:
            npoints -- number of sample points to take; if not given, draws
                       all known points.
            channel -- 0 or 1 (if stereo).  default: 0
            plotjoined -- whether to just draw dots or draw lines between sample points

        OUTPUT:
            a plot object that can be shown.
        """

        domain = self.domain(npoints = npoints)
        values = self.values(npoints=npoints, channel = channel)
        points = zip(domain, values)

        L = list_plot(points, plotjoined=plotjoined, **kwds)
        L.xmin(0)
        L.xmax(domain[-1])
        return L

    def plot_fft(self, npoints=None, channel=0, **kwds):
        v = self.vector(npoints=npoints)
        w = v.fft()
        z = [abs(x) for x in w]
        twopi = 2*math.pi
        data = zip(srange(0, twopi, twopi/len(z)),  z)
        L = list_plot(data, plotjoined=True, **kwds)
        L.xmin(0)
        L.xmax(twopi)
        return L

    def plot_raw(self, npoints=None, channel=0, plotjoined=True, **kwds):
        npoints = self._normalize_npoints(npoints)
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
            start -- the time index from which to begin the slice (in seconds)
            stop -- the time index from which to end the slice (in seconds)

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
                    channel_data = channels_sliced,
                    name = self._name)
