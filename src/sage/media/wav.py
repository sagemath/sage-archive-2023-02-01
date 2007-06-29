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
            width     -- the number of bytes per frame
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
            self._channels = wv.getnchannels()
            self._width = wv.getsampwidth()
            self._framerate = wv.getframerate()
            self._nframes = wv.getnframes()
            self._bytes = wv.readframes(self._nframes)
            wv.close()
        elif kwds:
            self._channels = kwds['channels']
            self._width = kwds['width']
            self._framerate = kwds['framerate']
            self._nframes = kwds['nframes']
            self._bytes = kwds['bytes']
        else:
            raise ValueError, "Must give a filename"

    def getchannels(self):
        return self._channels

    def getsampwidth(self):
        return self._width

    def getframerate(self):
        return self._framerate

    def getnframes(self):
        return self._nframes

    def readframes(self, nframes):
        return self._bytes[:nframes]

    def getlength(self):
        return float(self._nframes) / float(self._framerate)

    def _repr_(self):
        return "Wave file of length %s seconds"% self.getlength()

    def plot(self, npoints=None, **kwds):
        """
        Plots the audio data.
        """
        if npoints == None:
            npoints = self._nframes
        # figure out on what intervals to sample the data
        seconds = float(self._nframes) / float(self._width)
        sample_step = seconds / float(npoints)

        domain = [float(n * sample_step) / float(self._framerate) for n in range(npoints)]
        # now, how many of the frames do we sample?
        frame_skip = self._nframes / npoints
        # the values of the function at each point in the domain
        #values = [self._bb for b in [frame_skip*i for i in range(npoints)]]
        offset = 2**(self._width*8 -1)
        values = [ord(self._bytes[i]) - offset for i in [frame_skip*n for n in range(npoints)]]
        # now scale the values
        values = [float(s) / float(offset) for s in values]
        points = zip(domain, values)

        return list_plot(points, **kwds)

    def __getitem__(self, i):
        return ord(self._bytes[i])

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

    def __copy__(self, start, stop):
        return Wave(channels = self._channels,
                    width = self._width,
                    framerate = self._framerate,
                    bytes = self._bytes[start:stop],
                    nframes = stop - start)
