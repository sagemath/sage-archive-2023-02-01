# -*- encoding: utf-8 -*-
r"""
Video Output Types

This module defines the rich output types for video formats.
"""

#*****************************************************************************
#       Copyright (C) 2015 Martin von Gagern <Martin.vGagern@gmx.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os

from sage.repl.rich_output.output_basic import OutputBase
from sage.repl.rich_output.buffer import OutputBuffer


class OutputVideoAny(OutputBase):

    def __init__(self, video, ext, mimetype, attrs = {}):
        """
        Video of any common video format

        If a backend claims support for this class, then it should
        accept files in common video formats and at least present
        video controls to the user, or open a video player.
        Due to the large number of video container formats, codecs,
        possible bit rate requirements and so on, it might well be
        that the video still won't be played, but the user should at
        least see a useful message about what is going on.

        INPUT:

        - ``video`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`.
          The video data.
        - ``ext`` -- string. The file name extension for this video format.
        - ``mimetype`` -- string. The MIME type of the video format.
        - ``attrs`` -- dict. Attributes for a ``<video>`` tag in HTML.
          Keys are strings, and values either strings or boolean values.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputVideoAny
            sage: OutputVideoAny.example()  # indirect doctest
            OutputVideoAny container
        """
        assert isinstance(video, OutputBuffer)
        self.video = video
        self.ext = ext
        self.mimetype = mimetype
        self.attrs = attrs

    @classmethod
    def example(cls):
        r"""
        Construct a sample video output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputVideoAny`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputVideoAny
            sage: OutputVideoAny.example()
            OutputVideoAny container
            sage: OutputVideoAny.example().video
            buffer containing 5540 bytes
            sage: OutputVideoAny.example().ext
            '.ogv'
            sage: OutputVideoAny.example().mimetype
            'video/ogg'
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.ogv')
        return cls(OutputBuffer.from_file(filename), '.ogv', 'video/ogg',
                   {'controls': True})

    def html_fragment(self, url, link_attrs=''):
        r"""
        Construct a HTML fragment for embedding this video

        INPUT:

        - ``url`` -- string. The URL where the data of this video can be found.

        - ``link_attrs`` -- string. Can be used to style the fallback link
          which is presented to the user if the video is not supported.

        EXAMPLES::
            sage: from sage.repl.rich_output.output_catalog import OutputVideoAny
            sage: print(OutputVideoAny.example().html_fragment
            ....:       ('foo', 'class="bar"').replace('><','>\n<'))
            <video controls="controls">
            <source src="foo" type="video/ogg" />
            <p>
            <a target="_new" href="foo" class="bar">Download video/ogg video</a>
            </p>
            </video>
        """
        attrs = dict((k, (k if v is True else v))
                     for k, v in self.attrs.iteritems()
                     if v is not False)
        attrs = ''.join(' {}="{}"'.format(k, v) for k, v in attrs.iteritems())
        return ('<video{attrs}>'
                '<source src="{url}" type="{mimetype}" /><p>'
                '<a target="_new" href="{url}" {link_attrs}>'
                'Download {mimetype} video</a></p></video>'
        ).format(url=url,
                 mimetype=self.mimetype,
                 attrs=attrs,
                 link_attrs=link_attrs,
        )
