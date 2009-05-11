r"""
Texture/material support for 3D Graphics objects and plotting.
This is a very rough common interface for Tachyon, x3d, and obj (mtl).

AUTHOR:
    -- Robert Bradshaw (2007-07-07) Initial version.

"""
from sage.structure.sage_object import SageObject

from sage.plot.misc import Color


uniq_c = 0

from sage.plot.misc import colors

def is_Texture(x):
    return isinstance(x, Texture_class)

def Texture(id=None, **kwds):
    if isinstance(id, Texture_class):
        return id
    if kwds.has_key('texture'):
        t = kwds['texture']
        if is_Texture(t):
            return t
        else:
            raise TypeError, "texture keyword must be a texture object"
    if isinstance(id, dict):
        kwds = id
        if kwds.has_key('rgbcolor'):
            kwds['color'] = kwds['rgbcolor']
        id = None
    elif isinstance(id, Color):
        kwds['color'] = id.rgb()
        id = None
    elif isinstance(id, str) and colors.has_key(id):
        kwds['color'] = id
        #kwds = {"color": id}
        id = None
    elif isinstance(id, tuple):
        kwds['color'] = id
        id = None
    if id is None:
        global uniq_c
        uniq_c += 1
        id = "texture%s" % uniq_c
    return Texture_class(id, **kwds)

def parse_color(info, base=None):
    if isinstance(info, Color):
        return info.rgb()
    elif isinstance(info, str):
        try:
            return colors[info]
        except KeyError:
            raise ValueError, "unknown color '%s'"%info
    else:
        return (float(info*base[0]), float(info*base[1]), float(info*base[2]))


class Texture_class(SageObject):
    """
    We create a translucent texture:

        sage: from sage.plot.plot3d.texture import Texture
        sage: t = Texture(opacity=0.6)
        sage: t
        Texture(texture..., 6666ff)
        sage: t.opacity
        0.600000000000000
        sage: t.jmol_str('obj')
        'color obj translucent 0.4 [102,102,255]'
        sage: t.mtl_str()
        'newmtl texture2\nKa 0.2 0.2 0.5\nKd 0.4 0.4 1.0\nKs 0.0 0.0 0.0\nillum 1\nNs 1\nd 0.600000000000000'
        sage: t.tachyon_str()
        'Texdef texture2\n  Ambient 0.333333333333 Diffuse 0.666666666667 Specular 0.0 Opacity 0.600000000000000\n   Color 0.4 0.4 1.0\n   TexFunc 0'
        sage: t.x3d_str()
        "<Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance>"
    """
    def __init__(self, id, color=(.4, .4, 1), opacity=1, ambient=0.5, diffuse=1, specular=0, shininess=1, name=None, **kwds):

        self.id = id
        if name is None and isinstance(color, str):
            name = color
        self.name = name

        if not isinstance(color, tuple):
            color = parse_color(color)
        else:
            if len(color) == 4:
                opacity = color[3]
            color = (float(color[0]), float(color[1]), float(color[2]))

        self.color = color
        self.opacity = opacity
        self.shininess = shininess

        if not isinstance(ambient, tuple):
            ambient = parse_color(ambient, color)
        self.ambient = ambient

        if not isinstance(diffuse, tuple):
            diffuse = parse_color(diffuse, color)
        self.diffuse = diffuse

        if not isinstance(specular, tuple):
            specular = parse_color(specular, color)
        self.specular = specular

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture('yellow')
            Texture(texture..., yellow, ffff00)
            sage: Texture((1,1,0), opacity=.5)
            Texture(texture..., ffff00)
        """
        if self.name is not None:
            return "Texture(%s, %s, %s)" % (self.id, self.name, self.hex_rgb())
        else:
            return "Texture(%s, %s)" % (self.id, self.hex_rgb())

    def hex_rgb(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.texture import Texture
            sage: Texture('red').hex_rgb()
            'ff0000'
            sage: Texture((1, .5, 0)).hex_rgb()
            'ff7f00'
        """
        return "%02x%02x%02x" % tuple(int(255*s) for s in self.color)

    def tachyon_str(self):
        total_color = float(sum(self.ambient) + sum(self.diffuse) + sum(self.specular))
        if total_color == 0:
            total_color = 1
        return "Texdef %s\n" % self.id + \
         "  Ambient %s Diffuse %s Specular %s Opacity %s\n" % \
                (sum(self.ambient)/total_color,
                 sum(self.diffuse)/total_color,
                 sum(self.specular)/total_color,
                 self.opacity) + \
        "   Color %s %s %s\n" % (self.color[0], self.color[1], self.color[2]) + \
        "   TexFunc 0"

    def x3d_str(self):
        return "<Appearance><Material diffuseColor='%s %s %s' shininess='%s' specularColor='%s %s %s'/></Appearance>" % \
                (self.color[0], self.color[1], self.color[2], self.shininess, self.specular[0], self.specular[0], self.specular[0])

    def mtl_str(self):
        return "\n".join(["newmtl %s" % self.id,
                   "Ka %s %s %s" % self.ambient,
                   "Kd %s %s %s" % self.diffuse,
                   "Ks %s %s %s" % self.specular,
                   "illum %s" % (2 if sum(self.specular) > 0 else 1),
                   "Ns %s" % self.shininess,
                   "d %s" % self.opacity, ])

    def jmol_str(self, obj):
        """
        EXAMPLES:
            sage: sum([dodecahedron(center=[2.5*x, 0, 0], color=(1, 0, 0, x/10)) for x in range(11)]).show(aspect_ratio=[1,1,1], frame=False, zoom=2)
        """
        translucent = "translucent %s" % float(1-self.opacity) if self.opacity < 1 else ""
        return "color %s %s [%s,%s,%s]" % (obj, translucent,
                int(255*self.color[0]), int(255*self.color[1]), int(255*self.color[2]))



