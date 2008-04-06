r"""
Texture/material support for 3D Graphics objects and plotting.
This is a very rough common interface for Tachyon, x3d, and obj (mtl).

AUTHOR:
    -- Robert Bradshaw (2007-07-07) Initial version.

"""
from sage.structure.sage_object import SageObject


uniq_c = 0

colors = {
    "red"   : (1.0,0.0,0.0),
    "orange": (1.0,.5,0.0),
    "yellow": (1.0,1.0,0.0),
    "green" : (0.0,1.0,0.0),
    "blue"  : (0.0,0.0,1.0),
    "purple": (.5,0.0,1.0),
    "white" : (1.0,1.0,1.0),
    "black" : (0.0,0.0,0.0),
    "grey"  : (.5,.5,.5)
}

def Texture(id=None, **kwds):
    if isinstance(id, Texture_class):
        return id
    if isinstance(id, dict):
        kwds = id
        id = None
    elif isinstance(id, str) and colors.has_key(id):
        kwds = {"color": id}
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
        if isinstance(info, str):
            try:
                return colors[info]
            except KeyError:
                raise # TODO: parse hex?
        else:
            return (float(info*base[0]), float(info*base[1]), float(info*base[2]))


class Texture_class(SageObject):

    def __init__(self, id, color=(.5, .5, .5), opacity=1, ambient=0.5, diffuse=1, specular=0, shininess=1):
        self.id = id

        if not isinstance(color, tuple):
            color = parse_color(color)
        else:
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


    def tachyon_str(self):
        total_color = float(sum(self.ambient) + sum(self.diffuse) + sum(self.specular))
        if total_color == 0:
            total_color = 1
        return "    Texture Ambient %s Diffuse %s Specular %s Opacity %s\n" % \
                (sum(self.ambient)/total_color,
                 sum(self.diffuse)/total_color,
                 sum(self.specular)/total_color,
                 self.opacity) + \
        "       Color %s %s %s\n" % (self.color[0], self.color[1], self.color[2]) + \
        "       TexFunc 0"

    def x3d_str(self):
        return "<Appearance><Material diffuseColor='%s %s %s' shininess='%s' specularColor='%s %s %s'/></Appearance>" % \
                (self.color[0], self.color[1], self.color[2], self.shininess, self.specular[0], self.specular[0], self.specular[0])

    def mtl_str(self):
        return "\n".join(["newmtl %s" % self.id,
                   "Ka %s %s %s" % self.ambient,
                   "Kd %s %s %s" % self.diffuse,
                   "Ks %s %s %s" % self.specular,
                   "illum %s" % (2 if sum(self.specular) > 0 else 1),
                   "Tr %s" % self.opacity,
                   "Ns %s" % self.shininess ])

