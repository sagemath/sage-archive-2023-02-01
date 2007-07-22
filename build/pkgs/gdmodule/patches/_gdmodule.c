/******************************************************************
Copyright 1995 Richard Jones, Bureau of Meteorology Australia.
richard@bofh.asn.au

Current maintainer is
    Chris Gonnerman <chris.gonnerman@newcenturycomputers.net>
Please direct all questions and problems to me.

This module is a python wrapper for the GD library (version 1.8.3)

version 0.56
Revised 03/10/2005 by Chris Gonnerman
    -- added support for GIF files.
    -- the minimum gd library version for this gdmodule version is
       2.0.23.

version 0.55
Revised 03/10/2005 by Chris Gonnerman
    -- corrected error in the gd.image() constructor, pointed out
       by Betty Li.  maybe this time I got it right.
    -- the minimum gd library version for this gdmodule version is
       2.0.22.

version 0.54
Revised 03/03/2005 by Chris Gonnerman
    -- corrected error in the gd.image() constructor, pointed out
       by Betty Li.
    -- corrected yet another obvious error reported by Greg Hewgill,
       regarding an incorrect call to gdImagePolygon when fillcolor
       is not equal to 1.
    -- implemented the saveAlpha and alphaBlending features as
       suggested by Christopher Stone.
    -- fixed memory management issue regarding image deallocation
       reported by Matti Jagula.
    -- fixed memory management issue in image_filledpolygon
       reported by Sadruddin Rejeb.
    -- the minimum gd library version for this gdmodule version is
       2.0.22.

version 0.53
Revised 06/10/2004 by Chris Gonnerman
    -- corrected obvious error in image_settile(), as pointed out
       by Dan Mosedale
    -- applied some patches provided by John Hunter.  Several are
       for Windows compatibility, but this should be considered
       a "work in progress" in this version.
    -- corrected another obvious error (and an ancient one too)
       reported by Greg Hewgill, regarding the foreground color
       being used rather than the fill color in the image_polygon
       function.
    -- the minimum gd library version for this gdmodule version is
       2.0.22.

version 0.52
Revised 02/03/2004 by Chris Gonnerman
    -- incorporated new functions from the matplotlib project
       (provided by John Hunter)
    -- corrected error in PyObject_HEAD_INIT() call noted by
       Stefan R Kuzminski (thanks!)

version 0.51
Revised 09/24/2003 by Chris Gonnerman
    -- fixed memory management bug reported by Jack Diederich.

version 0.50
Revised 09/13/2003 by Chris Gonnerman
    -- documentation change only: the GD library must be version
       2.0.8 or higher (or you must choose version 0.40 or earlier
       of this module, which supports GD 1.8.3 or higher).

version 0.42
Revised 06/10/2003 by Chris Gonnerman
    -- implemented patch by Gregory P. Smith to support writing to
       arbitrary objects, superceding the patch in 0.41 from
       Andreas Rottman.  gddemo.py contains Rottman's example
       (which works with Smith's patched version) as well as an
       example (by me) of reading an image from a URL.  The
       documentation was updated to reflect these changes.

version 0.41
Revised 05/29/2003 by Chris Gonnerman
    -- implemented big patch by Andreas Rottmann to support writing
       images to memory via StringIO/CStringIO.  Currently the only
       documentation is an example in gddemo.py.
    -- implemented patch by Bob Galloway to remove the "specialness"
       of negative coordinates in the string_ttf/string_ft methods.
    -- implemented patch by Nathan Robertson to enable MacOSX fink
       builds.
    -- fixed bugs in the proxy class with regard to passing of
       _gd.image types.

version 0.40
Revised 09/18/2002 by Chris Gonnerman
    -- updated to deal with incomplete library installs

version 0.30
Revised 06/12/2002 by Chris Gonnerman
    -- initial conversion to two-tier design

version 0.26
Revised 12/10/2001 by Chris Gonnerman
    -- made unicode optional via define Py_UNICODEOBJECT_H

version 0.25
Revised 09/15/2001 by Chris Gonnerman
    -- implemented patch by Mike Romberg:
        adds additional public functions, available
        only if GD 2.0.1+ is installed.

version 0.24
Revised 08/08/2001 by Chris Gonnerman
    -- implemented patch by Mike Romberg:
        supports GCC 3.0 and GD 2.0.1
        adds public C function makeGDImage()
            It allows C or C++ code which uses the gd library
            directly to make an imageobject instance.

version 0.23
Revised 11/22/2000 by Chris Gonnerman
    -- included the patch from Tanimoto Osamu
    -- updated support to GD version 1.8.3
    -- cleaned up code, added ANSI prototypes

******************************************************************/

#include <Python.h>
#include "gd.h"
#include "gdfonts.h"
#include "gdfontl.h"
#include "gdfontmb.h"
#include "gdfontt.h"
#include "gdfontg.h"
#include <string.h>
#include <errno.h>

#ifdef HAVE_LIBTTF
#define HAVE_LIBFREETYPE
#endif

/* missing from gd.h */
/* gdImagePtr gdImageCreateFromXpm(char *filename); */

#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif



static PyObject *ErrorObject;
extern int errno;

/*
** Declarations for objects of type image
*/

typedef struct i_o {
    PyObject_HEAD
    gdImagePtr imagedata;
    int multiplier_x,origin_x;
    int multiplier_y,origin_y;
    struct i_o *current_brush;
    struct i_o *current_tile;
} imageobject;


// 2.0.22 replaces gdFontTinyRep with gdFontGetTiny, and so on
typedef gdFontPtr (*ptrGetFontFunction)();
typedef struct {
    char *name;
    ptrGetFontFunction func;
} fontstruct;

static fontstruct fonts[] = {
  {"gdFontTiny",&gdFontGetTiny},
  {"gdFontSmall",&gdFontGetSmall},
  {"gdFontMediumBold",&gdFontGetMediumBold},
  {"gdFontLarge",&gdFontGetLarge},
  {"gdFontGiant",&gdFontGetGiant},
  {NULL,NULL}
};

staticforward PyTypeObject Imagetype;

#define is_imageobject(v)        ((v)->ob_type == &Imagetype)

#define MIN(x,y) ((x)<(y)?(x):(y))
#define X(x) ((x)*self->multiplier_x+self->origin_x)
#define Y(y) ((y)*self->multiplier_y+self->origin_y)
#define W(x) ((x)*self->multiplier_x)
#define H(y) ((y)*self->multiplier_y)

static imageobject *newimageobject(PyObject *args);

/*
** Support Functions
*/

static PyObject *write_file(imageobject *img, PyObject *args, char fmt)
{
    char *filename;
    PyObject *fileobj;
    FILE *fp = NULL;
    int closeme = 0, use_fileobj_write = 0;
    int arg1 = -1, arg2 = -1;
    int filesize = 0;
    void *filedata = NULL;

    if(PyArg_ParseTuple(args, "O!|ii", &PyFile_Type, &fileobj, &arg1, &arg2)) {
        fp = PyFile_AsFile(fileobj);
    } else if(PyErr_Clear(), PyArg_ParseTuple(args, "z|ii", &filename, &arg1, &arg2)) {
        if((fp = fopen(filename, "wb"))) {
            closeme = 1;
        } else {
            PyErr_SetFromErrno(PyExc_IOError);
            return NULL;
        }
    } else if (PyErr_Clear(), PyArg_ParseTuple(args, "O|ii", &fileobj, &arg1, &arg2)) {
        /* if we're passed a random object that has a write method we will
         * attempt to call object.write(filedata) */
        if (!PyObject_HasAttrString(fileobj, "write")) {
            PyErr_SetString(ErrorObject, "first argument must be a file, string or object with a write method");
            return NULL;
        }
        use_fileobj_write = 1;
    }
    else
        return NULL;

    switch(fmt) {
    case 'f' : /* gif */
#ifdef HAVE_LIBGIF
        if (use_fileobj_write) {
            filedata = gdImageGifPtr(img->imagedata, &filesize);
        } else {
            gdImageGif(img->imagedata, fp);
        }
#else
        PyErr_SetString(PyExc_NotImplementedError,
                    "GIF Support Not Available");
        return NULL;
#endif
        break;
    case 'p' : /* png */
#ifdef HAVE_LIBPNG
        if (use_fileobj_write) {
            filedata = gdImagePngPtr(img->imagedata, &filesize);
        } else {
            gdImagePng(img->imagedata, fp);
        }
#else
        PyErr_SetString(PyExc_NotImplementedError,
                    "PNG Support Not Available");
        return NULL;
#endif
        break;
    case 'j' : /* jpeg */
#ifdef HAVE_LIBJPEG
        if (use_fileobj_write) {
            filedata = gdImageJpegPtr(img->imagedata, &filesize, arg1);
        } else {
            gdImageJpeg(img->imagedata, fp, arg1);
        }
#else
        PyErr_SetString(PyExc_NotImplementedError,
                    "PNG Support Not Available");
        return NULL;
#endif
        break;
    case 'g' : /* gd */
        if (use_fileobj_write) {
            filedata = gdImageGdPtr(img->imagedata, &filesize);
        } else {
            gdImageGd(img->imagedata, fp);
        }
        break;
    case 'G' : /* gd2 */
        if(arg1 == -1) arg1 = 0;
        if(arg2 != GD2_FMT_RAW && arg2 != GD2_FMT_COMPRESSED)
            arg2 = GD2_FMT_COMPRESSED;
        if (use_fileobj_write) {
            filedata = gdImageGd2Ptr(img->imagedata, arg1, arg2, &filesize);
        } else {
            gdImageGd2(img->imagedata, fp, arg1, arg2);
        }
        break;
    case 'w' : /* wbmp */
        if(arg1 == -1)
            arg1 = 0;
        if (use_fileobj_write) {
	  //filedata = gdImageWBMPPtr(img->imagedata, &filesize, arg1);
        } else {
            gdImageWBMP(img->imagedata, arg1, fp);
        }
        break;
    }

    if (use_fileobj_write || filedata) {
        PyObject *noerr;
        noerr = PyObject_CallMethod(fileobj, "write", "s#", filedata, filesize);
        gdFree(filedata);
        if (noerr == NULL)
            return NULL;
    } else if (closeme) {
        fclose(fp);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/*
** Methods for the image type
*/

imageobject *makeGDImage(gdImagePtr imagedata)
    {
    gdImagePtr newimg = gdImageCreate(gdImageSX(imagedata),
      gdImageSY(imagedata));
    imageobject *rval = 0;
    gdImageCopy(newimg, imagedata, 0, 0, 0, 0, gdImageSX(imagedata),
      gdImageSY(imagedata));

    if (!(rval = PyObject_NEW(imageobject, &Imagetype)))
        return NULL;

    rval->current_tile = rval->current_brush = NULL;
    rval->origin_x = rval->origin_y = 0;
    rval->multiplier_x = rval->multiplier_y = 1;
    rval->imagedata = newimg;
    return rval;
    }


/*** I/O Methods ***/

static PyObject *image_writegif(imageobject *self, PyObject *args)
{
    return write_file(self, args, 'f');
}


static PyObject *image_writepng(imageobject *self, PyObject *args)
{
    return write_file(self, args, 'p');
}

#ifdef HAVE_LIBJPEG
static PyObject *image_writejpeg(imageobject *self, PyObject *args)
{
    return write_file(self, args, 'j');
}
#endif

static PyObject *image_writegd(imageobject *self, PyObject *args)
{
    return write_file(self, args, 'g');
}


static PyObject *image_writegd2(imageobject *self, PyObject *args)
{
    return write_file(self, args, 'G');
}


static PyObject *image_writewbmp(imageobject *self, PyObject *args)
{
    return write_file(self, args, 'w');
}


/*** Drawing Methods ***/

static PyObject *image_setpixel(imageobject *self, PyObject *args)
{
    int x,y,color;

    if(!PyArg_ParseTuple(args, "(ii)i", &x, &y, &color))
        return NULL;

    gdImageSetPixel(self->imagedata, X(x), Y(y), color);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_line(imageobject *self, PyObject *args)
{
    int sx,sy,ex,ey,color;

    if(!PyArg_ParseTuple(args, "(ii)(ii)i", &sx, &sy, &ex, &ey, &color))
        return NULL;
    gdImageLine(self->imagedata, X(sx), Y(sy), X(ex), Y(ey), color);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *image_lines(imageobject *self, PyObject *args)
{
    int color,i,N;
    long sx,sy,ex,ey;
    PyObject *seq;
    PyObject *lastTup, *thisTup;

    if(!PyArg_ParseTuple(args, "Oi", &seq, &color))
        return NULL;

    seq = PySequence_Fast(seq, NULL);
    N = PySequence_Length(seq);
    if (N<2) {
      PyErr_SetString(PyExc_ValueError,
                    "lines() requires sequence of len(2) or greater");
      return NULL;
    }

    lastTup = PySequence_GetItem(seq, 0);
    sx = X(PyInt_AsLong(PySequence_GetItem(lastTup, 0)));
    sy = Y(PyInt_AsLong(PySequence_GetItem(lastTup, 1)));

    for (i=0; i<N; ++i) {
      thisTup = PySequence_GetItem(seq, i);
      ex = X(PyInt_AsLong(PySequence_GetItem(thisTup, 0)));
      ey = Y(PyInt_AsLong(PySequence_GetItem(thisTup, 1)));
      gdImageLine(self->imagedata, sx, sy, ex, ey, color);
      sx = ex;
      sy = ey;
   }

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_polygon(imageobject *self, PyObject *args)
{
    PyObject *point, *points;
    gdPointPtr gdpoints;
    int size, color, i, fillcolor = -1;

    if(!PyArg_ParseTuple(args, "O!i|i", &PyTuple_Type, &points, &color, &fillcolor)) {
        PyErr_Clear();
        if(PyArg_ParseTuple(args, "O!i|i", &PyList_Type, &points, &color, &fillcolor))
            points=PyList_AsTuple(points);
        else
            return NULL;
    }

    size = PyTuple_Size(points);

    gdpoints = (gdPointPtr)calloc(size,sizeof(gdPoint));

    for(i=0; i<size; i++) {
        point = PyTuple_GET_ITEM(points,i);
        gdpoints[i].x = X(PyInt_AS_LONG((PyIntObject *)PyTuple_GET_ITEM(point,0)));
        gdpoints[i].y = Y(PyInt_AS_LONG((PyIntObject *)PyTuple_GET_ITEM(point,1)));
    }

    if(fillcolor != -1)
        gdImageFilledPolygon(self->imagedata, gdpoints, size, fillcolor);

    gdImagePolygon(self->imagedata, gdpoints, size, color);

    free(gdpoints);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_rectangle(imageobject *self, PyObject *args)
{
    int tx,ty,bx,by,t,color,fillcolor,fill=0;

    if(PyArg_ParseTuple(args, "(ii)(ii)ii", &tx, &ty, &bx, &by, &color, &fillcolor))
        fill=1;
    else if(PyErr_Clear(), !PyArg_ParseTuple(args, "(ii)(ii)i", &tx, &ty, &bx, &by, &color))
        return NULL;

    tx = X(tx); ty = Y(ty);
    bx = X(bx); by = Y(by);

    if(tx > bx) {
        t = tx;
        tx = bx;
        bx = t;
    }

    if(ty > by) {
        t = ty;
        ty = by;
        by = t;
    }

    if(fill)
        gdImageFilledRectangle(self->imagedata, tx, ty, bx, by, fillcolor);

    gdImageRectangle(self->imagedata, tx, ty, bx, by, color);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_filledpolygon(imageobject *self, PyObject *args)
{
    PyObject *point, *points;
    gdPointPtr gdpoints;
    int size, color, i;

    if(!PyArg_ParseTuple(args, "O!i", &PyTuple_Type, &points, &color)) {
        PyErr_Clear();
        if(PyArg_ParseTuple(args, "O!i", &PyList_Type, &points, &color))
            points=PyList_AsTuple(points);
        else
            return NULL;
    }

    size = PyTuple_Size(points);
    gdpoints = (gdPointPtr)calloc(size,sizeof(gdPoint));

    for(i = 0; i < size; i++) {
        point = PyTuple_GET_ITEM(points,i);
        gdpoints[i].x = X(PyInt_AS_LONG((PyIntObject *)PyTuple_GET_ITEM(point,0)));
        gdpoints[i].y = Y(PyInt_AS_LONG((PyIntObject *)PyTuple_GET_ITEM(point,1)));
    }
    gdImageFilledPolygon(self->imagedata, gdpoints, size, color);
    free(gdpoints);

    Py_DECREF(points);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_filledrectangle(imageobject *self, PyObject *args)
{
    int tx,ty,bx,by,t,color;

    if(!PyArg_ParseTuple(args, "(ii)(ii)i", &tx, &ty, &bx, &by, &color))
        return NULL;
    tx = X(tx); ty = Y(ty);
    bx = X(bx); by = Y(by);
    if(tx > bx) {t=tx;tx=bx;bx=t;}
    if(ty > by) {t=ty;ty=by;by=t;}
    gdImageFilledRectangle(self->imagedata, tx, ty, bx, by, color);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_arc(imageobject *self, PyObject *args)
{
    int cx, cy, w, h, s, e, color, i;

    if(!PyArg_ParseTuple(args, "(ii)(ii)iii", &cx, &cy, &w, &h, &s, &e, &color))
        return NULL;
    if(e<s) {i=e;e=s;s=i;}
    gdImageArc(self->imagedata, X(cx), Y(cy), W(w), H(h), s, e, color);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *image_filledarc(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
    PyErr_SetString(PyExc_NotImplementedError,
                    "filledArc() requires gd 2.0 or later");
    return NULL;
#else
    int cx, cy, w, h, s, e, color, i, style;

    if(!PyArg_ParseTuple(args, "(ii)(ii)iiii", &cx, &cy, &w, &h, &s,
                         &e, &color, &style))
        return NULL;
    if(e<s) {i=e;e=s;s=i;}
    gdImageFilledArc(self->imagedata, X(cx), Y(cy), W(w), H(h), s, e,
                     color, style);
    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_filledellipse(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
    PyErr_SetString(PyExc_NotImplementedError,
                    "filledEllipse() requires gd 2.0 or later");
    return NULL;
#else
    int cx, cy, w, h, color;

    if(!PyArg_ParseTuple(args, "(ii)(ii)i", &cx, &cy, &w, &h, &color))
        return NULL;
    gdImageFilledEllipse(self->imagedata, X(cx), Y(cy), W(w), H(h), color);
    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_filltoborder(imageobject *self, PyObject *args)
{
    int x,y,border,color;

    if(!PyArg_ParseTuple(args, "(ii)ii", &x,&y,&border,&color))
        return NULL;
    gdImageFillToBorder(self->imagedata, X(x),Y(y),border,color);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_fill(imageobject *self, PyObject *args)
{
    int x,y,color;

    if(!PyArg_ParseTuple(args, "(ii)i", &x,&y,&color))
        return NULL;
    gdImageFill(self->imagedata, X(x),Y(y),color);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_setbrush(imageobject *self, PyObject *args)
{
    imageobject *brush;
    char *filename, *type; /* dummies */

    if(PyArg_ParseTuple(args, "z|z", &filename, &type))
        brush = (imageobject *)newimageobject(args);
    else if(PyErr_Clear(), PyArg_ParseTuple(args, "O!", &Imagetype, &brush))
        Py_INCREF(brush);
    else
        return NULL;
    if(self->current_brush){
        Py_DECREF(self->current_brush);
    }
    self->current_brush = brush;
    gdImageSetBrush(self->imagedata, brush->imagedata);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_setantialiased(imageobject *self, PyObject *args)
{
    int c;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;

    gdImageSetAntiAliased(self->imagedata, c);

    Py_INCREF(Py_None);
    return Py_None;
}



static PyObject *image_setclip(imageobject *self, PyObject *args)
{
    int tx,ty,bx,by,t;

    if(!PyArg_ParseTuple(args, "(ii)(ii)", &tx, &ty, &bx, &by))
        return NULL;

    tx = X(tx); ty = Y(ty);
    bx = X(bx); by = Y(by);

    if(tx > bx) {
        t = tx;
        tx = bx;
        bx = t;
    }

    if(ty > by) {
        t = ty;
        ty = by;
        by = t;
    }

    gdImageSetClip(self->imagedata, tx, ty, bx, by);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_settile(imageobject *self, PyObject *args)
{
    imageobject *tile;
    char *filename, *type; /* dummies */

    if(PyArg_ParseTuple(args, "z|z", &filename, &type))
        tile = (imageobject *)newimageobject(args);
    else
        if(PyErr_Clear(), PyArg_ParseTuple(args, "O!", &Imagetype, &tile))
            Py_INCREF(tile);
        else
            return NULL;

    if(self->current_tile) {
        Py_DECREF(self->current_tile);
    }

    self->current_tile = tile;
    gdImageSetTile(self->imagedata, tile->imagedata);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_setstyle(imageobject *self, PyObject *args)
{
    PyObject *style;
    int size, i, *stylearray;

    if(!PyArg_ParseTuple(args, "O!", &PyTuple_Type, &style)) {
        PyErr_Clear();
        if(PyArg_ParseTuple(args, "O!", &PyList_Type, &style))
            style=PyList_AsTuple(style);
        else
            return NULL;
    }

    size = PyTuple_Size(style);

    stylearray = (int *)calloc(size,sizeof(int));

    for(i=0; i<size; i++)
        stylearray[i] = PyInt_AS_LONG((PyIntObject *)PyTuple_GET_ITEM(style,i));

    gdImageSetStyle(self->imagedata, stylearray, size);
    free(stylearray);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *image_setthickness(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
    PyErr_SetString(PyExc_NotImplementedError,
                    "setThickness() requires gd 2.0 or later");
    return NULL;
#else
    int t;

    if(!PyArg_ParseTuple(args, "i", &t))
        return NULL;
    gdImageSetThickness(self->imagedata, t);
    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_alphablending(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
    PyErr_SetString(PyExc_NotImplementedError,
                    "alphaBlending() requires gd 2.0 or later");
    return NULL;
#else
    int blending;

    if(!PyArg_ParseTuple(args, "i", &blending))
        return NULL;
    gdImageAlphaBlending(self->imagedata, blending);
    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_alpha(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
    PyErr_SetString(PyExc_NotImplementedError,
                    "alpha() requires gd 2.0 or later");
    return NULL;
#else
    int color;

    if(!PyArg_ParseTuple(args, "i", &color))
        return NULL;
    return Py_BuildValue("i", gdImageAlpha(self->imagedata, color));
#endif
}


static PyObject *image_getpixel(imageobject *self, PyObject *args)
{
    int x,y;

    if(!PyArg_ParseTuple(args, "(ii)", &x,&y))
        return NULL;
    return Py_BuildValue("i",gdImageGetPixel(self->imagedata, X(x),Y(y)));
}


static PyObject *image_boundssafe(imageobject *self, PyObject *args)
{
    int x,y;

    if(!PyArg_ParseTuple(args, "(ii)", &x,&y))
        return NULL;
    return Py_BuildValue("i",gdImageBoundsSafe(self->imagedata, X(x),Y(y)));
}


static PyObject *image_blue(imageobject *self, PyObject *args)
{
    int c;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;
    return Py_BuildValue("i",gdImageBlue(self->imagedata,c));
}


static PyObject *image_green(imageobject *self, PyObject *args)
{
    int c;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;
    return Py_BuildValue("i",gdImageGreen(self->imagedata,c));
}


static PyObject *image_red(imageobject *self, PyObject *args)
{
    int c;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;
    return Py_BuildValue("i",gdImageRed(self->imagedata,c));
}


static PyObject *image_size(imageobject *self)
{
    return Py_BuildValue("(ii)",gdImageSX(self->imagedata),gdImageSY(self->imagedata));
}


static PyObject *image_char(imageobject *self, PyObject *args)
{
    int x,y,font,color;
    char c;

    if(!PyArg_ParseTuple(args, "i(ii)ii", &font,&x,&y,&c,&color))
        return NULL;
    gdImageChar(self->imagedata, fonts[font].func(), X(x), Y(y), c, color);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_charup(imageobject *self, PyObject *args)
{
    int x,y,font,color;
    char c;

    if(!PyArg_ParseTuple(args, "i(ii)si", &font,&x,&y,&c,&color))
        return NULL;
    gdImageCharUp(self->imagedata, fonts[font].func(), X(x), Y(y), c, color);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_get_bounding_rect(imageobject *self, PyObject *args)
{
#ifndef HAVE_LIBFREETYPE
    PyErr_SetString(PyExc_NotImplementedError,
                    "Freetype Support Not Available");
    return NULL;
#else
    double ptsize, angle;
    char *fontname, *str, *rc;
    int x, y, brect[8];

    if(!PyArg_ParseTuple(args, "sdd(ii)s",
        &fontname, &ptsize, &angle, &x, &y, &str))
            return NULL;

    rc = gdImageStringTTF(NULL, brect, 0,
        fontname, ptsize, angle, x, y, str);

    if(rc != NULL){
        PyErr_SetString(PyExc_ValueError, rc);
        return NULL;
    }

    return Py_BuildValue("(iiiiiiii)",brect[0],brect[1],brect[2],brect[3],
                 brect[4],brect[5],brect[6],brect[7]);
#endif
}


static PyObject *image_string_ft(imageobject *self, PyObject *args)
{
#ifndef HAVE_LIBFREETYPE
    PyErr_SetString(PyExc_NotImplementedError,
                    "Freetype Support Not Available");
    return NULL;
#else
#if GD2_VERS <= 1
    PyErr_SetString(PyExc_NotImplementedError,
                    "string_ft() requires gd 2.0 or later");
    return NULL;
#else
    int x,y,fg;
    double ptsize,angle;
    char *fontname,*str,*rc;
    int brect[8];

    if(!PyArg_ParseTuple(args, "sdd(ii)si",
        &fontname, &ptsize, &angle, &x,&y,&str,&fg))
            return NULL;
    rc = gdImageStringFT(NULL, brect, 0, fontname, ptsize, angle,
                  0, 0, str); /* calculate the bounding rect */
    if(rc != NULL){
      PyErr_SetString(PyExc_ValueError, rc);
      return NULL;
    }
    rc = gdImageStringTTF(self->imagedata, brect, fg, fontname,
            ptsize, angle, x, y, str);
    if(rc != NULL){
        PyErr_SetString(PyExc_ValueError, rc);
        return NULL;
    }

    return Py_BuildValue("(iiiiiiii)",
        brect[0], brect[1], brect[2], brect[3],
        brect[4], brect[5], brect[6], brect[7]);
#endif
#endif
}


static PyObject *image_string_ttf(imageobject *self, PyObject *args)
{
#ifndef HAVE_LIBFREETYPE
    PyErr_SetString(PyExc_NotImplementedError,
                    "Freetype Support Not Available");
    return NULL;
#else
    int x,y,fg;
    double ptsize,angle;
    char *fontname,*str,*rc;
    int brect[8];

    if(!PyArg_ParseTuple(args, "sdd(ii)si",
        &fontname, &ptsize, &angle, &x,&y,&str,&fg))
            return NULL;
    rc = gdImageStringTTF(NULL, brect, 0, fontname, ptsize, angle,
                  0, 0, str); /* calculate the bounding rect */
    if(rc != NULL){
      PyErr_SetString(PyExc_ValueError, rc);
      return NULL;
    }
    rc = gdImageStringTTF(self->imagedata, brect, fg, fontname,
            ptsize, angle, x, y, str);
    if(rc != NULL){
        PyErr_SetString(PyExc_ValueError, rc);
        return NULL;
    }

    return Py_BuildValue("(iiiiiiii)",
        brect[0], brect[1], brect[2], brect[3],
        brect[4], brect[5], brect[6], brect[7]);
#endif
}


static PyObject *image_string(imageobject *self, PyObject *args)
{
    int x,y,font,color;
    char *str;

    if(!PyArg_ParseTuple(args, "i(ii)si", &font,&x,&y,&str,&color))
        return NULL;
    gdImageString(self->imagedata, fonts[font].func(), X(x), Y(y), str, color);

    Py_INCREF(Py_None);
    return Py_None;
}

#ifdef Py_UNICODEOBJECT_H
static PyObject *image_string16(imageobject *self, PyObject *args)
{
    int x,y,font,color;
    Py_UNICODE *ustr;

    if(!PyArg_ParseTuple(args, "i(ii)ui", &font,&x,&y,&ustr,&color))
        return NULL;
    gdImageString16(self->imagedata, fonts[font].func(), X(x), Y(y), ustr,
                    color);

    Py_INCREF(Py_None);
    return Py_None;
}
#endif

static PyObject *image_stringup(imageobject *self, PyObject *args)
{
    int x,y,font,color;
    char *str;

    if(!PyArg_ParseTuple(args, "i(ii)si", &font,&x,&y,&str,&color))
        return NULL;
    gdImageStringUp(self->imagedata, fonts[font].func(), X(x), Y(y), str, color);

    Py_INCREF(Py_None);
    return Py_None;
}

#ifdef Py_UNICODEOBJECT_H
static PyObject *image_stringup16(imageobject *self, PyObject *args)
{
    int x,y,font,color;
    Py_UNICODE *ustr;

    if(!PyArg_ParseTuple(args, "i(ii)ui", &font,&x,&y,&ustr,&color))
        return NULL;
    gdImageStringUp16(self->imagedata, fonts[font].func(), X(x), Y(y),
                      ustr, color);

    Py_INCREF(Py_None);
    return Py_None;
}
#endif

static PyObject *image_colorallocate(imageobject *self, PyObject *args)
{
    int r,g,b;

    if(!PyArg_ParseTuple(args, "(iii)", &r, &g, &b))
        return NULL;
    return(Py_BuildValue("i",gdImageColorAllocate(self->imagedata, r, g, b)));
}

static PyObject *image_colorallocatealpha(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "colorAllocateAlpha() requires gd 2.0 or later");
    return NULL;
#else
    int r,g,b,a;

    if(!PyArg_ParseTuple(args, "(iiii)", &r, &g, &b, &a))
        return NULL;
    return(Py_BuildValue("i",gdImageColorAllocateAlpha(self->imagedata,
      r, g, b, a)));
#endif
}

static PyObject *image_colorclosest(imageobject *self, PyObject *args)
{
    int r,g,b;

    if(!PyArg_ParseTuple(args, "(iii)", &r, &g, &b))
        return NULL;
    return(Py_BuildValue("i",gdImageColorClosest(self->imagedata, r, g, b)));
}

static PyObject *image_colorclosestalpha(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "colorClosestAlpha() requires gd 2.0 or later");
    return NULL;
#else
    int r,g,b,a;

    if(!PyArg_ParseTuple(args, "(iiii)", &r, &g, &b, &a))
        return NULL;
    return(Py_BuildValue("i",gdImageColorClosestAlpha(self->imagedata, r, g, b,
      a)));
#endif
}

static PyObject *image_colorclosestHWB(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "colorClosestHWB() requires gd 2.0 or later");
    return NULL;
#else
    int r,g,b;

    if(!PyArg_ParseTuple(args, "(iii)", &r, &g, &b))
        return NULL;
    return(Py_BuildValue("i",gdImageColorClosestHWB(self->imagedata, r,
      g, b)));
#endif
}

static PyObject *image_colorexact(imageobject *self, PyObject *args)
{
    int r,g,b;

    if(!PyArg_ParseTuple(args, "(iii)", &r, &g, &b))
        return NULL;
    return(Py_BuildValue("i",gdImageColorExact(self->imagedata, r, g, b)));
}

static PyObject *image_colorresolve(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "colorResolve() requires gd 2.0 or later");
    return NULL;
#else
    int r,g,b;

    if(!PyArg_ParseTuple(args, "(iii)", &r, &g, &b))
        return NULL;
    return(Py_BuildValue("i",gdImageColorResolve(self->imagedata, r,
      g, b)));
#endif
}

static PyObject *image_colorresolvealpha(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "colorResolveAlpha() requires gd 2.0 or later");
    return NULL;
#else
    int r,g,b,a;

    if(!PyArg_ParseTuple(args, "(iiii)", &r, &g, &b, &a))
        return NULL;
    return(Py_BuildValue("i",gdImageColorResolveAlpha(self->imagedata, r,
      g, b, a)));
#endif
}

static PyObject *image_colorstotal(imageobject *self)
{
    return Py_BuildValue("i",gdImageColorsTotal(self->imagedata));
}


static PyObject *image_colorcomponents(imageobject *self, PyObject *args)
{
    int c,r,g,b;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;
    r=gdImageRed(self->imagedata,c);
    g=gdImageGreen(self->imagedata,c);
    b=gdImageBlue(self->imagedata,c);
    return Py_BuildValue("(iii)",r,g,b);
}


static PyObject *image_getclip(imageobject *self)
{
    int x1, y1, x2, y2;
    gdImageGetClip(self->imagedata, &x1, &y1, &x2, &y2);
    return Py_BuildValue("(ii)(ii)", x1, y1, x2, y2);
}

static PyObject *image_getinterlaced(imageobject *self)
{
    return Py_BuildValue("i",gdImageGetInterlaced(self->imagedata));
}


static PyObject *image_gettransparent(imageobject *self)
{
    return Py_BuildValue("i",gdImageGetTransparent(self->imagedata));
}


static PyObject *image_colordeallocate(imageobject *self, PyObject *args)
{
    int c;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;
    gdImageColorDeallocate(self->imagedata, c);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_colortransparent(imageobject *self, PyObject *args)
{
    int c;

    if(!PyArg_ParseTuple(args, "i", &c))
        return NULL;
    gdImageColorTransparent(self->imagedata, c);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_copyto(imageobject *self, PyObject *args)
{
    imageobject *dest;
    int dx,dy,sx,sy,w,h,dw,dh;

    dx=dy=sx=sy=0;
    w = gdImageSX(self->imagedata);
    h = gdImageSY(self->imagedata);
    if(!PyArg_ParseTuple(args, "O!|(ii)(ii)(ii)", &Imagetype, &dest, &dx, &dy, &sx, &sy, &w, &h))
        return NULL;
    dw = gdImageSX(dest->imagedata);
    dh = gdImageSY(dest->imagedata);
    gdImageCopy(dest->imagedata, self->imagedata, X(dx), Y(dy), X(sx), Y(sy), W(w), H(h));

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_copyresizedto(imageobject *self, PyObject *args)
{
    imageobject *dest;
    int dx,dy,sx,sy,dw,dh,sw,sh;

    dx=dy=sx=sy=0;
    sw = gdImageSX(self->imagedata);
    sh = gdImageSY(self->imagedata);
    if(PyArg_ParseTuple(args, "O!|(ii)(ii)", &Imagetype, &dest, &dx, &dy, &sx, &sy))
    {
        dw = gdImageSX(dest->imagedata);
        dh = gdImageSY(dest->imagedata);
    }
    else if(PyErr_Clear(), !PyArg_ParseTuple(args, "O!|(ii)(ii)(ii)(ii)", &Imagetype, &dest, &dx, &dy, &sx, &sy, &dw, &dh, &sw, &sh))
        return NULL;
    gdImageCopyResized(dest->imagedata, self->imagedata, X(dx), Y(dy), X(sx), Y(sy), W(dw), H(dh), W(sw), H(sh));

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_copyresampledto(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "copyResampledTo() requires gd 2.0 or later");
    return NULL;
#else
    imageobject *dest;
    int dx,dy,sx,sy,dw,dh,sw,sh;

    dx=dy=sx=sy=0;
    sw = gdImageSX(self->imagedata);
    sh = gdImageSY(self->imagedata);
    if(PyArg_ParseTuple(args, "O!|(ii)(ii)", &Imagetype, &dest, &dx, &dy,
      &sx, &sy))
    {
        dw = gdImageSX(dest->imagedata);
        dh = gdImageSY(dest->imagedata);
    }
    else if(PyErr_Clear(), !PyArg_ParseTuple(args, "O!|(ii)(ii)(ii)(ii)",
      &Imagetype, &dest, &dx, &dy, &sx, &sy, &dw, &dh, &sw, &sh))
        return NULL;
    gdImageCopyResampled(dest->imagedata, self->imagedata, X(dx), Y(dy),
      X(sx), Y(sy), W(dw), H(dh), W(sw), H(sh));

    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_copymergeto(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "copyMergeTo() requires gd 2.0 or later");
    return NULL;
#else
    imageobject *dest;
    int dx,dy,sx,sy,w,h,dw,dh,pct;

    dx=dy=sx=sy=0;
    w = gdImageSX(self->imagedata);
    h = gdImageSY(self->imagedata);
    pct = 100;
    if(!PyArg_ParseTuple(args, "O!|(ii)(ii)(ii)i", &Imagetype, &dest, &dx, &dy, &sx, &sy, &w, &h, &pct))
        return NULL;
    dw = gdImageSX(dest->imagedata);
    dh = gdImageSY(dest->imagedata);
    gdImageCopyMerge(dest->imagedata, self->imagedata, X(dx), Y(dy), X(sx), Y(sy), W(w), H(h), pct);

    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_copymergegrayto(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "copyMergeGrayTo() requires gd 2.0 or later");
    return NULL;
#else
    imageobject *dest;
    int dx,dy,sx,sy,w,h,dw,dh,pct;

    dx=dy=sx=sy=0;
    w = gdImageSX(self->imagedata);
    h = gdImageSY(self->imagedata);
    pct = 100;
    if(!PyArg_ParseTuple(args, "O!|(ii)(ii)(ii)i", &Imagetype, &dest, &dx, &dy, &sx, &sy, &w, &h, &pct))
        return NULL;
    dw = gdImageSX(dest->imagedata);
    dh = gdImageSY(dest->imagedata);
    gdImageCopyMergeGray(dest->imagedata, self->imagedata, X(dx), Y(dy), X(sx), Y(sy), W(w), H(h), pct);

    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_copypaletteto(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "copyPaletteTo() requires gd 2.0 or later");
    return NULL;
#else
    imageobject *dest;
    if(!PyArg_ParseTuple(args, "O!", &Imagetype, &dest))
        return NULL;

    gdImagePaletteCopy(dest->imagedata,  self->imagedata);
    Py_INCREF(Py_None);
    return Py_None;
#endif
}

static PyObject *image_compare(imageobject *self, PyObject *args)
{
#if GD2_VERS <= 1
 PyErr_SetString(PyExc_NotImplementedError,
   "compare() requires gd 2.0 or later");
    return NULL;
#else
    imageobject *dest;
    if(!PyArg_ParseTuple(args, "O!", &Imagetype, &dest))
        return NULL;

    return Py_BuildValue("i", gdImageCompare(dest->imagedata,
      self->imagedata));
#endif
}

static PyObject *image_interlace(imageobject *self, PyObject *args)
{
    int i;

    if(!PyArg_ParseTuple(args, "i", &i))
        return NULL;
    gdImageInterlace(self->imagedata, i);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_savealpha(imageobject *self, PyObject *args)
{
    int i;

    if(!PyArg_ParseTuple(args, "i", &i))
        return NULL;
    gdImageSaveAlpha(self->imagedata, i);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_origin(imageobject *self, PyObject *args)
{
    if(!PyArg_ParseTuple(args, "(ii)|ii",&self->origin_x,&self->origin_y,
        &self->multiplier_x,&self->multiplier_y))
        return NULL;

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *image_getorigin(imageobject *self)
{
    return Py_BuildValue("((ii)ii)",self->origin_x,self->origin_y,self->multiplier_x,self->multiplier_y);
}


static struct PyMethodDef image_methods[] = {

 {"writeGif",    (PyCFunction)image_writegif,    1,
    "writeGif(f)\n"
    "write the image to f as a GIF, where f is either an open file object or a\n"
    "file name."},

 {"writePng",    (PyCFunction)image_writepng,    1,
    "writePng(f)\n"
    "write the image to f as a PNG, where f is either an open file object or a\n"
    "file name."},
#ifdef HAVE_LIBJPEG
 {"writeJpeg",    (PyCFunction)image_writejpeg,    1,
    "writeJpeg(f,quality)\n"
    "write the image to f as a JPEG, where f is either an open file object or a\n"
    "file name.  quality is an optional integer value giving the JPEG quality from\n"
    "0 to 95%; leave out for default quality."},
#endif
 {"writeWbmp",    (PyCFunction)image_writewbmp,    1,
    "writeWbmp(f,fg)\n"
    "write the image to f as a WBMP, where f is either an open file object or a\n"
    "file name.  fg is a foreground color index which will be set in the output\n"
    "file; see the gd documentation for a further explanation of this value."},

 {"writeGd",    (PyCFunction)image_writegd,     1,
    "writeGd(f)\n"
    "write the image to f as a GD file, where f is either an open file object\n"
    "or a file name."},

 {"writeGd2",    (PyCFunction)image_writegd2,    1,
    "writeGd(f,chunksize,fmt)\n"
    "write the image to f as a GD2 file, where f is either an open file object\n"
    "or a file name.  see the gd documentation for an explanation of the chunksize\n"
    "and fmt arguments.  defaults will be used if not given."},

 {"setPixel",    (PyCFunction)image_setpixel,    1,
    "setPixel((x,y), color)\n"
    "set the pixel at (x,y) to color"},

 {"line",    (PyCFunction)image_line,    1,
    "line((x1,y1), (x2,y2), color)\n"
    "draw a line from (x1,y1) to (x2,y2) in color"},

 {"lines",    (PyCFunction)image_lines,    1,
    "line(seq, color)\n"
    "seq is a list of x,y tuples.  Draw a connected line in color"},


 {"polygon",    (PyCFunction)image_polygon,    1,
    "polygon(((x1,y1), (x2,y2), ..., (xn, yn)), color[, fillcolor])\n"
    "draw a polygon using the list or tuple of points (minimum 3) in color,\n"
    "optionally filled with fillcolor"},

 {"rectangle",    (PyCFunction)image_rectangle, 1,  "rectangle((x1,y1), (x2,y2), color[, fillcolor])\n"
    "draw a rectangle with upper corner (x1,y1), lower corner (x2,y2) in color,\n"
    "optionally filled with fillcolor"},

 {"filledPolygon",    (PyCFunction)image_filledpolygon, 1,
    "filledPolygon(((x1,y1), (x2,y2), ..., (xn, yn)), color)\n"
    "draw a filled polygon using the list or tuple of points (minimum 3) in color"},

 {"filledRectangle",    (PyCFunction)image_filledrectangle, 1,
    "filledRectangle((x1,y1), (x2,y2), color)\n"
    "draw a rectangle with upper corner (x1,y1), lower corner (x2,y2) in color"},

 {"arc",    (PyCFunction)image_arc,    1,
    "arc((x,y), (w,h), start, end, color)\n"
    "draw an ellipse centered at (x,y) with width w, height h from start\n"
    "degrees to end degrees in color."},

 {"filledArc", (PyCFunction)image_filledarc, 1,
  "filledArc((x,y), (w,h), start, end, color, style)\n"
  "draw an ellipse centered at (x,y) with width w, height h from start\n"
  "degrees to end degrees in color.  The style argument is a bitwise OR\n"
  "of the following (gdArc, gdChord, gdPie, gdNoFill, gdEdged)."},

 {"filledEllipse",    (PyCFunction)image_filledellipse,    1,
    "filledEllipse((x,y), (w,h), color)\n"
    "draw a filled ellipse centered at (x,y) with width w,\n"
    "height h in color."},

 {"fillToBorder",    (PyCFunction)image_filltoborder,    1,
    "fillToBorder((x,y), border, color)\n"
    "flood from point (x,y) to border color in color"},

 {"fill",    (PyCFunction)image_fill,    1,
    "fill((x,y), color)\n"
    "flood from point (x,y) in color for those pixels with the same color\n"
    "as the starting point"},

 {"setBrush",    (PyCFunction)image_setbrush,    1,
    "setBrush(image)\n"
    "set the drawing brush to <image> (use gdBrushed when drawing)"},

 {"setAntiAliased",    (PyCFunction)image_setantialiased, 1,  "setAntiAliased(color)\n"
    "Set the foreground color to be used when drawing antialiased lines.  Use the gdAntiAliased in place of the color when drawing"},

 {"setClip",    (PyCFunction)image_setclip, 1,  "setClip((x1,y1), (x2,y2))\n"
    "Set the image clipping rectangle"},


 {"setTile",    (PyCFunction)image_settile,    1,
    "setTile(image)\n"
    "set the fill tile to <image> (use gdTiled when filling)"},

 {"setStyle",    (PyCFunction)image_setstyle,    1,
    "setStyle(tuple)\n"
    "set the line bit-style to tuple of colors (use gdStyled when drawing)"},

 {"setThickness",    (PyCFunction)image_setthickness,    1,
    "setThickness(thickness)\n"
    "set the width of lines, polyons, and arcs."},

 {"alphaBlending", (PyCFunction)image_alphablending, 1,
    "alphaBlending(blending)\n"
    "toggles alpha channel blending for trucolor images."},

 {"alpha", (PyCFunction)image_alpha, 1,
    "alpha(color)\n"
    "returns the alpha component of the specified color."},

 {"blue", (PyCFunction)image_blue, 1,
    "blue(color)\n"
    "returns the blue component of the specified color."},

 {"green", (PyCFunction)image_green, 1,
    "green(color)\n"
    "returns the green component of the specified color."},

 {"red", (PyCFunction)image_red, 1,
    "red(color)\n"
    "returns the red component of the specified color."},

 {"getPixel",    (PyCFunction)image_getpixel,    1,
    "getPixel((x,y))\n"
    "color index of image at (x,y)"},

 {"boundsSafe",    (PyCFunction)image_boundssafe,    1,
    "boundsSafe((x,y))\n"
    "returns true if(x,y) is within image"},

 {"size",    (PyCFunction)image_size,    1,
    "size()\n"
    "return the 2-tuple size of image"},

 {"char",    (PyCFunction)image_char,1,
    "char(font, (x,y), c, color)\n"
    "draw string c at (x,y) using one of the pre-defined gdmodule fonts (gdFont*)"},

 {"charUp",    (PyCFunction)image_charup,1,
    "charUp(font, (x,y), c, color)\n"
    "vertically draw string c at (x,y) using one of the pre-defined gdmodule\n"
    "fonts (gdFont*)"},

 {"get_bounding_rect",    (PyCFunction)image_get_bounding_rect,1,
    "get_bounding_rect(font, ptsize, angle,(x,y), s)\n"
    "get bounding rect of string s using the truetype font"},

 {"string_ttf",    (PyCFunction)image_string_ttf,1,
    "string_ttf(font, ptsize, angle, (x,y), s, color)\n"
    "draw string s at (x,y) using the truetype font"},

 {"string_ft",    (PyCFunction)image_string_ft,1,
    "string_ft(font, ptsize, angle, (x,y), s, color)\n"
    "draw string s at (x,y) using the truetype font"},

 {"string",    (PyCFunction)image_string,1,
    "string(font, (x,y), s, color)\n"
    "draw string s at (x,y) using one of the pre-defined gdmodule fonts (gdFont*)"},

#ifdef Py_UNICODEOBJECT_H
 {"string16",    (PyCFunction)image_string16,1,
    "string(font, (x,y), s, color)\n"
    "draw unicode string s at (x,y) using one of the pre-defined gdmodule fonts (gdFont*)"},
#endif

 {"stringUp",    (PyCFunction)image_stringup,1,
    "stringUp(font, (x,y), s, color)\n"
    "vertically draw string s at (x,y) using one of the pre-defined gdmodule\n"
    "fonts (gdFont*)"},

#ifdef Py_UNICODEOBJECT_H
 {"stringUp16",    (PyCFunction)image_stringup16,1,
    "stringUp16(font, (x,y), s, color)\n"
    "vertically draw unicode string s at (x,y) using one of the pre-defined gdmodule\n"
    "fonts (gdFont*)"},
#endif

 {"colorAllocate",    (PyCFunction)image_colorallocate,1,
    "colorAllocate((r,g,b))\n"
    "allocate a color index to (r,g,b) (returns -1 if unable to)"},

 {"colorAllocateAlpha",    (PyCFunction)image_colorallocatealpha,1,
    "colorAllocateAlpha((r,g,b,a))\n"
    "allocate a color index to (r,g,b,a) (returns -1 if unable to)"},

 {"colorClosest",    (PyCFunction)image_colorclosest,    1,
    "colorClosest((r,g,b))\n"
    "return the color index closest to (r,g,b) (returns -1 if unable to)"},

 {"colorClosestAlpha",    (PyCFunction)image_colorclosestalpha,    1,
    "colorClosestAlpha((r,g,b,a))\n"
    "return the color index closest to (r,g,b,a) (returns -1 if unable to)"},

 {"colorClosestHWB",    (PyCFunction)image_colorclosestHWB,    1,
    "colorClosestHWB((r,g,b))\n"
    "return the color index closest in hue whiteness and blackness to (r,g,b) (returns -1 if unable to)"},

 {"colorExact",    (PyCFunction)image_colorexact,    1,
    "colorExact((r,g,b))\n"
    "return an exact color index match for (r,g,b) (returns -1 if unable to)"},

 {"colorResolve", (PyCFunction)image_colorresolve, 1,
  "colorResolve((r,g,b))\n"
  "returns the exact color (after possibly allocating) or the closes match."},

 {"colorResolveAlpha", (PyCFunction)image_colorresolvealpha, 1,
  "colorResolve((r,g,b,a))\n"
  "returns the exact color (after possibly allocating) or the closes match."},

 {"colorsTotal",    (PyCFunction)image_colorstotal,    1,
    "colorsTotal()\n"
    "returns the number of colors currently allocated"},

 {"colorComponents",    (PyCFunction)image_colorcomponents,    1,
    "colorComponents(color)\n"
    "returns a 3-tulple of the (r,g,b) components of color"},

 {"getClip",    (PyCFunction)image_getclip,    1,
    "getClip()\n"
    "returns (x1,y1), (x2,y2) of the image clip rectange"},

 {"getInterlaced",    (PyCFunction)image_getinterlaced,    1,
    "getInterlaced()\n"
    "returns true if the image is interlaced"},

 {"getTransparent",    (PyCFunction)image_gettransparent,    1,
    "getTransparent()\n"
    "returns transparent color index or -1"},

 {"colorDeallocate",    (PyCFunction)image_colordeallocate,    1,
    "colorDeallocate(color)\n"
    "deallocate color from the image palette"},

 {"colorTransparent",    (PyCFunction)image_colortransparent,    1,
    "colorTransparent(color)\n"
    "set the transparent color to color"},

 {"copyTo",    (PyCFunction)image_copyto,    1,
    "copyTo(image[, (dx,dy)[, (sx,sy)[, (w,h)]]])\n"
    "copy from (sx,sy), width sw and height sh to destination image (dx,dy)"},

 {"copyResizedTo",    (PyCFunction)image_copyresizedto,1,
    "copyResizedTocopyTo(image[, (dx,dy)[, (sx,sy)[, (dw,dh)[, (sw,sh)]]]])\n"
    "copy from (sx,sy), width sw and height sh to destination image (dx,dy), \n"
    "width dw and height dh"},

 {"copyResampledTo",    (PyCFunction)image_copyresampledto,1,
    "copyResampledTo(image[, (dx,dy)[, (sx,sy)[, (dw,dh)[, (sw,sh)]]]])\n"
    "copy from (sx,sy), width sw and height sh to destination image (dx,dy), \n"
    "width dw and height dh using smooth pixel interpolation."},

 {"copyMergeTo",    (PyCFunction)image_copymergeto,    1,
    "copyMergeTo(image[, (dx,dy)[, (sx,sy)[, (w,h)]], pct])\n"
    "copy from (sx,sy), width sw and height sh to destination image (dx,dy)\n"
    "If pct is 100 then this is equivalent to copyTo().  If pct is less\n"
    "than 100 the images are merged." },

 {"copyMergeGrayTo",    (PyCFunction)image_copymergegrayto,    1,
    "copyMergeGrayTo(image[, (dx,dy)[, (sx,sy)[, (w,h)]], pct])\n"
    "copy from (sx,sy), width sw and height sh to destination image (dx,dy)\n"
    "If pct is 100 then this is equivalent to copyTo().  If pct is less\n"
    "than 100 the images are merged.  The hue is preserved by first\n"
    "converting the destination pixels to gray." },

 {"copyPaletteTo", (PyCFunction)image_copypaletteto, 1,
    "copyPaletteTo(image)\n"
  "copy the palette from one image to another."},

 {"compare", (PyCFunction)image_compare, 1,
    "compare(image)\n"
  "compares this image with another.  Returns a bitmask whose values\n"
  "are composed of CMP_IMAGE, CMP_NUM_COLORS, CMP_COLOR, CMP_SIZE_X,\n"
  "CMP_SIZE_Y, CMP_TRANSPARENT, CMP_BACKGROUND, CMP_INTERLACE, CMP_TRUECOLOR"},

 {"interlace",    (PyCFunction)image_interlace,    1,
    "interlace()\n"
    "set the interlace bit"},

 {"origin",    (PyCFunction)image_origin,    1,
    "origin((x,y)[,xmult,ymult])\n"
    "set the origin of the image to (x,y) and multiply all x, y, width and\n"
    "height factors by xmult and ymult (typically either 1 or -1)"},

 {"getOrigin",    (PyCFunction)image_getorigin,    1,
    "getOrigin()\n"
    "returns the origin parameters ((x,y),xmult,ymult)"},

 {"saveAlpha",    (PyCFunction)image_savealpha,    1,
    "saveAlpha(n)\n"
    "if n = 1, alpha channel information will be saved"
    " (assuming the format supports it, i.e. PNG)"},

 {NULL,        NULL}        /* sentinel */
};


/*
** Table of file types understood to gd "constructors"
*/

static struct {
    char *ext;
    gdImagePtr (*func)(FILE*);

} ext_table[] = {

#ifdef HAVE_LIBGIF
    {"gif",  gdImageCreateFromGif},
#endif
#ifdef HAVE_LIBPNG
    {"png",  gdImageCreateFromPng},
#endif
#ifdef HAVE_LIBJPEG
    {"jpeg", gdImageCreateFromJpeg},
    {"jpg",  gdImageCreateFromJpeg},
    {"jfif", gdImageCreateFromJpeg},
#endif
    {"gd",   gdImageCreateFromGd},
    {"gd2",  gdImageCreateFromGd2},
    {"xbm",  gdImageCreateFromXbm}, /* this is internal to gd */

    {NULL, NULL}
};

static struct {
    char *ext;
    gdImagePtr (*func)(gdIOCtxPtr);

} ext_table_ctx[] = {

#ifdef HAVE_LIBPNG
    {"png",  gdImageCreateFromPngCtx},
#endif
#ifdef HAVE_LIBJPEG
    {"jpeg", gdImageCreateFromJpegCtx},
    {"jpg",  gdImageCreateFromJpegCtx},
    {"jfif", gdImageCreateFromJpegCtx},
#endif
    {"gd",   gdImageCreateFromGdCtx},
    {"gd2",  gdImageCreateFromGd2Ctx},

    {NULL, NULL}
};


/*
** Code to act as a gdIOCtx wrapper for a python file object
** (read-only; write support would not be difficult to add)
*/

struct PyFileIfaceObj_gdIOCtx {
    gdIOCtx ctx;
    PyObject *fileIfaceObj;
    PyObject *strObj;        // our reference to the data string we're currently
};

int PyFileIfaceObj_IOCtx_GetC(gdIOCtx *ctx)
{
    struct PyFileIfaceObj_gdIOCtx *pctx = (struct PyFileIfaceObj_gdIOCtx *)ctx;
    if (pctx->strObj) {
        Py_DECREF(pctx->strObj);
        pctx->strObj = NULL;
    }
    pctx->strObj = PyObject_CallMethod(pctx->fileIfaceObj, "read", "i", 1);
    if (!pctx->strObj || !PyString_Check(pctx->strObj)) {
        return EOF;
    }
    if (PyString_GET_SIZE(pctx->strObj) == 1) {
        return (int)(unsigned char)PyString_AS_STRING(pctx->strObj)[0];
    }
    return EOF;
}

int PyFileIfaceObj_IOCtx_GetBuf(gdIOCtx *ctx, void *data, int size)
{
    int err;
    char *value;
    struct PyFileIfaceObj_gdIOCtx *pctx = (struct PyFileIfaceObj_gdIOCtx *)ctx;
    if (pctx->strObj) {
        Py_DECREF(pctx->strObj);
        pctx->strObj = NULL;
    }
    pctx->strObj = PyObject_CallMethod(pctx->fileIfaceObj, "read", "i", size);
    if (!pctx->strObj) {
        return 0;
    }
    err = PyString_AsStringAndSize(pctx->strObj, &value, &size);
    if (err < 0) {
        /* this throws away the python exception since the gd library
         * won't pass it up properly.  gdmodule should create its own
         * since the "file" couldn't be read properly.  */
        PyErr_Clear();
        return 0;
    }
    memcpy(data, value, size);
    return size;
}

void PyFileIfaceObj_IOCtx_Free(gdIOCtx *ctx)
{
    struct PyFileIfaceObj_gdIOCtx *pctx = (struct PyFileIfaceObj_gdIOCtx *)ctx;
    if (pctx->strObj) {
        Py_DECREF(pctx->strObj);
        pctx->strObj = NULL;
    }
    if (pctx->fileIfaceObj) {
        Py_DECREF(pctx->fileIfaceObj);
        pctx->fileIfaceObj = NULL;
    }
    /* NOTE: we leave deallocation of the ctx structure itself to outside
     * code for memory allocation symmetry.  This function is safe to
     * call multiple times (gd should call it + we call it to be safe). */
}

struct PyFileIfaceObj_gdIOCtx * alloc_PyFileIfaceObj_IOCtx(PyObject *fileIfaceObj)
{
    struct PyFileIfaceObj_gdIOCtx *pctx;
    pctx = calloc(1, sizeof(struct PyFileIfaceObj_gdIOCtx));
    if (!pctx)
        return NULL;
    pctx->ctx.getC = PyFileIfaceObj_IOCtx_GetC;
    pctx->ctx.getBuf = PyFileIfaceObj_IOCtx_GetBuf;
    pctx->ctx.gd_free = PyFileIfaceObj_IOCtx_Free;
    Py_INCREF(fileIfaceObj);
    pctx->fileIfaceObj = fileIfaceObj;
    return pctx;
}

void free_PyFileIfaceObj_IOCtx(struct PyFileIfaceObj_gdIOCtx *pctx)
{
    if (!pctx)
        return;
    assert(pctx->ctx.gd_free != NULL);
    pctx->ctx.gd_free((gdIOCtxPtr)pctx);  // free its internal resources
    free(pctx);
}
/*
** Code to create the imageobject
*/

static imageobject *newimageobject(PyObject *args)
{
    imageobject *self, *srcimage;
    int xdim=0, ydim=0, i, trueColor=0;
    char *filename,*ext=0;
    FILE *fp;
    PyObject *readObj;

    if(!(self = PyObject_NEW(imageobject, &Imagetype)))
        return NULL;

    self->current_tile = self->current_brush = NULL;
    self->origin_x = self->origin_y = 0;
    self->multiplier_x = self->multiplier_y = 1;
    self->imagedata = NULL;

    if(PyArg_ParseTuple(args, "")) {
        PyErr_SetString(PyExc_ValueError,
            "image size or source filename required");
        Py_DECREF(self);
        return NULL;
    } else if(PyErr_Clear(),
        PyArg_ParseTuple(args, "O!|(ii)i",
            &Imagetype, &srcimage, &xdim, &ydim, &trueColor))
    {
        if(!xdim) xdim = gdImageSX(srcimage->imagedata);
        if(!ydim) ydim = gdImageSY(srcimage->imagedata);
#if GD2_VERS <= 1
        trueColor = 0;
#endif
        if (!trueColor){
            if(!(self->imagedata = gdImageCreate(xdim,ydim))) {
                Py_DECREF(self);
                return NULL;
            }
        }
#if GD2_VERS > 1
        else
            if(!(self->imagedata = gdImageCreateTrueColor(xdim,ydim))) {
                Py_DECREF(self);
                return NULL;
            }
#endif
        if(xdim == gdImageSX(srcimage->imagedata) &&
            ydim == gdImageSY(srcimage->imagedata))
            gdImageCopy(self->imagedata,srcimage->imagedata,0,0,0,0,xdim,ydim);
        else
            gdImageCopyResized(self->imagedata,srcimage->imagedata,
                0,0,0,0,xdim,ydim,gdImageSX(srcimage->imagedata),
                gdImageSY(srcimage->imagedata));
    } else if(PyErr_Clear(),
        PyArg_ParseTuple(args, "(ii)|i", &xdim, &ydim, &trueColor))
    {
        if(!xdim || !ydim) {
            PyErr_SetString(PyExc_ValueError, "dimensions cannot be 0");
            Py_DECREF(self);
            return NULL;
        }
#if GD2_VERS <= 1
        trueColor = 0;
#endif
        if (!trueColor){
            if(!(self->imagedata = gdImageCreate(xdim,ydim))) {
                Py_DECREF(self);
                return NULL;
            }
        }
#if GD2_VERS > 1
        else
            if(!(self->imagedata = gdImageCreateTrueColor(xdim,ydim))) {
                Py_DECREF(self);
                return NULL;
            }
#endif
    } else if(PyErr_Clear(), PyArg_ParseTuple(args, "s|s", &filename, &ext)) {

        if(!ext) {
            if(!(ext = strrchr(filename,'.'))) {
                PyErr_SetString(PyExc_IOError,
                    "need an extension to determine file type (.png|.jpeg|.jpg|.gd|.gd2|.xpm|.xbm)");
                Py_DECREF(self);
                return NULL;
            }
            ext++;
        }

        if(strcmp(ext, "xpm") == 0) {
#ifndef HAVE_LIBXPM
            PyErr_SetString(PyExc_NotImplementedError,
                    "XPM Support Not Available");
            Py_DECREF(self);
            return NULL;
#else
            /* gdImageCreateFromXpm takes a filename rather than a FILE * */

            if(!(self->imagedata = gdImageCreateFromXpm(filename))) {
                PyErr_SetString(PyExc_IOError, "corrupt or invalid image file");
                Py_DECREF(self);
                return(NULL);
            }

            return self;
#endif
        }

        if(!(fp = fopen(filename,"rb"))) {
            PyErr_SetFromErrno(PyExc_IOError);
            Py_DECREF(self);
            return(NULL);
        }

        for(i = 0; ext_table[i].ext != NULL; i++) {

            if(strcmp(ext, ext_table[i].ext) == 0) {

                if(!(self->imagedata = ext_table[i].func(fp))) {
                    fclose(fp);
                    PyErr_SetString(PyExc_IOError,
                        "corrupt or invalid image file (may be unsupported)");
                    Py_DECREF(self);
                    return(NULL);
                }

                fclose(fp);
                return self;
            }
        }

        /* if we fall out, we didn't find the file type */

        PyErr_SetString(PyExc_IOError,
            "unsupported file type (only .gif|.png|.jpeg|.jpg|.gd|.gd2|.xbm|.xpm accepted)");

        Py_DECREF(self);
        return(NULL);

    } else if(PyErr_Clear(),
        PyArg_ParseTuple(args, "Oz", &readObj, &ext))
    {
        struct PyFileIfaceObj_gdIOCtx *ourIOCtx = NULL;

        if (!PyObject_HasAttrString(readObj, "read")) {
            PyErr_SetString(ErrorObject, "non-Image objects must have a read() method");
            Py_DECREF(self);
            return NULL;
        }

        ourIOCtx = alloc_PyFileIfaceObj_IOCtx(readObj);
        if (!ourIOCtx) {
            PyErr_NoMemory();
            Py_DECREF(self);
            return(NULL);
        }

        for(i = 0; ext_table_ctx[i].ext != NULL; i++) {

            if(strcmp(ext, ext_table_ctx[i].ext) == 0) {

                if(!(self->imagedata = ext_table_ctx[i].func((gdIOCtxPtr)ourIOCtx))) {
                    free_PyFileIfaceObj_IOCtx(ourIOCtx);
                    PyErr_SetString(PyExc_IOError,
                        "corrupt or invalid image data (may be unsupported)");
                    Py_DECREF(self);
                    return(NULL);
                }

                free_PyFileIfaceObj_IOCtx(ourIOCtx);
                return self;
            }
        }

        /* if we fall out, we didn't find the file type */

        PyErr_SetString(PyExc_IOError,
            "unsupported file type (only png, jpeg, gd, & gd2 can be read from an object)");

        free_PyFileIfaceObj_IOCtx(ourIOCtx);

        Py_DECREF(self);
        return(NULL);
    } else {
        PyErr_SetString(PyExc_ValueError, "invalid argument list");
        Py_DECREF(self);
        return(NULL);
    }

    return self;
}


static void image_dealloc(imageobject *self)
{
    if(self->current_brush) {
        Py_DECREF(self->current_brush);
    }

    if(self->current_tile) {
        Py_DECREF(self->current_tile);
    }

    if(self->imagedata) gdImageDestroy(self->imagedata);

    PyObject_DEL(self);
}


static PyObject *image_getattr(PyObject *self, char *name)
{
    return Py_FindMethod(image_methods, self, name);
}


static PyObject *image_print(PyObject *self, FILE *fp, int flags)
{
    gdImagePtr im;

    im = ((imageobject *)self)->imagedata;
    fprintf(fp,"<%dx%d image object at 0x%lx>",gdImageSX(im),gdImageSY(im),
      (unsigned long)self);
    return 0;
}


static PyTypeObject Imagetype = {
    PyObject_HEAD_INIT(NULL)
    0,                                  /*ob_size*/
    "image",                            /*tp_name*/
    sizeof(imageobject),                /*tp_basicsize*/
    0,                                  /*tp_itemsize*/
    /* methods */
    (destructor)image_dealloc,          /*tp_dealloc*/
    (printfunc)image_print,             /*tp_print*/
    (getattrfunc)image_getattr,         /*tp_getattr*/
    0,                                  /*tp_setattr*/
    0,                                  /*tp_compare*/
    0,                                  /*tp_repr*/
    0,                                  /*tp_as_number*/
    0,                                  /*tp_as_sequence*/
    0,                                  /*tp_as_mapping*/
    0,                                  /*tp_hash*/
    0,                                  /*tp_call */
    0,                                  /*tp_str */
    0,0,0,0,                            /*tp_xxx1-4 */
};


static PyObject *gd_image(PyObject *self, PyObject *args)
{
    return (PyObject *)newimageobject(args);
}


static PyObject *gd_fontSSize(PyObject *self, PyObject *args)
{
    int font;
    char *str;
    int len;

    if(!PyArg_ParseTuple(args, "is", &font, &str))
        return NULL;

    if(font < 0 && font >= (sizeof(fonts)-1)) {
        PyErr_SetString(PyExc_ValueError, "Font value not valid");
        return NULL;
    }

    len = strlen(str);

    return Py_BuildValue("(ii)", len*(fonts[font].func()->w),
                 fonts[font].func()->h);
}


/*
** List of methods defined in the module
*/

static struct PyMethodDef gd_methods[] = {
    {"image", gd_image, 1,
        "image(image[,(w,h)] | file | file,type | (w,h))\n"
        "create GD image from file of type png, jpeg, gd, gd, gd2, xbm, or xpm.\n"
        "the existing image, optionally resized to width w and height h\n"
        "or blank with width w and height h"},
    {"fontstrsize", gd_fontSSize, 1,
        "fontstrsize(font, string)\n"
        "return a tuple containing the size in pixels of the given string in the\n"
        "given font"},
    {NULL,        NULL}        /* sentinel */
};


/* Initialization function for the module (*must* be called init_gd) */

void DLLEXPORT init_gd(void)
{
    PyObject *m, *d, *v;
    int i=0;

    /* Create the module and add the functions */
    m = Py_InitModule("_gd", gd_methods);

    /* Add some symbolic constants to the module */
    d = PyModule_GetDict(m);
    ErrorObject = PyString_FromString("gd.error");
    PyDict_SetItemString(d, "error", ErrorObject);

    /* add in the two font constants */
    while(fonts[i].name) {
        v = Py_BuildValue("i", i);
        PyDict_SetItemString(d, fonts[i++].name, v);
    }

    /* and the standard gd constants */

    v = Py_BuildValue("i", gdAntiAliased);
    PyDict_SetItemString(d, "gdAntiAliased", v);

    v = Py_BuildValue("i", gdBrushed);
    PyDict_SetItemString(d, "gdBrushed", v);

    v = Py_BuildValue("i", gdMaxColors);
    PyDict_SetItemString(d, "gdMaxColors", v);

    v = Py_BuildValue("i", gdMaxColors);
    PyDict_SetItemString(d, "gdMaxColors", v);

    v = Py_BuildValue("i", gdStyled);
    PyDict_SetItemString(d, "gdStyled", v);

    v = Py_BuildValue("i", gdStyledBrushed);
    PyDict_SetItemString(d, "gdStyledBrushed", v);

    v = Py_BuildValue("i", gdDashSize);
    PyDict_SetItemString(d, "gdDashSize", v);

    v = Py_BuildValue("i", gdTiled);
    PyDict_SetItemString(d, "gdTiled", v);

    v = Py_BuildValue("i", gdTransparent);
    PyDict_SetItemString(d, "gdTransparent", v);

#if GD2_VERS > 1
    v = Py_BuildValue("i", gdArc);
    PyDict_SetItemString(d, "gdArc", v);

    v = Py_BuildValue("i", gdChord);
    PyDict_SetItemString(d, "gdChord", v);

    v = Py_BuildValue("i", gdPie);
    PyDict_SetItemString(d, "gdPie", v);

    v = Py_BuildValue("i", gdNoFill);
    PyDict_SetItemString(d, "gdNoFill", v);

    v = Py_BuildValue("i", gdEdged);
    PyDict_SetItemString(d, "gdEdged", v);

    v = Py_BuildValue("i", GD_CMP_IMAGE);
    PyDict_SetItemString(d, "CMP_IMAGE", v);

    v = Py_BuildValue("i", GD_CMP_NUM_COLORS);
    PyDict_SetItemString(d, "CMP_NUM_COLORS", v);

    v = Py_BuildValue("i", GD_CMP_COLOR);
    PyDict_SetItemString(d, "CMP_COLOR", v);

    v = Py_BuildValue("i", GD_CMP_SIZE_X);
    PyDict_SetItemString(d, "CMP_SIZE_X", v);

    v = Py_BuildValue("i", GD_CMP_SIZE_Y);
    PyDict_SetItemString(d, "CMP_SIZE_Y", v);

    v = Py_BuildValue("i", GD_CMP_TRANSPARENT);
    PyDict_SetItemString(d, "CMP_TRANSPARENT", v);

    v = Py_BuildValue("i", GD_CMP_BACKGROUND);
    PyDict_SetItemString(d, "CMP_BACKGROUND", v);

    v = Py_BuildValue("i", GD_CMP_INTERLACE);
    PyDict_SetItemString(d, "CMP_INTERLACE", v);

    v = Py_BuildValue("i", GD_CMP_TRUECOLOR);
    PyDict_SetItemString(d, "CMP_TRUECOLOR", v);
#endif

    /* Check for errors */
    if(PyErr_Occurred())
        Py_FatalError("can't initialize module gd");
}


/* end of file. */
