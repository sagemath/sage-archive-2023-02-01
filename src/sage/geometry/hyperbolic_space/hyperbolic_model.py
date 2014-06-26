# -*- coding: utf-8 -*-
r"""
Hyperbolic Models

In this module, a hyperbolic model is a collection of data that allow
the user to implement new models of hyperbolic space with minimal effort.
The data include facts about the underlying set (such as whether the
model is bounded), facts about the metric (such as whether the model is
conformal), facts about the isometry group (such as whether it is a
linear or projective group), and more.  Generally speaking, any data
or method that pertains to the model itself -- rather than the points,
geodesics, or isometries of the model -- is implemented in this module.

Abstractly, a model of hyperbolic space is a connected, simply connected
manifold equipped with a complete Riemannian metric of constant curvature
`-1`.  This module records information sufficient to enable computations
in hyperbolic space without explicitly specifying the underlying set or
its Riemannian metric.  Although, see the
`SageManifolds <http://sagemanifolds.obspm.fr/>`_ project if 
you would like to take this approach.

This module implements the abstract base class for a model of hyperbolic
space of arbitrary dimension.  It also contains the implementations of
specific models of hyperbolic geometry.

AUTHORS:

- Greg Laun (2013): Initial version.

EXAMPLES:

We illustrate how the classes in this module encode data by comparing
the upper half plane (UHP), Poincaré disk (PD) and hyperboloid (HM)
models.  First we import::

    sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP as U
    sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelPD as P
    sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelHM as H

We note that the UHP and PD models are bounded while the HM model is
not::

   sage: U.bounded and P.bounded
   True
   sage: H.bounded
   False

The isometry groups of UHP and PD are projective, while that of HM is
linear:
    
    sage: U.isometry_group_is_projective
    True
    sage: H.isometry_group_is_projective
    False

The models are responsible for determining if the coordinates of points
and the matrix of linear maps are appropriate for constructing points
and isometries in hyperbolic space:
  
    sage: U.point_in_model(2 + I)
    True
    sage: U.point_in_model(2 - I)
    False
    sage: U.point_in_model(2)
    False
    sage: U.bdry_point_in_model(2)
    True

"""

#***********************************************************************
#
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_import import lazy_import
from sage.functions.other import imag, real
from sage.rings.all import CC, RR
from sage.rings.integer import Integer
from sage.symbolic.pynac import I
from sage.rings.infinity import infinity
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.matrix.all import matrix

lazy_import('sage.misc.misc', 'attrcall')
lazy_import('sage.modules.free_module_element', 'vector')
lazy_import('sage.functions.other','sqrt')

lazy_import('sage.geometry.hyperbolic_space.model_factory', 'ModelFactory')

#####################################################################
## Some helper functions

def SL2R_to_SO21(A):
    r"""
    Given a matrix in `SL(2, \RR)` return its irreducible representation in
    `O(2,1)`.

    Note that this is not the only homomorphism, but it is the only one
    that works in the context of the implemented 2D hyperbolic geometry
    models.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_model import SL2R_to_SO21
        sage: A = SL2R_to_SO21(identity_matrix(2))
        sage: J =  matrix([[1,0,0],[0,1,0],[0,0,-1]]) #Lorentzian Gram matrix
        sage: norm(A.transpose()*J*A - J) < 10**-4
        True
    """
    a,b,c,d = (A/A.det().sqrt()).list()
    B = matrix(3, [a*d + b*c, a*c - b*d, a*c + b*d, a*b - c*d,
                   Integer(1)/Integer(2)*a**2 - Integer(1)/Integer(2)*b**2 -
                   Integer(1)/Integer(2)*c**2 + Integer(1)/Integer(2)*d**2,
                   Integer(1)/Integer(2)*a**2 + Integer(1)/Integer(2)*b**2 -
                   Integer(1)/Integer(2)*c**2 - Integer(1)/Integer(2)*d**2,
                   a*b + c*d, Integer(1)/Integer(2)*a**2 -
                   Integer(1)/Integer(2)*b**2 + Integer(1)/Integer(2)*c**2 -
                   Integer(1)/Integer(2)*d**2, Integer(1)/Integer(2)*a**2 +
                   Integer(1)/Integer(2)*b**2 + Integer(1)/Integer(2)*c**2 +
                   Integer(1)/Integer(2)*d**2])
    B = B.apply_map(attrcall('real')) # Kill ~0 imaginary parts
    if A.det() > 0:
        return B
    else:
        # Orientation-reversing isometries swap the nappes of
        #  the lightcone.  This fixes that issue.
        return -B

def SO21_to_SL2R(M):
    r"""
    A homomorphism from `SO(2, 1)` to `SL(2, \RR)`.

    Note that this is not the only homomorphism, but it is the only one
    that works in the context of the implemented 2D hyperbolic geometry
    models.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_model import SO21_to_SL2R
        sage: (SO21_to_SL2R(identity_matrix(3)) - identity_matrix(2)).norm() < 10**-4
        True
    """
    ####################################################################
    # SL(2,R) is the double cover of SO (2,1)^+, so we need to choose  #
    # a lift.  I have formulas for the absolute values of each entry   #
    # a,b ,c,d of the lift matrix(2,[a,b,c,d]), but we need to choose  #
    # one entry to be positive.  I choose d for no particular reason,  #
    # unless d = 0, then we choose c > 0.  The basic strategy for this #
    # function is to find the linear map induced by the SO(2,1)        #
    # element on the Lie algebra sl(2, R).  This corresponds to the    #
    # Adjoint action by a matrix A or -A in SL(2,R).  To find which    #
    # matrix let X,Y,Z be a basis for sl(2,R) and look at the images   #
    # of X,Y,Z as well as the second and third standard basis vectors  #
    # for 2x2 matrices (these are traceless, so are in the Lie         #
    # algebra).  These corresponds to AXA^-1 etc and give formulas     #
    # for the entries of A.                                            #
    ####################################################################
    (m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9) = M.list()
    d = sqrt(Integer(1)/Integer(2)*m_5 - Integer(1)/Integer(2)*m_6 -
             Integer(1)/Integer(2)*m_8 + Integer(1)/Integer(2)*m_9)
    if M.det() > 0: #EPSILON?
        det_sign = 1
    elif M.det() < 0: #EPSILON?
        det_sign = -1
    if d > 0: #EPSILON?
        c = (-Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/d
        b = (-Integer(1)/Integer(2)*m_2 + Integer(1)/Integer(2)*m_3)/d
        ad = det_sign*1 + b*c # ad - bc = pm 1
        a = ad/d
    else: # d is 0, so we make c > 0
        c = sqrt(-Integer(1)/Integer(2)*m_5 - Integer(1)/Integer(2)*m_6 +
                  Integer(1)/Integer(2)*m_8 + Integer(1)/Integer(2)*m_9)
        d = (-Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/c
            #d = 0, so ad - bc = -bc = pm 1.
        b = - (det_sign*1)/c
        a = (Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/b
    A = matrix(2,[a,b,c,d])
    return A

def mobius_transform(A, z):
    r"""
    Given a matrix ``A`` in `GL(2, \CC)` and a point ``z`` in the complex
    plane return the mobius transformation action of ``A`` on ``z``.

    INPUT:

    - ``A`` -- a `2 \times 2` invertible matrix over the complex numbers
    - ``z`` -- a complex number or infinity

    OUTPUT:

    - a complex number or infinity

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform as mobius_transform
        sage: mobius_transform(matrix(2,[1,2,3,4]),2 + I)
        2/109*I + 43/109
        sage: y = var('y')
        sage: mobius_transform(matrix(2,[1,0,0,1]),x + I*y)
        x + I*y

    The matrix must be square and `2`x`2`::

        sage: mobius_transform(matrix([[3,1,2],[1,2,5]]),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring

        sage: mobius_transform(identity_matrix(3),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring

    The matrix can be symbolic or can be a matrix over the real
    or complex numbers, but must be invertible::

        sage: (a,b,c,d) = var('a,b,c,d');
        sage: mobius_transform(matrix(2,[a,b,c,d]),I)
        (I*a + b)/(I*c + d)

        sage: mobius_transform(matrix(2,[0,0,0,0]),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring
    """
    if A.ncols() == 2 and A.nrows() == 2 and A.det() != 0:
            (a,b,c,d) = A.list()
            if z == infinity:
                if c == 0:
                    return infinity
                return a/c
            if a*d - b*c < 0:
                w = z.conjugate() # Reverses orientation
            else:
                w = z
            if c*z + d == 0:
                return infinity
            else:
                return (a*w + b)/(c*w + d)
    else:
        raise TypeError("A must be an invertible 2x2 matrix over the"
                        " complex numbers or a symbolic ring")

def PD_preserve_orientation(A):
    r"""
    For a PD isometry, determine if it preserves orientation.
    This test is more more involved than just checking the sign
    of the determinant, and it is used a few times in this file.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_model import PD_preserve_orientation as orient
        sage: orient(matrix(2, [-I, 0, 0, I]))
        True
        sage: orient(matrix(2, [0, I, I, 0]))
        False
    """
    return bool(A[1][0] == A[0][1].conjugate() and A[1][1] == A[0][0].conjugate()
                and abs(A[0][0]) - abs(A[0][1]) != 0)


#####################################################################
## The actual classes

class HyperbolicModel(UniqueRepresentation):
    r"""
    Abstract base class for Hyperbolic Models.
    """
    name = "Abstract Hyperbolic Model"
    short_name = "Abstract"
    bounded = False
    conformal = False
    dimension = 0
    isometry_group = None
    isometry_group_is_projective = False
    pt_conversion_dict = {}
    isom_conversion_dict = {}

    @classmethod
    def point_in_model(cls, p): #Abstract
        r"""
        Return ``True`` if the point is in the given model and ``False``
        otherwise.

        INPUT:

        - any object that can converted into a complex number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: UHP.point_in_model(I)
            True
            sage: UHP.point_in_model(-I)
            False
        """
        return True

    @classmethod
    def point_test(cls, p): #Abstract
        r"""
        Test whether a point is in the model.  If the point is in the
        model, do nothing.  Otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.point_test(2 + I)
            sage: HyperbolicModelUHP.point_test(2 - I)
            Traceback (most recent call last):
            ...
            ValueError: -I + 2 is not a valid point in the UHP model
        """
        if not (cls.point_in_model(p) or cls.bdry_point_in_model(p)):
            error_string = "{0} is not a valid point in the {1} model"
            raise ValueError(error_string.format(p, cls.short_name))

    @classmethod
    def bdry_point_in_model(cls, p): #Abstract
        r"""
        Return ``True`` if the point is on the ideal boundary of hyperbolic
        space and ``False`` otherwise.

        INPUT:

        - any object that can converted into a complex number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: UHP.bdry_point_in_model(I)
            False
        """
        return True

    @classmethod
    def bdry_point_test(cls, p): #Abstract
        r"""
        Test whether a point is in the model.  If the point is in the
        model, do nothing; otherwise raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.bdry_point_test(2)
            sage: HyperbolicModelUHP.bdry_point_test(1 + I)
            Traceback (most recent call last):
            ...
            ValueError: I + 1 is not a valid boundary point in the UHP model
        """
        if not cls.bounded or not cls.bdry_point_in_model(p):
            error_string = "{0} is not a valid boundary point in the {1} model"
            raise ValueError(error_string.format(p, cls.short_name))

    @classmethod
    def isometry_in_model(cls, A): #Abstract
        r"""
        Return ``True`` if the input matrix represents an isometry of the
        given model and ``False`` otherwise.

        INPUT:

        - a matrix that represents an isometry in the appropriate model

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: UHP.isometry_in_model(identity_matrix(2))
            True

            sage: UHP.isometry_in_model(identity_matrix(3))
            False
        """
        return True

    @classmethod
    def isometry_act_on_point(cls, A, p): #Abtsract
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: I2 = identity_matrix(2)
            sage: p = UHP.random_point().coordinates()
            sage: bool(norm(HyperbolicModelUHP.isometry_act_on_point(I2, p) - p) < 10**-9)
            True
        """
        return A * vector(p)

    @classmethod
    def isometry_test(cls, A): #Abstract
        r"""
        Test whether an isometry is in the model.

        If the isometry is in the model, do nothing. Otherwise, raise
        a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.isometry_test(identity_matrix(2))
            sage: HyperbolicModelUHP.isometry_test(matrix(2, [I,1,2,1]))
            Traceback (most recent call last):
            ...
            ValueError:
            [I 1]
            [2 1] is not a valid isometry in the UHP model.
        """
        if not cls.isometry_in_model(A):
            error_string = "\n{0} is not a valid isometry in the {1} model."
            raise ValueError(error_string.format(A, cls.short_name))

    @classmethod
    def point_to_model(cls, coordinates, model_name): #Abstract
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: UHP.point_to_model(I, 'UHP')
            I
            sage: UHP.point_to_model(I, 'PD')
            0
            sage: UHP.point_to_model(3 + I, 'KM')
            (6/11, 9/11)
            sage: UHP.point_to_model(3 + I, 'HM')
            (3, 9/2, 11/2)
            sage: UHP.point_to_model(I, 'PD')
            0
            sage: PD.point_to_model(0, 'UHP')
            I
            sage: UHP.point_to_model(I, 'UHP')
            I
            sage: KM.point_to_model((0, 0), 'UHP')
            I
            sage: KM.point_to_model((0, 0), 'HM')
            (0, 0, 1)
            sage: HM.point_to_model(vector((0,0,1)), 'UHP')
            I
            sage: HM.point_to_model(vector((0,0,1)), 'KM')
            (0, 0)

        It is an error to try to convert a boundary point to a model
        that doesn't support boundary points::

            sage: UHP.point_to_model(infinity, 'HM')
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented for the HM model
        """
        cls.point_test(coordinates)
        model = ModelFactory.find_model(model_name)
        if (not model.bounded) and cls.bdry_point_in_model(coordinates):
            raise NotImplementedError("boundary points are not implemented for"
                                      " the {0} model".format(model_name))
        return cls.pt_conversion_dict[model_name](coordinates)

    @classmethod
    def isometry_to_model(cls, A, model_name): #Abstract
        r"""
        Convert ``A`` from the current model to the model specified in
        ``model_name``.

        INPUT:

        - ``A`` -- a matrix in the current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: A = matrix(2,[I, 0, 0, -I])
            sage: PD.isometry_to_model(A, 'UHP')
            [ 0  1]
            [-1  0]

            sage: PD.isometry_to_model(A, 'HM')
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: PD.isometry_to_model(A, 'KM')
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: B = diagonal_matrix([-1, -1, 1])
            sage: HM.isometry_to_model(B, 'UHP')
            [ 0 -1]
            [ 1  0]

            sage: HM.isometry_to_model(B, 'PD')
            [-I  0]
            [ 0  I]

            sage: HM.isometry_to_model(B, 'KM')
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]
        """
        cls.isometry_test(A)
        return cls.isom_conversion_dict[model_name](A)


class HyperbolicModelUHP(HyperbolicModel, UniqueRepresentation):
    r"""
    Upper Half Plane model.
    """
    name = "Upper Half Plane Model"
    short_name = "UHP"
    bounded = True
    conformal = True
    dimension = 2
    isometry_group = "PSL(2, \\Bold{R})"
    isometry_group_is_projective = True
    pt_conversion_dict = {
        'UHP' : lambda p :  p,
        'PD' : lambda p :  (p - I)/(Integer(1) - I*p),
        'HM' : lambda p :  vector((real(p)/imag(p),
                                   (real(p)**2 + imag(p)**2 - 1)/(2*imag(p)),
                                   (real(p)**2 + imag(p)**2 + 1)/(2*imag(p)))),
        'KM' : lambda p :   ((2*real(p))/(real(p)**2 + imag(p)**2 + 1),
                             (real(p)**2 + imag(p)**2 - 1)/(real(p)**2 +
                                                            imag(p)**2 + 1))
        }
    isom_conversion_dict = {
        'UHP': lambda A : A,
        'PD' : lambda A : matrix(2,[1,-I,-I,1]) * A * matrix(2,[1,I,I,1])/Integer(2),
        'HM' : SL2R_to_SO21,
        'KM' : SL2R_to_SO21
        }

    @classmethod
    def point_in_model(cls, p): #UHP
        r"""
        Check whether a complex number lies in the open upper half
        plane. In the UHP.model_name_name, this is hyperbolic space.

        EXAMPLES::

            sage: UHP.point_in_model(1 + I)
            True
            sage: UHP.point_in_model(infinity)
            False
            sage: UHP.point_in_model(CC(infinity))
            False
            sage: UHP.point_in_model(RR(infinity))
            False
            sage: UHP.point_in_model(1)
            False
            sage: UHP.point_in_model(12)
            False
            sage: UHP.point_in_model(1 - I)
            False
            sage: UHP.point_in_model(-2*I)
            False
        """
        return bool(imag(CC(p)) > 0)

    @classmethod
    def bdry_point_in_model(cls, p): #UHP
        r"""
        Check whether a complex number is a real number or ``\infty``.
        In the ``UHP.model_name_name``, this is the ideal boundary of
        hyperbolic space.

        EXAMPLES::

            sage: UHP.bdry_point_in_model(1 + I)
            False
            sage: UHP.bdry_point_in_model(infinity)
            True
            sage: UHP.bdry_point_in_model(CC(infinity))
            True
            sage: UHP.bdry_point_in_model(RR(infinity))
            True
            sage: UHP.bdry_point_in_model(1)
            True
            sage: UHP.bdry_point_in_model(12)
            True
            sage: UHP.bdry_point_in_model(1 - I)
            False
            sage: UHP.bdry_point_in_model(-2*I)
            False
        """
        im = abs(imag(CC(p)).n())
        return bool( (im < EPSILON) or (p == infinity) )

    @classmethod #UHP
    def isometry_act_on_point(cls, A, p): #UHP
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: I2 = identity_matrix(2)
            sage: p = UHP.random_point().coordinates()
            sage: bool(norm(HyperbolicModelUHP.isometry_act_on_point(I2, p) - p) < 10**-9)
            True
        """
        return mobius_transform(A, p)

    @classmethod
    def isometry_in_model(cls, A): #UHP
        r"""
        Check that ``A`` acts as an isometry on the upper half plane.
        That is, ``A`` must be an invertible `2 \times 2` matrix with real
        entries.

        EXAMPLES::

            sage: A = matrix(2,[1,2,3,4])
            sage: UHP.isometry_in_model(A)
            True
            sage: B = matrix(2,[I,2,4,1])
            sage: UHP.isometry_in_model(B)
            False
        """
        return bool(A.ncols() == 2 and A.nrows() == 2 and
            sum([k in RR for k in A.list()]) == 4 and
            abs(A.det()) > -EPSILON)

    @classmethod
    def point_to_model(cls, coordinates, model_name): #UHP
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: UHP.point_to_model(I, 'UHP')
            I
            sage: UHP.point_to_model(I, 'PD')
            0
            sage: UHP.point_to_model(3 + I, 'KM')
            (6/11, 9/11)
            sage: UHP.point_to_model(3 + I, 'HM')
            (3, 9/2, 11/2)

        It is an error to try to convert a boundary point to a model
        that doesn't support boundary points::

            sage: UHP.point_to_model(infinity, 'HM')
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented for the HM model
        """
        p = coordinates
        if (cls.bdry_point_in_model(p) and not
                ModelFactory.find_model(model_name).bounded):
            raise NotImplementedError("boundary points are not implemented for"
                                      " the {0} model".format(model_name))
        if p == infinity:
            return {
                'UHP' : p,
                'PD' : I,
                'KM' : (0, 1)
                }[model_name]
        return cls.pt_conversion_dict[model_name](coordinates)

    @classmethod # UHP
    def isometry_to_model(cls, A, model_name): # UHP
        r"""
        Convert ``A`` from the current model to the model specified in
        ``model_name``.

        INPUT:

        - ``A`` -- a matrix in the current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.isometry_to_model(matrix(2,[0, 1, 1, 0]),'PD')
            [0 I]
            [I 0]
        """
        cls.isometry_test(A)
        if A.det() < 0 and model_name == 'PD':
            return cls.isom_conversion_dict[model_name](I * A)
        return cls.isom_conversion_dict[model_name](A)


class HyperbolicModelPD(HyperbolicModel, UniqueRepresentation):
    r"""
    Poincaré Disk Model.
    """
    name = "Poincare Disk Model" # u"Poincaré Disk Model"
    short_name = "PD"
    bounded = True
    conformal = True
    dimension = 2
    isometry_group = "PU (1, 1)"
    isometry_group_is_projective = True
    pt_conversion_dict =  {
            'PD': lambda p :  p,
            'UHP': lambda p :  (p + I)/(Integer(1) + I*p),
            'HM' : lambda p :  vector((
                2*real(p)/(1 - real(p)**2 - imag(p)**2),
                2*imag(p)/(1 - real(p)**2 - imag(p)**2),
                (real(p)**2 + imag(p)**2 + 1)/(1 -real(p)**2 - imag(p)**2)
                )),
            'KM' : lambda p :  (
                2*real(p)/(Integer(1) + real(p)**2 +imag(p)**2),
                2*imag(p)/(Integer(1) + real(p)**2 + imag(p)**2)
                )
            }
    isom_conversion_dict = {
            'PD' : lambda A : A,
            'UHP': lambda A :  (matrix(2,[1,I,I,1])*A*
                matrix(2,[1,-I,-I,1])/Integer(2)),
            'KM' : lambda A : SL2R_to_SO21( matrix(2,[1,I,I,1]) * A *
                                             matrix(2,[1,-I,-I,1])/Integer(2)),
            'HM' : lambda A : SL2R_to_SO21( matrix(2,[1,I,I,1]) * A *
                                             matrix(2,[1,-I,-I,1])/Integer(2))
            }

    @classmethod
    def point_in_model(cls, p): #PD
        r"""
        Check whether a complex number lies in the open unit disk.

        EXAMPLES::

            sage: PD.point_in_model(1.00)
            False

            sage: PD.point_in_model(1/2 + I/2)
            True

            sage: PD.point_in_model(1 + .2*I)
            False
        """
        return bool(abs(CC(p)) < 1)


    @classmethod
    def bdry_point_in_model(cls, p): #PD
        r"""
        Check whether a complex number lies in the open unit disk.

        EXAMPLES::

            sage: PD.bdry_point_in_model(1.00)
            True

            sage: PD.bdry_point_in_model(1/2 + I/2)
            False

            sage: PD.bdry_point_in_model(1 + .2*I)
            False
        """
        return  bool(abs(abs(CC(p))- 1) < EPSILON)


    @classmethod
    def isometry_act_on_point(cls, A, p): #PD
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelPD
            sage: I2 = identity_matrix(2)
            sage: q = PD.random_point().coordinates()
            sage: bool(norm(HyperbolicModelPD.isometry_act_on_point(I2, q) - q) < 10**-9)
            True
        """
        _image = mobius_transform(A, p)
        if not PD_preserve_orientation(A):
            return mobius_transform(I*matrix(2,[0,1,1,0]), _image)
        return _image


    @classmethod
    def isometry_in_model(cls, A): #PD
        r"""
        Check if the given matrix ``A`` is in the group `U(1,1)`.

        EXAMPLES::

            sage: z = [CC.random_element() for k in range(2)]; z.sort(key=abs)
            sage: A = matrix(2,[z[1], z[0],z[0].conjugate(),z[1].conjugate()])
            sage: PD.isometry_in_model(A)
            True
        """
        # alpha = A[0][0]
        # beta = A[0][1]
        # Orientation preserving and reversing
        return PD_preserve_orientation(A) or PD_preserve_orientation(I*A)

    @classmethod
    def point_to_model(cls, coordinates, model_name): #PD
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: PD.point_to_model(0, 'UHP')
            I

            sage: PD.point_to_model(I, 'UHP')
            +Infinity

            sage: PD.point_to_model(-I, 'UHP')
            0
        """
        if model_name == 'UHP' and coordinates == I:
            return infinity
        return super(HyperbolicModelPD, cls).point_to_model(coordinates,
                                                            model_name)

    @classmethod # PD
    def isometry_to_model(cls, A, model_name): # PD
        r"""
        Convert ``A`` from the current model to the model specified in
        ``model_name``.

        INPUT:

        - ``A`` -- a matrix in the current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES:

        We check that orientation-reversing isometries behave as they
        should::

            sage: PD.isometry_to_model(matrix(2,[0,I,I,0]),'UHP')
            [ 0 -1]
            [-1  0]
        """
        cls.isometry_test(A)
        # Check for orientation-reversing isometries.
        if (not PD_preserve_orientation(A) and model_name == 'UHP'):
            return cls.isom_conversion_dict[model_name](I*A)
        return cls.isom_conversion_dict[model_name](A)

class HyperbolicModelKM(HyperbolicModel, UniqueRepresentation):
    r"""
    Klein Model.
    """
    name = "Klein Disk Model"
    short_name = "KM"
    bounded = True
    conformal = False
    dimension = 2
    isometry_group_is_projective = True
    isometry_group = "PSO(2, 1)"
    pt_conversion_dict = {
            'UHP' : lambda p:  -p[0]/(p[1] - 1) +\
                I*(-(sqrt(-p[0]**2 -p[1]**2 + 1) - p[0]**2 - p[1]**2 +
                     1)/((p[1] - 1)*sqrt(-p[0]**2 - p[1]**2 + 1) + p[1] - 1)),
            'PD' : lambda p :  (p[0]/(1 + (1 - p[0]**2 - p[1]**2).sqrt()) +  \
                               I*p[1]/(1 + (1 - p[0]**2 - p[1]**2).sqrt())),
            'KM' : lambda p :  p,
            'HM' : lambda p :  vector((2*p[0],2*p[1], 1 + p[0]**2 +
                                       p[1]**2))/(1 - p[0]**2 - p[1]**2)
            }
    isom_conversion_dict = {
            'UHP' : SO21_to_SL2R,
            'PD' : lambda A :   matrix(2,[1,-I,-I,1]) * SO21_to_SL2R(A) *\
                matrix(2,[1,I,I,1])/Integer(2),
            'KM' : lambda A :  A,
            'HM' : lambda A :  A
            }

    @classmethod
    def point_in_model(cls, p): #KM
        r"""
        Check whether a point lies in the open unit disk.

        EXAMPLES::

            sage: KM.point_in_model((1,0))
            False

            sage: KM.point_in_model((1/2 , 1/2))
            True

            sage: KM.point_in_model((1 , .2))
            False
        """
        return len(p) == 2 and bool(p[0]**2 + p[1]**2 < 1)

    @classmethod
    def bdry_point_in_model(cls, p): #KM
        r"""
        Check whether a point lies in the unit circle, which corresponds
        to the ideal boundary of the hyperbolic plane in the Klein model.

        EXAMPLES::

            sage: KM.bdry_point_in_model((1,0))
            True

            sage: KM.bdry_point_in_model((1/2 , 1/2))
            False

            sage: KM.bdry_point_in_model((1 , .2))
            False
        """
        return len(p) == 2 and bool(abs(p[0]**2 + p[1]**2 - 1) < EPSILON)


    @classmethod #KM
    def isometry_act_on_point(cls, A, p): #KM
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelKM
            sage: I3 = identity_matrix(3)
            sage: v = vector(KM.random_point().coordinates())
            sage: bool(norm(HyperbolicModelKM.isometry_act_on_point(I3, v) - v) < 10**-9)
            True
        """
        v = A*vector((list(p) + [1]))
        if v[2] == 0:
            return infinity
        return v[0:2]/v[2]

    @classmethod
    def isometry_in_model(cls, A): #KM
        r"""
        Check if the given matrix ``A`` is in the group `SO(2,1)`.

        EXAMPLES::

            sage: A = matrix(3,[1, 0, 0, 0, 17/8, 15/8, 0, 15/8, 17/8])
            sage: KM.isometry_in_model(A)
            True
        """
        from sage.geometry.hyperbolic_space.hyperbolic_constants import LORENTZ_GRAM
        return bool((A*LORENTZ_GRAM*A.transpose() - LORENTZ_GRAM).norm()**2 <
                EPSILON)

    @classmethod
    def point_to_model(cls, coordinates, model_name): #KM
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: KM.point_to_model((0, 0), 'UHP')
            I

            sage: KM.point_to_model((0, 0), 'HM')
            (0, 0, 1)

            sage: KM.point_to_model((0,1), 'UHP')
            +Infinity
        """
        if model_name == 'UHP' and tuple(coordinates) == (0,1):
            return infinity
        return super(HyperbolicModelKM, cls).point_to_model(coordinates,
                                                            model_name)


class HyperbolicModelHM(HyperbolicModel, UniqueRepresentation):
    r"""
    Hyperboloid Model.
    """
    name = "Hyperboloid Model"
    short_name = "HM"
    bounded = False
    conformal = True
    dimension = 2
    isometry_group = "SO(2, 1)"
    pt_conversion_dict = {
            'UHP' : lambda p : -((p[0]*p[2] + p[0]) +
                                 I*(p[2] +1))/((p[1] - 1)*p[2] - p[0]**2 -
                                               p[1]**2 + p[1] - 1),
            'PD' : lambda p : p[0]/(1 + p[2]) + I* (p[1]/(1 + p[2])),
            'KM' : lambda p : (p[0]/(1 + p[2]), p[1]/(1 + p[2])),
            'HM' : lambda p : p
            }
    isom_conversion_dict =  {
            'UHP' : SO21_to_SL2R,
            'PD' : lambda A : matrix(2,[1,-I,-I,1]) * SO21_to_SL2R(A) *\
                matrix(2,[1,I,I,1])/Integer(2),
            'KM' : lambda A : A,
            'HM' : lambda A : A
            }

    @classmethod
    def point_in_model(cls, p): #HM
        r"""
        Check whether a complex number lies in the hyperboloid.

        EXAMPLES::

            sage: HM.point_in_model((0,0,1))
            True

            sage: HM.point_in_model((1,0,sqrt(2)))
            True

            sage: HM.point_in_model((1,2,1))
            False
        """
        return len(p) == 3 and bool(p[0]**2 + p[1]**2 - p[2]**2 + 1 < EPSILON)

    @classmethod
    def bdry_point_in_model(cls, p):  #HM
        r"""
        Return ``False`` since the Hyperboloid model has no boundary points.

        EXAMPLES::

            sage: HM.bdry_point_in_model((0,0,1))
            False

            sage: HM.bdry_point_in_model((1,0,sqrt(2)))
            False

            sage: HM.bdry_point_in_model((1,2,1))
            False
        """
        return False


    @classmethod
    def isometry_in_model(cls, A):  #HM
        r"""
        Test that the matrix ``A`` is in the group `SO(2,1)^+`.

        EXAMPLES::

           sage: A = diagonal_matrix([1,1,-1])
           sage: HM.isometry_in_model(A)
           True
        """
        from sage.geometry.hyperbolic_space.hyperbolic_constants import LORENTZ_GRAM
        return bool((A*LORENTZ_GRAM*A.transpose() - LORENTZ_GRAM).norm()**2 < EPSILON)

