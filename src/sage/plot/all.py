from graphics import show_default, Graphics
from plot import plot, graphics_array, list_plot, parametric_plot, polar_plot
from plot import plot_loglog, plot_semilogx, plot_semilogy
from plot import list_plot_loglog, list_plot_semilogx, list_plot_semilogy
from line import line, line2d
from arrow import arrow, arrow2d
from bar_chart import bar_chart
from bezier_path import bezier_path
from scatter_plot import scatter_plot
from disk import disk
from point import point, points, point2d
from matrix_plot import matrix_plot
from plot_field import plot_vector_field, plot_slope_field
from text import text
from polygon import polygon, polygon2d
from circle import circle
from ellipse import ellipse
from contour_plot import contour_plot, implicit_plot, region_plot
from density_plot import density_plot

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.complex_plot",["complex_plot"])

from arc import arc

from animate import animate

from plot3d.tachyon import Tachyon

from colors import Color, hue, rainbow, colors, colormaps

from step import plot_step_function

from hyperbolic_arc import hyperbolic_arc
from hyperbolic_triangle import hyperbolic_triangle
