# -*- python -*-
# encoding: utf-8
"""
plot.py

Created by Michael Brickenstein on 2008-10-17.
Copyright (c) 2008 The PolyBoRi Team.
"""

import sys
import os
from polybori.PyPolyBoRi import Ring, Polynomial, BooleSet
from subprocess import Popen, PIPE


graph_template="""
digraph polynomial{
graph [
ordering="out"
#if highlight_monomial
, label = "${display_monomial(highlight_monomial)}"
#end
, fontsize=${fontsize}
];

#for n in nodes
${identifier(n)}[label="${label(n)}", shape="${shape(n)}"];
#end

#for n in non_constant_nodes
${identifier(n)} -> ${identifier(n.else_branch())} [style="dashed", color="${color_else}", arrowhead="vee", penwidth="${penwidth_else(n)}"];
${identifier(n)} -> ${identifier(n.then_branch())} [color="${color_then}", arrowhead="vee", penwidth="${penwidth_then(n)}"];
#end
}
"""

graph_template_jinja="""
digraph polynomial{
{% if landscape %}
rankdir=LR;
{% endif %}
graph [ ordering="out"
{% if highlight_monomial %}
, label = "{{display_monomial(highlight_monomial)}}"
{% endif %}
, fontsize={{fontsize}}
];

{% for n in nodes %}
{{identifier(n)}}[label="{{label(n)}}", shape="{{shape(n)}}"];
{% endfor %}

{% for n in non_constant_nodes %}
{{identifier(n)}} -> {{identifier(n.else_branch())}} [style="dashed", color="{{color_else}}", arrowhead="vee", penwidth="{{penwidth_else(n)}}"];
{{identifier(n)}} -> {{identifier(n.then_branch())}} [color="{{color_then}}", arrowhead="vee", penwidth="{{penwidth_then(n)}}"];
{% endfor %}
}
"""


ELSE="else"
THEN="then"

def render_genshi(data_dict):
    from genshi.template import TextTemplate
    tmpl = TextTemplate(graph_template)
    stream = tmpl.generate(**data_dict)
    return str(stream)

def render_jinja(data_dict):
    try:
        from jinja2 import Environment
    except:
        from jinja import Environment
    env = Environment()
    tmpl = env.from_string(graph_template_jinja)
    return tmpl.render(**data_dict)

def monomial_path_in_zdd(mon, graph):
    res=[]
    if not mon in BooleSet(graph):
        raise ValueError
    graph_nav=BooleSet(graph).navigation()
    mon_nav=BooleSet(mon).navigation()
    while not mon_nav.constant():
        while graph_nav.value()<mon_nav.value():
            res.append((graph_nav, ELSE))
            graph_nav=graph_nav.else_branch()
        assert mon_nav.value()==graph_nav.value()
        res.append((graph_nav,THEN))
        mon_nav=mon_nav.then_branch()
        graph_nav=graph_nav.then_branch()
    while not graph_nav.constant():
        res.append((graph_nav, ELSE))
        graph_nav=graph_nav.else_branch()
    return dict(res)
def plot(p, filename, colored=True,format="png",
    highlight_monomial=None, fontsize=14,
    template_engine='jinja', landscape=False
    ):
    """plots ZDD structure to <filename> in format <format>

    EXAMPLES:

    >>> r=Ring(1000)
    >>> x = Variable = VariableFactory(r)
    >>> plot(x(1)+x(0),"/dev/null", colored=True)
    >>> plot(x(1)+x(0),"/dev/null", colored=False)
    """
    THICK_PEN=5
    highlight_path=dict()
    if highlight_monomial:
        highlight_path=monomial_path_in_zdd(highlight_monomial, p)
    def display_monomial(m):
        return unicode(m).replace("*",u"⋅")
    def penwidth_else(n):
        if n in highlight_path and highlight_path[n]==ELSE:
            return THICK_PEN
        return 1

    def penwidth_then(n):
        if n in highlight_path and highlight_path[n]==THEN:
            return THICK_PEN
        return 1
    if not colored:
        color_then="black"
        color_else="black"
    else:
        color_then="red"
        color_else="blue"


    def find_navs(nav):
        if not nav in nodes:
            nodes.add(nav)
            if not nav.constant():
                find_navs(nav.then_branch())
                find_navs(nav.else_branch())
    p=Polynomial(p)
    nodes=set()
    nav=p.navigation()
    find_navs(nav)
    non_constant_nodes=[n for n in nodes if not n.constant()]
    node_to_int=dict([(n,i) for (i,n) in enumerate(nodes)])


    r=p.ring()
    def identifier(n):
        return "n"+str(node_to_int[n])
    def label(n):
        if n.constant():
            if n.terminal_one():
                return "1"
            else:
                return "0"
        else:
            return str(r.variable(n.value()))
    def shape(n):
        if n.constant():
            return "box"
        else:
            return "ellipse"
    renderers=dict(genshi=render_genshi, jinja=render_jinja)

    dot_input=renderers[template_engine](locals())
    if isinstance(dot_input, unicode):
        dot_input=dot_input.encode('utf-8')
    process = Popen(["dot", "-T"+format, "-o",filename], stdin=PIPE, stdout=PIPE)

    process.stdin.write(dot_input)
    process.stdin.close()
    process.wait()


def main():
    r=Ring(1000)
    x = Variable = VariableFactory(r)
    from os import system
    from polybori.specialsets import all_monomials_of_degree_d, power_set
    full_set=list(power_set([Variable(i) for i in xrange(15)]))
    from random import Random
    generator=Random(123)
    random_set=sum(generator.sample(full_set,30))
    full_polynomial=list(all_monomials_of_degree_d(3, [Variable(i) for i in xrange(30)]))
    random_poly=sum(generator.sample(full_polynomial,30))
    polynomials=[
    x(1)*x(2)+x(3),
    (x(1)+1)*(x(2)+x(3)),
    (x(1)+1)*(x(2)+1)*(x(3)+1),
    x(1)*x(2)+x(2)*x(3)+x(1)*x(3)+x(1),
    x(0)+x(1)+x(2)+x(3)+x(4)+x(5),
    all_monomials_of_degree_d(3,[x(i) for i in xrange(10)]),
    power_set([x(i) for i in xrange(10)]),
    random_poly,
    random_set,
    Polynomial(all_monomials_of_degree_d(3,[x(i) for i in xrange(10)])) +
        Polynomial(power_set([x(i) for i in xrange(10)])),
    Polynomial(power_set([x(i) for i in xrange(10)]))+1
    ]
    for colored in [True,False]:
        if colored:
            colored_suffix="_colored"
        else:
            colored_suffix=""
        for format in ["png", "svg"]:
            for (i,p) in enumerate(polynomials):

                #dot_file=str(i) +colored_suffix+".dot"
                #f=open(dot_file, "w")
                #f.write(dot)
                #f.close()
                out_file=str(i)+colored_suffix+"."+format
                plot(p, out_file, colored=colored,format=format)
                #system("dot -Tpng -o "+png_file+" " + dot_file)

if __name__ == '__main__':
    main()
