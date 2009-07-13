************
Introduction
************

Explorer ce tutoriel en entier devrait vous prendre au maximum trois à
quatre heures. Vous pouvez le lire en version HTML ou PDF, ou encore
l'explorer interactivement à l'intérieur de Sage en cliquant sur
``Help`` puis sur ``Tutorial`` depuis le *notebook* (il est possible que
vous tombiez sur la version en anglais).

Sage est écrit en grande partie en Python, mais aucune connaissance de
Python n'est nécessaire pour lire ce tutoriel. Par la suite, vous
souhaiterez sans doute apprendre Python, et il existe pour cela de
nombreuses ressources libres d'excellente qualité, dont [PyT]_ et
[Dive]_. Mais si ce que vous voulez est découvrir rapidement Sage, ce
tutoriel est le bon endroit où commencer. Voici quelques exemples :

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223

    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]

    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)

    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]

    sage: E = EllipticCurve([1,2,3,4,5]);
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1

    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    1/(I*sqrt(3) + 5/9*sqrt(73) + 3/4)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{1}{I \, \sqrt{3} + \frac{5}{9} \, \sqrt{73} + \frac{3}{4}}

Installation
============

Si Sage n'est pas installé sur votre ordinateur, vous pouvez essayer
quelques commandes en ligne à l'adresse http://www.sagenb.org.

Des instructions pour installer Sage sur votre ordinateur sont
disponibles dans le guide d'installation (*Installation Guide*), dans
la section documentation de la page web principale de Sage [Sage]_.
Nous nous limiterons ici à deux remarques.

#. La version téléchargeable de Sage vient avec ses dépendances.
   Autrement dit, bien que Sage utilise Python, IPython, PARI, GAP,
   Singular, Maxima, NTL, GMP, etc., vous n'avez pas besoin de les
   installer séparément, ils sont fournis dans la distribution Sage. En
   revanche, pour utiliser certaines des fonctionnalités de Sage, par
   exemple Macaulay ou KASH, il vous faudra d'abord installer le paquet
   optionnel correspondant, ou du moins avoir le logiciel correspondant
   installé sur votre ordinateur. Macaulay et KASH sont disponibles pour
   SAGE sous forme de modules optionnels (pour une liste de tous les
   paquets optionnels disponibles, tapez ``sage -optional`` ou visitez
   la page *Download* (Téléchargement) du site web de Sage).

#. La version binaire pré-compilée de Sage (disponible sur le site web)
   est souvent plus facile et plus rapide à installer que la
   distribution en code source. Pour l'installer, décompressez
   l'archive et lancez simplement le programme ``sage``.

Les différentes manières d'utiliser Sage
========================================

Il y a plusieurs façons d'utiliser Sage.

-  **Interface graphique (« notebook ») :** voir la section sur le
   *Notebook* du manuel de référence, et :ref:`section-notebook`
   ci-dessous ;

-  **Ligne de commande :** voir :ref:`chapter-interactive_shell` ;

-  **Programmes :** en écrivant des programmes interprétés ou
   compilés en Sage (voir :ref:`section-loadattach` et :ref:`section-compile`) ;

-  **Scripts :** en écrivant des programmes Python indépendants qui font
   appel à la bibliothèque Sage (voir :ref:`section-standalone`).


Objectifs à long terme de Sage
===============================

-  **Étre utile :** le public visé par Sage comprend les étudiants  (du lycée
   au doctorat), les enseignants et les chercheurs en mathématiques.
   Le but est de fournir un logiciel qui permette d'explorer toutes
   sortes de constructions mathématiques et de faire des expériences
   avec, en algèbre, en géométrie, en arithmétique et théorie des
   nombres, en analyse, en calcul numérique, etc. Sage facilite
   l'expérimentation interactive avec des objets mathématiques.

-  **Être efficace :** c'est-à-dire rapide. Sage fait appel à des
   logiciels matures et soigneusement optimisés comme GMP, PARI, GAP et
   NTL, ce qui le rend très rapide pour certaines opérations.

-  **Être libre/open-source :** le code source doit être disponible
   librement et lisible, de sorte que les utilisateurs puissent
   comprendre ce que fait le système et l'étendre facilement. Tout
   comme les mathématiciens acquièrent une compréhension plus profonde
   d'un théorème en lisant sa preuve soigneusement, ou simplement en la
   parcourant, les personnes qui font des calculs devraient être en
   mesure de comprendre comment ceux-ci fonctionnent en lisant un code
   source documenté. Si vous publiez un article dans lequel vous
   utilisez Sage pour faire des calculs, vous avez la garantie que vos
   lecteurs auront accès librement à Sage et à son code source, et vous
   pouvez même archiver et redistribuer vous-même la version de Sage que
   vous utilisez.

-  **Être facile à compiler :** le code source de Sage devrait être
   facile à compiler pour les utilisateurs de Linux, d'OS X et de
   Windows. Cela rend le système plus flexible pour les utilisateurs qui
   souhaiteraient le modifier.

-  **Favoriser la coopération :** fournir des interfaces robustes à
   la plupart des autres systèmes de calcul formel, notamment PARI, GAP,
   Singular, Maxima, KASH, Magma, Maple et Mathematica. Sage cherche à
   unifier et étendre les logiciels existants.

-  **Être bien documenté :** tutoriel, guide du programmeur, manuel de
   référence, guides pratiques, avec de nombreux exemples et une
   discussion des concepts mathématiques sous-jacents.

-  **Être extensible  :** permettre de définir de nouveaux types de
   données ou des types dérivés de types existants, et d'utiliser du
   code écrit dans différents langages.

-  **Être convivial :** il doit être facile de comprendre quelles
   fonctionnalités sont disponibles pour travailler avec un objet donné,
   et de consulter la documentation et le code source. Également,
   arriver à un bon niveau d'assistance utilisateur.

.. [Dive] Dive into Python, Freely available online at
          http://diveintopython.org

.. [PyT] The Python Tutorial, http://www.python.org/

.. [Sage] Sage, http://www.sagemath.org
