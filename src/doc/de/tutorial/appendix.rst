******
Anhang
******

.. _section-precedence:

Binäre arithmetische Operatorrangfolge
========================================

Was ist ``3^2*4 + 2%5``? Der Wert (38) wird durch diese
"Operatorrangfolge-Tabelle" festgelegt. Die Tabelle unterhalb basiert
auf der Tabelle in Abschnitt § 5.15 des *Python Language Reference
Manual* von G. Rossum und F. Drake. Die Operatoren sind hier in
aufsteigender Ordnung der Bindungstärke aufgelistet.


==========================  =================
Operatoren                  Beschreibung
==========================  =================
or                          Boolesches oder
and                         Boolesches und
not                         Boolesches nicht
in, not in                  Zugehörigkeit
is, is not                  Identitätstest
>, <=, >, >=, ==, !=        Vergleich
+, -                        Addition, Subtraktion
\*, /, %                    Multiplikation, Division, Restbildung
\*\*, ^                     Exponentiation
==========================  =================

Um also ``3^2*4 + 2%5`` zu berechnen klammert Sage den Ausdruck in
folgender Weise: ``((3^2)*4) + (2%5)``. Es wird daher zuerst ``3^2``,
was ``9`` ist, dann wird sowohl ``(3^2)*4`` als auch ``2%5`` berechnet,
und schließlich werden diese beiden Werte addiert.
