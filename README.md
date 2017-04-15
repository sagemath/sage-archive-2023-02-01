Mac Lane infrastructure for discrete (pseudo-)valuations which should work on an unmodified Sage 8.0 or later.

```
$ sage
sage: from mac_lane import *
sage: pAdicValuation(QQ, 2)
2-adic valuation
```

To run the included tests, execute `sage -tp --optional=sage,standalone mac_lane/`.
