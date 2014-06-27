from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer


class MatchingGame(SageObject):
    r"""
    """
    def __init__(self, generator):
        r"""
        EXAMPLES:

        quick test. ::

            sage: suitr_pref = {'J': ['A', 'D', 'C', 'B'],
            ....:               'K': ['A', 'B', 'C', 'D'],
            ....:               'L': ['B', 'D', 'C', 'A'],
            ....:               'M': ['C', 'A', 'B', 'D']}
            sage: reviewr_pref = {'A': ['L', 'J', 'K', 'M'],
            ....:                 'B': ['J', 'M', 'L', 'K'],
            ....:                 'C': ['K', 'M', 'L', 'J'],
            ....:                 'D': ['M', 'K', 'J', 'L']}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.suitors
            ['K', 'J', 'M', 'L']
            sage: m.reviewers
            ['A', 'C', 'B', 'D']

        works for numbers too? ::

            sage: suit = {0: [3, 4],
            ....:         1: [3, 4]}
            sage: revr = {3: [0, 1],
            ....:         4: [1, 0]}
            sage: g = MatchingGame([suit, revr])
        """
        self.suitors = []
        self.reviewers = []
        if type(generator) is Integer:
            for i in range(generator):
                self.add_suitor()
                self.add_reviewer()
        if type(generator[0]) is dict and type(generator[1]) is dict:
            self._dict_game(generator[0], generator[1])
        else:
            raise TypeError("generator must be integer or a list of 2 dictionaries.")

    def _dict_game(self, suitor_dict, reviwer_dict):
        for i in suitor_dict:
            self.add_suitor(i)
        for k in reviwer_dict:
            self.add_reviewer(k)

        for i in self.suitors:
            i.pref = suitor_dict[i.name]
        for k in self.reviewers:
            k.pref = reviwer_dict[k.name]

    def _repr_(self):
        r"""
        """
        pass

    def _latex(self):
        r"""
        """
        pass

    def _is_complete(self):
        r"""
        """
        if len(self.suitors) != len(self.reviewers):
            raise ValueError("Must have the same number of reviewers as suitors")

        for suitor in self.suitors:
            if sorted(suitor.pref) != sorted(self.reviewers):
                raise ValueError("Preferences incomplete")

        for reviewer in self.reviewers:
            if sorted(reviewer.pref) != sorted(self.suitors):
                raise ValueError("Preferences incomplete")

        return True

    def add_suitor(self, name=False):
        r"""
        """
        if name is False:
            name = len(self.suitors)
        new_suitor = _Player(name, 'suitor', len(self.reviewers))
        self.suitors.append(new_suitor)
        for r in self.reviewers:
            r.pref = [-1 for s in self.suitors]

    def add_reviewer(self, name=False):
        r"""
        """
        if name is False:
            name = len(self.reviewers)
        new_reviewer = _Player(name, 'reviewer', len(self.suitors))
        self.reviewers.append(new_reviewer)
        for s in self.suitors:
            s.pref = [-1 for r in self.reviewers]


class _Player():
    def __init__(self, name, player_type, len_pref):
        self.name = name
        self.type = player_type
        self.pref = [-1 for i in range(len_pref)]
        self.partner = False

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return repr(self.name)

