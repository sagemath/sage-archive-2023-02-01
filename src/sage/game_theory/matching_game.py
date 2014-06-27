from sage.structure.sage_object import SageObject


class MatchingGame(SageObject):
    r"""
    """
    def __init__(self, num_pairs=0):
        r"""
        """
        self.suitors = {}
        self.reviewers = {}
        for i in range(num_pairs):
            self.add_suitor()
            self.add_reviewer()

    def _repr_(self):
        r"""
        """
        pass

    def _latex(self):
        r"""
        """
        pass

    def _is_complete(self):
        pass

    def add_suitor(self, name=False):
        if name is False:
            name = len(self.suitors)
        new_suitor = _Player(name, 'suitor', len(self.reviewers))
        self.suitors[new_suitor] = new_suitor.pref
        for r in self.reviewers:
            self.reviewers[r] = [-1 for s in self.suitors]

    def add_reviewer(self, name=False):
        if name is False:
            name = len(self.reviewers)
        new_reviewer = _Player(name, 'reviewer', len(self.suitors))
        self.reviewers[new_reviewer] = new_reviewer.pref
        for s in self.suitors:
            self.suitors[s] = [-1 for r in self.reviewers]


class _Player():
    def __init__(self, name, player_type, len_pref):
        self.name = name
        self.type = player_type
        self.pref = [-1 for i in range(len_pref)]

    def __hash__(self):
        return hash(self.name)

