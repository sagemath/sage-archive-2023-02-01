from sage.structure.sage_object import SageObject


class Matching_Game(SageObject):
    r"""
    """
    def __init__(self, num_pairs=0):
        r"""
        """
        self.suitors = []
        self.reviewers = []
        for i in range(num_pairs):
            self.add_player('suitor')
            self.add_player('reviewer')

    def _repr_(self):
        r"""
        """
        pass

    def _latex(self):
        r"""
        """
        pass

    def _is_complete(self):
        if len(self.suitors) != len(self.reviewers):
            raise ValueError("must have the same number of suitors as reviewers.")

        suitor_check = [any(i < 0 for i in s.pref) for s in self.suitors]
        reviewer_check = [any(i < 0 for i in r.pref) for r in self.reviewers]

        if any(suitor_check) or any(reviewer_check):
            raise ValueError("players preferences not completed")

        for suitor in self.suitors:
            if sorted(suitor.pref) != [i for i in range(len(self.reviewers))]:
                raise ValueError("suitor preferences not complete")

        for reviewer in self.reviewers:
            if sorted(reviewer.pref) != [i for i in range(len(self.suitors))]:
                raise ValueError("reviewer preferences not complete")

    def add_player(self, type):
        if type == 'suitor':
            self.suitors.append(_Player('suitor'))
            for reviewer in self.reviewers:
                reviewer.preference = [-1 for i in range(len(self.suitors))]
        elif type == 'reviewer':
            self.reviewers.append(_Player('reviewer'))
            for suitor in self.suitors:
                suitor.preference = [-1 for i in range(len(self.reviewers))]
        else:
            raise ValueError("type must be suitor or wife")


class _Player():
    def __init__(self, type):
        self.type = type
        self.pref = []