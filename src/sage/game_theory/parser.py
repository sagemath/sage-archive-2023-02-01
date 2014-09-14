class Parser():
    r"""
    A class for parsing the outputs of different algorithms so that they can
    be used be used computationally
    """

    def __init__(self, mess, game=None):
        self.mess = mess
        self.game = game

    def format_lrs(self):
        from sage.misc.sage_eval import sage_eval
        p2_strategies = []
        p1_strategies = []
        for i in self.mess:
            if i.startswith('2'):
                nums = [sage_eval(k) for k in i.split()]
                p2_strategies.append(tuple(nums[1:-1]))
            elif i.startswith('1'):
                nums = [sage_eval(k) for k in i.split()]
                p1_strategies.append(tuple(nums[1:-1]))

        return [list(a) for a in zip(p1_strategies, p2_strategies)]
