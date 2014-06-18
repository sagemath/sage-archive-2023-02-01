class Parser():
    r"""
    A class for parsing the outputs of different algorithms so that they can
    be used be used computationaly
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
                p2_strategies.append(nums[1:-1])
            elif i.startswith('1'):
                nums = [sage_eval(k) for k in i.split()]
                p1_strategies.append(nums[1:-1])

        return zip(p1_strategies, p2_strategies)

    def format_gambit(self):
        nice_stuff = []
        for gambitstrategy in self.mess:
            gambitstrategy = eval(str(gambitstrategy)[str(gambitstrategy).index("["): str(gambitstrategy).index("]") + 1])
            profile = [gambitstrategy[:len(self.game.players[int(0)].strategies)]]
            for player in list(self.game.players)[1:]:
                previousplayerstrategylength = len(profile[-1])
                profile.append(gambitstrategy[previousplayerstrategylength: previousplayerstrategylength + len(player.strategies)])
            nice_stuff.append(profile)

        return nice_stuff
