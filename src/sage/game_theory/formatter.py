class Formatter():

    def __init__(self, messy_crap, game=None):
        self.messy_crap = messy_crap
        self.game = game

    def format_lrs(self):
        from fractions import Fraction
        p2_strategies = []
        p1_strategies = []
        for i in self.messy_crap:
            if i.startswith('2'):
                nums = [float(Fraction(k)) for k in i.split()]
                p2_strategies.append(nums[1:-1])
            elif i.startswith('1'):
                nums = [float(Fraction(k)) for k in i.split()]
                p1_strategies.append(nums[1:-1])

        return zip(p1_strategies, p2_strategies)

    def format_gambit(self):
        nice_stuff = []
        for gambitstrategy in self.messy_crap:
            gambitstrategy = eval(str(gambitstrategy)[str(gambitstrategy).index("["): str(gambitstrategy).index("]") + 1])
            profile = [gambitstrategy[:len(self.game.players[int(0)].strategies)]]
            for player in list(self.game.players)[1:]:
                previousplayerstrategylength = len(profile[-1])
                profile.append(gambitstrategy[previousplayerstrategylength: previousplayerstrategylength + len(player.strategies)])
            nice_stuff.append(profile)

        return nice_stuff
