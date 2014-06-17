class Formatter():

    def __init__(self, messy_crap, game=None):
        self.messy_crap = messy_crap
        self.game = game

    def format_lrs():
        pass

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
