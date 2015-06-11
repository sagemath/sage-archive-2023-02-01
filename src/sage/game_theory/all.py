from sage.misc.lazy_import import lazy_import

import catalog as game_theory

lazy_import('sage.game_theory.cooperative_game', 'CooperativeGame')
lazy_import('sage.game_theory.normal_form_game', 'NormalFormGame')
lazy_import('sage.game_theory.matching_game', 'MatchingGame')
