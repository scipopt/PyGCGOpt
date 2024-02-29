from pygcgopt import Model, Score

import pytest

import os

class PyScore(Score):
    def scorecalculate(self, partialdec):
        nmasterconss = float(partialdec.getNMasterconss())
        nconss = float(partialdec.getNConss())
        score = 1 - nmasterconss/nconss

        return {"scorevalue": score}

@pytest.mark.parametrize("lp_file", [
    "instances_bpp/N1C1W4_M.BPP.lp", "instances_bpp/N1C2W2_O.BPP.lp"
])
def test_score(lp_file):
    dirname = os.path.dirname(__file__)
    lp_file = os.path.join(dirname, lp_file)

    m = Model()

    proxyScore = PyScore()
    m.includeScore(proxyScore, "pyscore", "python", "Python score test")
    assert "pyscore" in m.listScores()

    m.readProblem(lp_file)
    m.detect()

    partdecs = m.listDecompositions()

    for partdec in partdecs:
        assert partdec.getScore(proxyScore) == 1 - float(partdec.getNMasterconss())/float(partdec.getNConss())
