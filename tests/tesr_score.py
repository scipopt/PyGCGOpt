from pygcgopt import Model, Score

class PyScore(Score):
    def scorecalculate(self, partialdecid):
        return {"scorevalue": 0.5}

def test_pyscore():

    m = Model()

    proxyScore = PyScore()
    m.includeScore(proxyScore, "pyscore", "python", "Python score test")

    m.readProblem("/Users/lentz/Desktop/BR7_56.lp.gz")

    m.optimize()

if __name__ == "__main__":
   test_pyscore()
