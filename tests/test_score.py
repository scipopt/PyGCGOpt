from pygcgopt import Model, Score

class PyScore(Score):
    def scorecalculate(self, partialdec):
        nconss = float(partialdec.getNConss())
        nmasterconss = float(partialdec.getNMasterconss())
        score = 1 - nmasterconss/nconss

        return {"scorevalue": score}

def test_pyscore():

    m = Model()

    proxyScore = PyScore()
    m.includeScore(proxyScore, "pyscore", "python", "Python score test")

    m.readProblem("/Users/lentz/Desktop/BR7_56.lp.gz")
    m.detect()

    partdecs = m.listDecompositions()

    for i in partdecs:
        print(str(i.getScore(proxyScore)))

    return m.listScores()

if __name__ == "__main__":
   print(test_pyscore())
