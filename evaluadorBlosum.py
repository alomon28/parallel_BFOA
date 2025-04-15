import blosum as bl

class evaluadorBlosum():
    
    def __init__(self):
        matrix = bl.BLOSUM(62)
        self.matrix = matrix
        self.cache = self._precompute_matrix()

    def _precompute_matrix(self):
        chars = list(self.matrix.keys()) + ['-']
        cache = {}
        for a in chars:
            for b in chars:
                if a == '-' or b == '-':
                    cache[(a, b)] = -8
                else:
                    cache[(a, b)] = self.matrix[a][b]
        return cache

    def getScore(self, A, B):
        return self.cache.get((A, B), -8)


