import blosum as bl

class evaluadorBlosum():
    
    def __init__(self):
        matrix = bl.BLOSUM(62)
        
        self.matrix = matrix
        
    def showMatrix(self):
        print(self.matrix)
        
    def getScore(self, A, B):
        if A == "-" and B == "-":
            return -12  # penalización más alta por par gap-gap
        elif A == "-" or B == "-":
            return -6   # penalización menor que antes por un solo gap
        else:
            return self.matrix[A][B]
    
    pass


