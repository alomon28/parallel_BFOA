import copy
import math
from multiprocessing import Manager, Pool, managers
from pickle import FALSE, TRUE
from evaluadorBlosum import evaluadorBlosum
import numpy
from fastaReader import fastaReader
import random
from copy import copy
import copy
import concurrent.futures


class bacteria():
    

    def __init__(self, numBacterias):
        # manager = Manager()
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))

    def resetListas(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
        
        
  
    def cuadra(self, numSec, poblacion):
        #ciclo para recorrer poblacion
        for i in range(len(poblacion)):
            #obtiene las secuencias de la bacteria
            bacterTmp = poblacion[i]
            # print("bacterTmp: ", bacterTmp)
            bacterTmp = list(bacterTmp)
            # print("bacterTmp: ", bacterTmp)
            bacterTmp = bacterTmp[:numSec]
            # obtiene el tama�o de la secuencia m�s larga
            maxLen = 0
            for j in range(numSec):
                if len(bacterTmp[j]) > maxLen:
                    maxLen = len(bacterTmp[j])
                    #rellena con gaps las secuencias m�s cortas
                    for t in range(numSec):
                        gap_count = maxLen - len(bacterTmp[t])
                        if gap_count > 0:
                            bacterTmp[t].extend(["-"] * gap_count)
                            #actualiza la poblacion
                            poblacion[i] = tuple(bacterTmp)
                            
            
        
        
        
        



    """metodo que recorre la matriz y elimina las columnas con gaps en todos los elementos"""
    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1
  
                
            
        """metodo para eliminar un elemento especifico en cada secuencia"""
    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]



    """metodo para saber si alguna columna de self.matrix tiene  gap en todos los elementos"""
    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True
    
    

    def tumbo(self, numSec, poblacion, numGaps):
        #inserta un gap en una posicion aleatoria de una secuencia aleatoria
        #recorre la poblacion
        for i in range(len(poblacion)):
            #obtiene las secuencias de la bacteria
            bacterTmp = poblacion[i]
            bacterTmp = list(bacterTmp)
            # bacterTmp = bacterTmp[:numSec]
            #ciclo para insertar gaps
            for j in range(numGaps):
                #selecciona secuencia
                seqnum = random.randint(0, len(bacterTmp)-1)
                #selecciona posicion
                pos = random.randint(0, len(bacterTmp[seqnum]))
                part1 = bacterTmp[seqnum][:pos]
                part2 = bacterTmp[seqnum][pos:]
                temp = part1 + ["-"] + part2
                bacterTmp[seqnum] = temp
                #actualiza la poblacion
                poblacion[i] = tuple(bacterTmp)
        
       
            
    def creaGranListaPares(self, poblacion):   
        # granListaPares = list(range(len(poblacion)))
        #ciclo para recorrer poblacion
        for i in range(len(poblacion)):  #recorre poblacion
            pares = list()
            bacterTmp = poblacion[i]
            bacterTmp = list(bacterTmp)
            #ciclo para recorrer secuencias
            for j in range(len(bacterTmp)):     #recorre secuencias de bacteria
                column = self.getColumn(bacterTmp, j)
                pares = pares + self.obtener_pares_unicos(column)
            self.granListaPares[i] = pares
            # print("Bacteria: ", i, " Pares: ", pares)
            
        # return self.granListaPares
    


    def evaluaFila(self, fila, num, cache):
        score = 0
        for par in fila:
            score += cache.get((par[0], par[1]), -8)
        self.blosumScore[num] = score

    
    def evaluaBlosum(self):
        evaluador = evaluadorBlosum()
        cache = evaluador.cache
        def wrapped_evaluaFila(args):
            fila, num = args
            return self.evaluaFila(fila, num, cache)
        with concurrent.futures.ThreadPoolExecutor() as executor:
            list(executor.map(wrapped_evaluaFila, [(copy.deepcopy(self.granListaPares[i]), i) for i in range(len(self.granListaPares))]))
    #Cambiamos a ThreadPoolExecutor para mantener acceso al objeto cache,
    # ya que multiprocessing.Pool no puede compartir objetos complejos fácilmente sin serialización especial.



    def getColumn(self, bacterTmp, colNum):
        column = []
        #obtiene las secuencias de la bacteria
        # bacterTmp = poblacion[bactNum]
        # bacterTmp = list(bacterTmp)
        #obtiene el caracter de cada secuencia en la columna
        for i in range(len(bacterTmp)):
            column.append(bacterTmp[i][colNum])
        return column
            
        
            
    

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)  

    #------------------------------------------------------------Atract y Repel lineal
    
  


    def compute_diff(self, args):
        indexBacteria, otherBlosumScore, self.blosumScore, d, w = args
        diff = (self.blosumScore[indexBacteria] - otherBlosumScore) ** 2.0
        self.NFE[indexBacteria] += 1
        return d * numpy.exp(w * diff)

    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        with Pool() as pool:
            args = [(indexBacteria, otherBlosumScore, self.blosumScore, d, w) for otherBlosumScore in self.blosumScore]
            results = pool.map(self.compute_diff, args)
            pool.close()  # Close the pool to prevent any more tasks from being submitted
            pool.join()   # Wait for the worker processes to exit
    
        total = sum(results)
    
        if atracTrue:
            self.tablaAtract[indexBacteria] = total
        else:
            self.tablaRepel[indexBacteria] = total
        

  
    def creaTablaAtract(self, poblacion, d, w):                   #lineal
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria,d, w, TRUE)
            # print("invocando indexBacteria numero: ", indexBacteria)
        # print("tablaAtract: ", self.tablaAtract)

    def creaTablaRepel(self, poblacion, d, w):                   #lineal
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria,d, w, FALSE)
            # print("invocando indexBacteria numero: ", indexBacteria)
        # print("tablaAtract: ", self.tablaAtract)
    
    def creaTablasAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        #invoca ambos metodos en paralelo
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(self.creaTablaAtract, poblacion, dAttr, wAttr)
            executor.submit(self.creaTablaRepel, poblacion, dRepel, wRepel)
            



            #-----------------------------------------------------------
            
    def creaTablaInteraction(self):
        #llena la tabla con la suma de atract y repel
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def creaTablaFitness(self):
        #llena la tabla con la suma de interaction y blosumScore
        for i in range(len(self.tablaInteraction)):
            valorBlsm = self.blosumScore[i]
            valorInteract = self.tablaInteraction[i]
            #suma ambos valores
            valorFitness =  valorBlsm + valorInteract
            
            self.tablaFitness[i] = valorFitness
    
    def getNFE(self):
        return sum(self.NFE)
        
        
    def obtieneBest(self, globalNFE):
        bestIdx = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] > self.tablaFitness[bestIdx]:
                bestIdx = i
        print("-------------------   Best: ", bestIdx, " Fitness: ", self.tablaFitness[bestIdx], "BlosumScore ",  self.blosumScore[bestIdx], "Interaction: ", self.tablaInteraction[bestIdx], "NFE: ", globalNFE)
        return bestIdx, self.tablaFitness[bestIdx]

    def replaceWorst(self, poblacion, best):
        worst = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] < self.tablaFitness[worst]:
                worst = i
        # print("Worst: ", worst,  "Blosum ",self.blosumScore[worst], "Fitness: ", self.tablaFitness[worst], "BlosumScore: ", self.blosumScore[worst], "Atract: ", self.tablaAtract[worst], "Repel: ", self.tablaRepel[worst], "Interaction: ", self.tablaInteraction[worst])
        #reemplaza la bacteria peor por una copia de la mejor
        poblacion[worst] = copy.deepcopy(poblacion[best])
        
        