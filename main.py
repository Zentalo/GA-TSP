import random
import math
import numpy as np

class Life(object):
    def __init__(self, gene=None):
        self.gene = gene
        self.fit = -1.0

    def getFit(self):
        return self.fit

    def setFit(self, f):
        self.fit = f

class TSP(object):
    def __init__(self, gene_length=8, cross_rate=0.9, mutate_rate=0.15, life_num=80, iter_max=30000):
        self.cross_rate = cross_rate  # 交叉率

        self.mutate_rate = mutate_rate  # 突变率
        self.mutate_times = 0         # 突变次数

        self.fit_sum = 0.0              # 当代适应和

        self.best_one = None            # 最优Gene

        self.lives = list()                # Gene全体
        self.life_num = life_num        # Gene总数
        self.gene_length = gene_length  # Gene长度

        self.iter = 0                   # 当前代数
        self.iter_max = iter_max        # 最大代数
        self.min_distance = list()          # 最小距离

        self.same = 1
        for i in range(self.life_num):
            self.lives.append(Life(self.init_life()))

        data = np.loadtxt('distance.txt', delimiter=',')
        matrix = data.reshape(8, 8)
        self.distance_matrix = matrix

    def init_life(self):
        ge = list(range(self.gene_length))
        random.shuffle(ge)
        return ge

    def distance(self, l):
        distance = 0
        for i in range(-1, self.gene_length - 1):
            p1 = l.gene[i]
            p2 = l.gene[i + 1]
            distance += self.distance_matrix[p1][p2]
        return distance

    def evolve(self):
        while (self.iter < self.iter_max ):
            # calculate the sum of Fitness
            self.fit_sum = .0
            self.best_one = Life()
            self.best_one.setFit(-1.0)
            for l in self.lives:
                l.setFit(1.0/self.distance(l))
                if l.getFit() > self.best_one.getFit():
                    self.best_one = l
                self.fit_sum += l.getFit()
            # evolution
            newlives = list()
            newlives.append(self.best_one)
            while(len(newlives)<self.life_num):
                newlives.append(self.birth())
            self.lives = newlives
            self.iter+=1
            self.min_distance.append(self.distance(self.best_one))
            if (self.iter % 50 == 0):
                print("Iter Times: %d,Best Distance: %d" % (self.iter, self.distance(self.best_one)))
                print(self.best_one.gene)
                print("Same:%d" %(self.same))
            if (self.iter>=2):
                if(self.min_distance[-1]==self.min_distance[-2]):
                    self.same+=1
                else:
                    self.same=1
            if (self.same>1000):
                return

    def birth(self):
        life1 = Life()
        life2 = Life()
        # select 1st
        r1 = random.uniform(0, self.fit_sum)
        for l in self.lives:
            r1 -= l.getFit()
            if r1 <= 0:
                life1 = l
        # select 2nd
        r2 = random.uniform(0, self.fit_sum)
        for l in self.lives:
            r2 -= l.getFit()
            if r2 <= 0:
                life2 = l
        # cross
        r3 = random.random()
        if r3 < self.cross_rate:
            gene = self.cross(life1, life2)
        else:
            gene = life1.gene
        # mutate
        r4 = random.random()
        if r4 < self.mutate_rate:
            gene = self.mutate(gene)
            self.mutate_times += 1
        return Life(gene)

    def cross(self, l1, l2):
        int1 = random.randint(0, self.gene_length-1)
        int2 = l2.gene[(l2.gene.index(int1)+1) % self.gene_length]
        index1 = l1.gene.index(int1)
        index2 = l1.gene.index(int2)
        new = list()
        if (index1 < index2):
            for i in range(0, index1+1):
                new.append(l1.gene[i])
            for i in range(index1+1, index2+1):
                new.append(l1.gene[index2+index1+1-i])
            for i in range(index2+1, self.gene_length):
                new.append(l1.gene[i])
        else:
            for i in range(0, index2+1):
                new.append(l1.gene[i])
            for i in range(index2+1,index1+1):
                new.append(l1.gene[index1+index2+1-i])
            for i in range(index1+1, self.gene_length):
                new.append(l1.gene[i])
        return new

    def mutate(self,ge):
        index1 = random.randint(0, self.gene_length-1)
        index2 = random.randint(0, self.gene_length-1)
        new = list()
        if (index1 < index2):
            for i in range(0, index1 + 1):
                new.append(ge[i])
            for i in range(index1 + 1, index2 + 1):
                new.append(ge[index2 + index1 + 1 - i])
            for i in range(index2 + 1, self.gene_length):
                new.append(ge[i])
        else:
            for i in range(0, index2 + 1):
                new.append(ge[i])
            for i in range(index2 + 1, index1 + 1):
                new.append(ge[index1 + index2 + 1 - i])
            for i in range(index1 + 1, self.gene_length):
                new.append(ge[i])
        return new

tsp=TSP()
tsp.evolve()

