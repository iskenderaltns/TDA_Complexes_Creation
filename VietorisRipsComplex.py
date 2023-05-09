import itertools
import math
import numpy as np
import matplotlib.pyplot as plt


class VietorisRipsComplex:
    def __init__(self, data, distance) -> None:
        # lets data X a set finite in R^d
        self.data = data
        self.distance = distance
        self.result = []
        self.binary_intersections = []
        self.triple_intersections = []

    @staticmethod
    def distance_btp(p1, p2):
        r0 = pow((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0.5)
        return r0

    def control_intersection_binary(self):
        combs = list(itertools.combinations(self.data, 2))
        unique_combs = [list(comb) for comb in combs if comb[0] != comb[1]]
        result_ = list(filter(lambda p: self.distance_btp(p[0], p[1]) < 2*self.distance, unique_combs))
        return result_

    def control_intersection_triple(self):
        self.binary_intersections = self.control_intersection_binary()
        result_ = []
        for p12 in self.binary_intersections:
            list_p1 = [p[1] for p in self.binary_intersections if p12[0] == p[0]]
            list_p2 = [p[1] for p in self.binary_intersections if p12[1] == p[0]]
            intersect = [p for p in list_p1 if p in list_p2]
            if intersect:
                k = [[p12[0], p12[1], intersect[i]] for i in range(len(intersect))]
                result_.extend(k)

        self.triple_intersections = result_
        return result_

    def draw_complex(self):
        self.control_intersection_triple()
        if self.triple_intersections:
            fig, ax = plt.subplots(figsize=(16, 9))
            for points in self.triple_intersections:
                for point in points:
                    color = np.random.rand(3)
                    circle = plt.Circle(point, self.distance, fill=False, edgecolor=color)
                    ax.add_artist(circle)
                x = [p[0] for p in points]
                y = [p[1] for p in points]
                ax.plot(x, y, 'bo')
                x.append(self.triple_intersections[0][0][0])
                y.append(self.triple_intersections[0][0][1])
                edge = plt.Line2D(x, y, color=np.random.rand(3, ))
                ax.add_artist(edge)


            for points in self.binary_intersections:
                for point in points:
                    color = np.random.rand(3)
                    circle = plt.Circle(point, self.distance, fill=False, edgecolor=color)
                    ax.add_artist(circle)
                x = [p[0] for p in points]
                y = [p[1] for p in points]
                ax.plot(x, y, 'bo')
                edge = plt.Line2D(x, y, color=np.random.rand(3, ))
                ax.add_artist(edge)

            x_cor = [p[0] for p in self.data]
            y_cor = [p[1] for p in self.data]
            ax.set_xlim(min(x_cor) - 2 * self.distance, max(x_cor) + 2 * self.distance)
            ax.set_ylim(min(y_cor) - 2 * self.distance, max(y_cor) + 2 * self.distance)

            self.drawData()
            self.drawSimplex()
            plt.show()
        else:
            self.drawData()
            self.drawSimplex()

    def drawSimplex(self):
        self.control_intersection_triple()
        fig, ax = plt.subplots(figsize=(16, 9))
        for points in self.binary_intersections:
            x = [p[0] for p in points]
            y = [p[1] for p in points]
            ax.plot(x, y, 'bo')
            edge = plt.Line2D(x, y, color=np.random.rand(3, ))
            ax.add_artist(edge)
        '''
        x = [i[0] for p in self.triple_intersections for i in p]
        y = [i[1] for p in self.triple_intersections for i in p]
        
        ax.fill(x, y, facecolor='grey', alpha=0.4)
        '''


        for i in self.triple_intersections:
            x = [p[0] for p in i]
            y = [p[1] for p in i]
            x.append(x[0])
            y.append(y[0])
            ax.fill(x, y, facecolor='grey', alpha=0.4)

        x_cor = [p[0] for p in self.data]
        y_cor = [p[1] for p in self.data]
        ax.plot(x_cor, y_cor, 'yo')
        ax.set_xlim(min(x_cor) - 2 * self.distance, max(x_cor) + 2 * self.distance)
        ax.set_ylim(min(y_cor) - 2 * self.distance, max(y_cor) + 2 * self.distance)

        plt.show()

    def drawData(self):
        x = [p[0] for p in self.data]
        y = [p[1] for p in self.data]
        fig, ax = plt.subplots(figsize=(16, 9))
        ax.plot(x, y, 'ro')
        ax.set_xlim(min(x) - 2 * self.distance, max(x) + 2 * self.distance)
        ax.set_ylim(min(y) - 2 * self.distance, max(y) + 2 * self.distance)
        plt.show()





data3 = [[4, 4], [1.691, 0], [6.309, 0], [11, 0], [1, -2], [3, 6], [4, 6], [5, 9], [1, 7], [7, 5], [6, 2], [4, 3], [3, 1], [-1, 0], [-2, 3], [-2, -1]]
data2 = [[4, 4], [1.691, 0], [6.309, 0], [11, 0]]
data1 = [[1, 2], [2, 3], [0, 0.5], [2.5, 1.3], [4, 8]]
a = VietorisRipsComplex(data3, math.sqrt(5.8))
a.draw_complex()


