import itertools
import math
import numpy as np
import matplotlib.pyplot as plt


class CechComplex:
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
        result_ = list(filter(lambda p: self.distance_btp(p[0], p[1]) < 2 * self.distance, unique_combs))
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

    @staticmethod
    def find_points_intersection(p1, p2, r1, r2):
        constant = (r1 ** 2) - (r2 ** 2) - (p1[0] ** 2) - (p1[1] ** 2) + (p2[0] ** 2) + (p2[1] ** 2)



        if p1[0] != p2[0]:
            m = (p2[1] - p1[1]) / (p2[0] - p1[0])
            a = 1 + (m ** 2)
            b = (-2 * ((constant / (2 * (p2[0] - p1[0]))) - p1[0]) * m) - (2 * p1[1])
            c = ((constant / (2 * (p2[0] - p1[0]))) - p1[0]) ** 2 + (p1[1] ** 2) - (r1 ** 2)
            if p1[1] != p2[1]:

                y1 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
                y2 = (-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
                x1 = ((constant - y1 * (2 * (p2[1] - p1[1]))) / (2 * (p2[0] - p1[0])))
                x2 = ((constant - y2 * (2 * (p2[1] - p1[1]))) / (2 * (p2[0] - p1[0])))

                return (x1, y1), (x2, y2)
            else:
                x = constant / (2 * (p2[0] - p1[0]))
                y1 = p1[1] + math.sqrt((r1 ** 2) - ((x - p1[0]) ** 2))
                y2 = p1[1] - math.sqrt((r1 ** 2) - ((x - p1[0]) ** 2))
                return (x, y1), (x, y2)
        else:
            y = constant / (2 * (p2[1] - p1[1]))
            x1 = p1[0] + math.sqrt((r1 ** 2) - ((y - p1[1]) ** 2))
            x2 = p1[0] - math.sqrt((r1 ** 2) - ((y - p1[1]) ** 2))
            return (x1, y), (x2, y)

    @staticmethod
    def find_intersection_two_lines(line1, line2):
        x1, y1, x2, y2 = line1
        x3, y3, x4, y4 = line2

        determinant = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

        if determinant == 0:
            return None

        if determinant != 0 and (x1 - x3) * (y1 - y4) == (x1 - x4) * (y1 - y3):
            return True

        x = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / determinant
        y = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / determinant

        return x, y

    @staticmethod
    def control_same_point(x, y, z):
        if abs(x - y) < 1e-12 and abs(y-z) < 1e-12 and abs(x - z) < 1e-12:
            return True
        else:
            return False


    def control_in_circle(self, points, x, y):
        for p in points:
            if (p[0] - x) ** 2 + (p[1] - y) ** 2 < self.distance ** 2:
                return True
            else:
                continue
        return False

    def control_is_not_empty(self):
        self.control_intersection_triple()
        point_intersection = []
        for points in self.triple_intersections:
            combs = list(itertools.combinations(points, 2))
            unique_combs = [list(comb) for comb in combs if comb[0] != comb[1]]
            k = [[self.find_points_intersection(i[0], i[1], self.distance, self.distance)] for i in unique_combs]
            combs1 = list(itertools.combinations(k, 2))
            unique_combs1 = [list(comb) for comb in combs1 if comb[0] != comb[1]]
            # print(unique_combs1)
            m = [self.find_intersection_two_lines((p[0][0][0][0], p[0][0][0][1], p[0][0][1][0], p[0][0][1][1]),
                                                  (p[1][0][0][0], p[1][0][0][1], p[1][0][1][0], p[1][0][1][1])) for p in
                 unique_combs1]
            if None in m or self.control_same_point(m[0][0], m[1][0], m[2][0]) is False or self.control_same_point(m[0][1], m[1][1], m[2][1]) is False:
                continue

            else:
                if self.control_in_circle(points, m[0][0], m[0][1]):
                    m.append(points)
                    point_intersection.append(m)
                else:
                    continue

        return point_intersection

    def draw_complex(self):
        com = self.control_is_not_empty()
        if com or self.binary_intersections:

            point_intersection = [p[0] for p in com]

            d = [i for p in com for i in p[-1]]

            fig, ax = plt.subplots(figsize=(16, 9))
            for p in d:
                color = np.random.rand(3)
                circle = plt.Circle(p, self.distance, color=color, fill=True, edgecolor=None)
                ax.add_artist(circle)


            x = [p[0] for p in d]
            y = [p[1] for p in d]
            ax.plot(x, y, 'bo')

            x1 = [p[0] for p in point_intersection]
            y1 = [p[1] for p in point_intersection]
            ax.plot(x1, y1, 'ro')

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
            ax.set_xlim(min(x_cor)-2*self.distance, max(x_cor)+2*self.distance)
            ax.set_ylim(min(y_cor)-2*self.distance, max(y_cor)+2*self.distance)

            self.drawData()
            self.drawSimplex()
            plt.show()
        else:
            return False

    def drawSimplex(self):
        com = self.control_is_not_empty()
        if com or self.binary_intersections:


            point_intersection = [p[0] for p in com]

            d = [i for p in com for i in p[-1]]

            fig, ax = plt.subplots(figsize=(16, 9))

            '''
            x = [p[0] for p in d]
            y = [p[1] for p in d]
            ax.plot(x, y, 'bo')

            x1 = [p[0] for p in point_intersection]
            y1 = [p[1] for p in point_intersection]
            ax.plot(x1, y1, 'ro')
            print(self.binary_intersections)
            
            '''


            for points in self.binary_intersections:

                x = [p[0] for p in points]
                y = [p[1] for p in points]
                ax.plot(x, y, 'bo')
                edge = plt.Line2D(x, y, color=np.random.rand(3, ))
                ax.add_artist(edge)

            g = []
            for i in range(0, len(d), 3):
                g.append(d[i:i + 3])

            for i in g:
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
        else:

            return False

    def drawData(self):
        x = [p[0] for p in self.data]
        y = [p[1] for p in self.data]
        fig, ax = plt.subplots(figsize=(16, 9))
        ax.plot(x, y, 'ro')
        ax.set_xlim(min(x) - 2 * self.distance, max(x) + 2 * self.distance)
        ax.set_ylim(min(y) - 2 * self.distance, max(y) + 2 * self.distance)
        plt.show()

    def result_for_alpha(self):
        point_intersection_three = self.control_is_not_empty()
        return point_intersection_three


data3 = [[4, 4], [1.691, 0], [6.309, 0], [11, 0], [1, -2]]
data2 = [[4, 4], [1.691, 0], [6.309, 0], [11, 0], [1, -2], [3, 6], [4, 6], [5, 9], [1, 7], [7, 5], [6, 2], [4, 3], [3, 1], [-1, 0], [-2, 3], [-2, -1]]
data1 = [[1, 2], [2, 3], [0, 0.5], [2.5, 1.3], [4, 8]]
data4 = [[0, 0], [0, 1], [1, 0], [1, 1], [0.5, 0.5], [0.5, 1.5], [1.5, 0.5], [1.5, 1.5], [1, 2], [2, 1]]
data5 = [[4, 4], [1.691, 0], [6.309, 0], [11, 0]]
a = CechComplex(data2, math.sqrt(5.8))
a.draw_complex()
