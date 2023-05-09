import itertools
import math
from matplotlib import pyplot as plt

class AlphaComplex:
    def __init__(self, data, distance) -> None:
        self.data = data
        self.distance = distance
        self.complex = None
        self.one_simplex = []
        self.two_simplex = []
        self.not_gabriel = []
        self.binary_intersections = []
        self.triple_intersections = []

    @staticmethod
    def create_circle_w_two_point(p1, p2):

        x0, y0 = (p1[0] + p2[0])/2, (p1[1] + p2[1])/2
        r = math.sqrt((x0 - p1[0])**2 + (y0 - p1[1])**2)

        return (x0, y0), r

    @staticmethod
    def create_circle_w_three_point(p1, p2, p3):
        x1, y1 = p1[0], p1[1]
        x2, y2 = p2[0], p2[1]
        x3, y3 = p3[0], p3[1]

        try:
            x0 = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (
                    y1 - y2)) / (
                         2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)))
            y0 = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (
                    x2 - x1)) / (
                         2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)))

            r = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

            return (x0, y0), r
        except ZeroDivisionError:
            return False, False

    def control_in_circle(self, circle_center, circle_r):

        c_x, c_y = circle_center
        r = circle_r

        for p in self.data:

            if c_x - r < p[0] < c_x + r and c_y - r < p[1] < c_y + r:
                if (c_x - p[0]) ** 2 + (c_y - p[1]) ** 2 + 1e-10 < r ** 2:
                    return False
                else:
                    continue
            else:
                continue
        return True

    def triple_point_selection(self):
        combs = list(itertools.combinations(self.data, 3))
        unique_combs = [list(comb) for comb in combs if comb[0] != comb[1] and comb[0] != comb[2] and comb[1] != comb[2]]
        return unique_combs

    def delaunay_tri(self):
        triple = self.triple_point_selection()
        complex_d = []
        circles = []
        for p3 in triple:
            c_center, r = self.create_circle_w_three_point(p3[0], p3[1], p3[2])
            if c_center is not False:
                if self.control_in_circle(c_center, r) is True:
                    complex_d.append(p3)
                    circles.append([c_center, r])
                else:
                    continue

        return complex_d, circles

    def control_gabriel(self, p3):
        combs = list(itertools.combinations(p3, 2))
        unique_combs = [list(comb) for comb in combs if comb[0] != comb[1]]
        not_gabriel = []
        for p2 in unique_combs:
            c_center, r = self.create_circle_w_two_point(p2[0], p2[1])
            if self.control_in_circle(c_center, r) is not True:
                not_gabriel.append(p2)
        if not_gabriel:
            return not_gabriel
        else:
            return False


    def alpha_complex(self):
        com = self.delaunay_tri()[0]
        for p3 in com:
            if self.control_gabriel(p3) is False:
                self.triple_intersections.append(p3)
            else:
                p2 = self.control_gabriel(p3)
                p3.remove(p2[0][0])
                p3.remove(p2[0][1])
                self.not_gabriel.append([p3[0], p2[0]])

    def draw_complex(self):
        self.alpha_complex()
        fig, ax = plt.subplots(figsize=(16, 9))
        for tri in self.triple_intersections:
            x = [p[0] for p in tri]
            y = [p[1] for p in tri]
            x.append(x[0])
            y.append(y[0])
            ax.fill(x, y, facecolor='gray', alpha=0.2)
            edge = plt.Line2D(x, y, color='black')
            ax.add_artist(edge)


        print(self.not_gabriel)
        for line in self.not_gabriel:
            x = [line[0][0], line[1][0][0]]
            y = [line[0][1], line[1][0][1]]
            edge = plt.Line2D(x, y, color='black')
            ax.add_artist(edge)
            x1 = [line[0][0], line[1][1][0]]
            y1 = [line[0][1], line[1][1][1]]
            edge1 = plt.Line2D(x1, y1, color='black')
            ax.add_artist(edge1)

        x_cor = [p[0] for p in self.data]
        y_cor = [p[1] for p in self.data]
        ax.plot(x_cor, y_cor, 'yo')
        ax.set_xlim(min(x_cor) - 2 * self.distance, max(x_cor) + 2 * self.distance)
        ax.set_ylim(min(y_cor) - 2 * self.distance, max(y_cor) + 2 * self.distance)

        plt.show()


data5 = [[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]]
a = AlphaComplex(data5, None)
a.draw_complex()
