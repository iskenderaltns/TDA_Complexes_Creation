import random
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt


class DelaunayComplex:
    def __init__(self, data, command) -> None:
        self.data = data
        self.complex = None
        self.command = command

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
        except ZeroDivisionError or RuntimeWarning:
            return False, False

    def control_in_circle(self, circle_center, circle_r):

        c_x, c_y = circle_center
        r = circle_r

        for p in self.data:

            if c_x - r < p[0] < c_x + r and c_y - r < p[1] < c_y + r:
                if (c_x - p[0]) ** 2 + (c_y - p[1]) ** 2 + 1e-8 < r ** 2:
                    return False
                else:
                    continue
            else:
                continue
        return True

    def triple_point_selection(self):
        combs = list(itertools.combinations(self.data, 3))
        unique_combs = [list(comb) for comb in combs if
                        comb[0] != comb[1] and comb[0] != comb[2] and comb[1] != comb[2]]
        return unique_combs

    def create_circle_control_tri(self):
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

    @staticmethod
    def find_intersection_two_lines(line1, line2):
        x1, y1 = line1[0][0], line1[0][1]
        x2, y2 = line1[1][0], line1[1][1]
        x3, y3 = line2[0][0], line2[0][1]
        x4, y4 = line2[1][0], line2[1][1]

        determinant = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

        if determinant == 0:
            return False

        if determinant != 0 and (x1 - x3) * (y1 - y4) == (x1 - x4) * (y1 - y3):
            return False

        x = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / determinant
        y = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / determinant
        if (x == x1 and y == y1) or (x == x2 and y == y2) or (x == x3 and y == y3) or (x == x4 and y == y4):
            return False

        return x,y

    @staticmethod
    def control_2(t1, t2):
        count = 0
        inter = []
        for p in t1:
            if p in t2:
                inter.append(p)
                count += 1
            else:
                continue
        if count == 2:
            return inter
        else:
            return False

    def control_is_not_empty(self):
        com, circ = self.create_circle_control_tri()

        circ_point = {}
        for i in range(len(circ)):
            if circ.count(circ[i]) == 4:
                if str(circ[i]) not in circ_point:
                    circ_point[str(circ[i])] = com[i]
                else:
                    circ_point[str(circ[i])].extend(com[i])

        lines_inter = []
        for circ in circ_point:
            lines = list(itertools.combinations(circ_point[circ], 2))
            comb_lines = list(itertools.combinations(lines, 2))
            intersection = []
            for twoline in comb_lines:
                if self.find_intersection_two_lines(twoline[0], twoline[1]) is not False:
                    inter = self.find_intersection_two_lines(twoline[0], twoline[1])
                    if inter not in intersection:
                        intersection.append(inter)
                        lines_inter.append(twoline)

        return lines_inter



    def draw_complex(self):

        com, circ = self.create_circle_control_tri()
        fig, ax = plt.subplots(figsize=(16, 9))

        lines_intersection = self.control_is_not_empty()
        intersection = [points for lines in lines_intersection for points in lines[0]]
        for point in com:
            count = 0
            for p in point:
                if p in intersection:
                    count += 1

            if count == 1 or count == 0:

                x = [p[0] for p in point]
                y = [p[1] for p in point]
                x.append(point[0][0])
                y.append(point[0][1])
                edge = plt.Line2D(x, y, color=np.random.rand(3, ))
                ax.add_artist(edge)

            else:
                x1 = [p[0] for p in point if p in intersection]
                y1 = [p[1] for p in point if p in intersection]
                x = [p[0] for p in point if p not in intersection]
                y = [p[1] for p in point if p not in intersection]
                z1 = [x1[0], x[0], x1[1]]
                z2 = [y1[0], y[0], y1[1]]
                edge = plt.Line2D(z1, z2, color=np.random.rand(3, ))
                ax.add_artist(edge)


        if self.command == 'all':

            for c in circ:
                circle = plt.Circle(c[0], c[1], fill=False, edgecolor='black')
                ax.add_artist(circle)

        x = [p[0] for p in self.data]
        y = [p[1] for p in self.data]

        ax.plot(x, y, 'bo')
        x_cor = [p[0] for p in self.data]
        y_cor = [p[1] for p in self.data]

        ax.set_xlim(min(x_cor) - 3, max(x_cor) + 3)
        ax.set_ylim(min(y_cor) - 3, max(y_cor) + 3)

        plt.show()

    def send_for_alpha(self):
        return self.create_circle_control_tri()[1]


data2 = [[0, 0], [0, 1.1], [1, 0], [1, 1]]
data1 = [[0, 0], [0, 1], [1, 0], [1, 1], [0.5, 0.5], [0.5, 1.5], [1.5, 0.5], [1.5, 1.5], [1, 2], [2, 1]]
a = DelaunayComplex(data2, '')
a.draw_complex()
