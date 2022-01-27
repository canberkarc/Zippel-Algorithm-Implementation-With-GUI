#!/usr/bin/python
# -*- coding: latin-1 -*-
import time
from pygame.locals import *
import random
import pygame
from sympy import *

x, y, z = symbols('x y z')
funcs = [x ** 4 * y + x ** 2 * y ** 2 * z + x * y * z + 3, x ** 2 * y + x ** 3 * y ** 3 * z + x * y * z + 5,
         5 * (x ** 4) * (y ** 2) * (z ** 3) + ( x ** 3 * y * z ** 2) + x + 8 * y * z * x ** 2 + 3 * x ** 4 * y ** 4 * z ** 4
        ]
random_num = 0
current_func = funcs[random_num]

'''
 dense interpolation algorithm 'D'. Gives polynomial that fits all the given points.
 polynomial is such that,
 f(p_i) = m_i
'''


def dense_interpolation(p, m):
    x = symbols('x')
    f = m[0] * x ** 0
    # print type(f)
    q = x - p[0]

    for i in range(1, len(p)):
        temp = q.subs(x, p[i])
        qpi_inv = q.subs(x, p[i]) ** (-1)
        f += qpi_inv * q * (m[i] - f.subs(x, p[i]))
        q = (x - p[i]) * q
        # print 'f: ', f, '\nq: ', q, '\nqpi_inv: ', qpi_inv, '\n', 'temp: ', temp, '\n'
    if f == nan:
        return sympify(1)
    else:
        return expand(f)


'''ORACLE function 'F' which gives output of function at given points'''


def generate_function_F(starting_points):
    func = current_func
    return func.subs({x: starting_points[0], y: starting_points[1], z: starting_points[2]})


'''sparse interpolation algorithm 'S' '''


def sparse_interpolation(set_of_variables, starting_points, degree_bound, function_F):
    S = [sympify(0)]  # Set S = {(0)} and p0 = a0
    # print(S)
    linear_equations_list = []

    '''iterate through each variable'''
    for i in range(0, len(set_of_variables)):
        r = []
        X = set_of_variables[i]
        P = []

        for j in range(0, degree_bound):
            temp = random.randint(1, 10000)
            r.append(temp)

            linear_equations_list = []
            t = len(S)
            skeletal = sympify(0)
            G = symbols('g0:%d' % (t + 1))
            for k in range(0, t):
                skeletal += G[k] * S[k]
            for k in range(0, t):
                if i == 0:
                    generated_funtion_output = generate_function_F([r[j]] + starting_points[1:len(starting_points)])
                    # print oracle_out, ' oracle'
                    P.append(generated_funtion_output)
                    sub_list = [(set_of_variables[tem], A[tem]) for tem in range(i)]
                    linear_equations_list.append(skeletal.subs(sub_list) - generated_funtion_output)
                else:
                    A = random.sample(xrange(10000), i)
                    # print 'A: ', A
                    or_lis = A + [r[j]] + starting_points[i + 1:len(starting_points)]
                    generated_funtion_output = generate_function_F(or_lis)
                    sub_list = [(set_of_variables[tem], A[tem]) for tem in range(i)]
                    linear_equations_list.append(skeletal.subs(sub_list) - generated_funtion_output)
            if i != 0:  # S7
                sol1 = linsolve(linear_equations_list, G)
                solution = next(iter(sol1))
                # print 'solution: ', len(solution[0].free_symbols)
                # print solution, '\n'
                temp_p = sympify(0)

                for mon in range(len(S)):
                    temp_p += S[mon] * solution[mon]
                P.append(temp_p)

        '''
        For i = 0,  pass the list p directly to the dense_interpolation algorithm 'D'.
        Store the coefficients to the list 'S' and the resulting polynomial in function_F

        If i is not 0, then for each monomial is 'S', pass the corresponding coefficients from
        list 'P' and from function_F to the dense_interpolation algorithm 'D'.
        Merge each monomial to its corresponding result from 'D' and simplify
        '''
        if i == 0:
            function_F = sympify(
                dense_interpolation([starting_points[0]] + r, [function_F] + P).subs(symbols('x'), set_of_variables[0]))
            S = function_F.as_coefficients_dict().keys()
            print(S)
        else:
            P = [function_F] + P
            r = [starting_points[i]] + r
            function_F = sympify(0)
            # print P, 'P\n', r, 'r'
            for mon in range(len(S)):
                p_i_temp = []
                # m_i_temp = []

                for pol in range(len(P)):
                    p_i_temp += [P[pol].as_coefficients_dict()[S[mon]]]
                interp_temp = dense_interpolation(r, p_i_temp).subs(symbols('x'), set_of_variables[i])
                function_F += S[mon] * interp_temp
            function_F = expand(function_F)
            S = function_F.as_coefficients_dict().keys()
    return function_F


def draw_text(surface, font, text, position, color):
    # draw user-defined text in pygame graphics surface
    lable = font.render(text, 1, color)
    surface.blit(lable, position)


########################################################################
# Global settings:
SIZE = 500  # size of display window for single instance
STATUS_HEIGHT = 80  # height of status bar in display window
STATUS_HEIGHT2 = 30  # height of status bar within instance subwindows
STATUS_HEIGHT3 = 45  # height of status bar at the bottom (Github info)
DELIM_WIDTH = 5  # width of delimiter of the circles output
CITY_RADIUS = 5  # radius of circle representing city
FONTSIZE = 20  # font size for control section buttons
VERBOSE = False  # level of chattiness
SAVEPLOT = True  # save plot of tour length vs iteration (True) or only display it (False)
SLEEP = 0  # delay (in seconds) after plotting new configuration
N = 200  # initial number of cities
SEED = None  # random seed
VERSION = "1.0"  # version
COLORS = {"WHITE": (255, 255, 255), "RED": (255, 0, 0), "GREEN": (0, 255, 0), "BLUE": (0, 0, 255), "BLACK": (0, 0, 0),
          "YELLOW": (255, 255, 0),
          "LIGHT_BLUE": (0, 125, 227), "GREY1": (120, 120, 120), "GREY2": (224, 224, 224), "LIGHTBLUE": (102, 178, 255),
          "LIGHTRED": (255, 153, 153), "LIGHTYELLOW": (255, 255, 153), "PINK": (255, 51, 255), "DARKBLUE": (0, 0, 153),
          "LAVENDER": (204, 153, 255), "LIGHTGREEN": (153, 255, 204), "BROWN": (102, 51, 0), "OLIVE": (153, 153, 0),
          "DARKGREY": (105, 105, 105)}

########################################################################
# Initialisation:
pygame.init()
helv20 = pygame.font.SysFont("Helvetica", 20)
helv24 = pygame.font.SysFont("Helvetica", 24)
helv100 = pygame.font.SysFont("Helvetica", 40)

# start clock:
# set display surface for pygame:
SWIDTH = 2 * SIZE + DELIM_WIDTH
SHEIGHT = SIZE + STATUS_HEIGHT + STATUS_HEIGHT2 + STATUS_HEIGHT3 + 70

surface = pygame.display.set_mode((SWIDTH, SHEIGHT))
surface.set_alpha(None)
pygame.display.set_caption("CSE 426 SYMBOLIC COMPUTATION SEMESTER PROJECT")


class Button:

    def __init__(self, width, height, text, color, tcolor):
        self.width = width
        self.height = height
        self.text = text
        self.color = color
        self.tcolor = tcolor

    def SetText(self, text):
        self.text = text

    def PlaceButton(self, surface, x, y):
        self.x = x
        self.y = y
        surface = self.DrawButton(surface, x, y)
        surface = self.ButtonText(surface, x, y)

    def DrawButton(self, surface, x, y):
        pygame.draw.rect(surface, self.color, (x, y, self.width, self.height), 0)
        return surface

    def ButtonText(self, surface, x, y):
        font_size = int(self.width // len(self.text))
        font = pygame.font.SysFont("Arial", FONTSIZE)
        text = font.render(self.text, 1, self.tcolor)
        surface.blit(text, ((x + self.width / 2) - text.get_width() / 2, (y + self.height / 2) - text.get_height() / 2))
        return surface

    def IsPressed(self, mouse):
        return mouse[0] > self.x and \
               mouse[1] > self.y and \
               mouse[0] <= self.x + self.width and \
               mouse[1] < self.y + self.height


########################################################################
# Main loop:
def mainloop():
    global random_num, current_func, funcs
    interpolated_polynomial = ""
    button_run = Button(100, 30, "CALCULATE", COLORS["RED"], COLORS["BLACK"])
    image = pygame.image.load('back.jpg')
    image = pygame.transform.scale(image, (1050, 750))

    while True:
        # Event handler:
        for event in pygame.event.get():
            # pygame event handler
            if event.type == MOUSEBUTTONDOWN:
                if button_run.IsPressed(pygame.mouse.get_pos()):
                    random_num = random.randint(0, 2)
                    current_func = funcs[random_num]

                    interpolated_polynomial = "Generating And Calculating Interpolated Polynomial..."
                    text_surface = helv24.render(str(interpolated_polynomial), True, "GREEN")
                    surface.blit(text_surface, (260, 550))
                    pygame.display.flip()
                    time.sleep(1)
                    set_of_variables = symbols('x0:3')  # set of variables {X1, X2, ... ,Xv}
                    degree_bound = 4
                    starting_points = [1, 2, 3]  # {a1,a2, ... ,av}
                    function_F = sympify(generate_function_F(starting_points))

                    interpolated_polynomial = sparse_interpolation(set_of_variables, starting_points, degree_bound,
                                                                   function_F)
            elif event.type == KEYDOWN:
                # key is pressed
                if event.key in [K_ESCAPE, K_q]:
                    pygame.quit()
                    return


            surface.fill(COLORS["WHITE"])

            surface.blit(image, (0, 0))

            draw_text(surface, helv100, "CSE 426 SYMBOLIC COMPUTATION SEMESTER PROJECT", (35, 100), COLORS["RED"])
            draw_text(surface, helv100, "MULTIVARIATE POLYNOMIAL INTERPOLATION", (130, 160), COLORS["BLACK"])
            draw_text(surface, helv100, "Instructor: Asst.Prof.Zafeirakis Zafeirakopoulos", (135, 220),
                      COLORS["DARKGREY"])
            draw_text(surface, helv20, "Input Polynomial: ", (135 - 80, 300), COLORS["RED"])
            text_surface = helv24.render(str(current_func), True, (0, 0, 0))
            surface.blit(text_surface, (210, 300))

            draw_text(surface, helv20, "Interpolated Polynomial: ", (135 - 80, 370), COLORS["RED"])
            print(interpolated_polynomial)
            text_surface = helv24.render(str(interpolated_polynomial), True, (0, 0, 0))
            surface.blit(text_surface, (345 - 100, 370))

            button_run.PlaceButton(surface, 440, 600)
            draw_text(surface, helv20,
                      "Members: Mohammad Ashraf Yawar 161044123, Can Duyar 171044075, Canberk Arici 171044062",
                      (135, 690), COLORS["RED"])

            pygame.display.flip()


if __name__ == "__main__":
    mainloop()
