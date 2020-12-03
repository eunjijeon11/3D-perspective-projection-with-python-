import pygame
import os
import math
from matrix import matrix_multiplication

def main():
    os.environ["SDL_VIDEO_CENTERED"]='1'
    black, white, blue, red = (20, 20, 20), (230, 230, 230), (0, 154, 255), (237, 28, 38)
    width, height = 1920, 1080

    points = []

    run = False

    molecule = input("molecule(중심원자, 원자2, ...): ")
    mtype = find_type(molecule)

    if mtype != "wrong":
        if mtype == "사면체":
            points.append([[2], [0], [0]])
            points.append([[4], [0], [0]])
            points.append([[2],[2],[0]])
            points.append([[2],[0],[2]])
            points.append([[2+-2/3*(3**(1/2))],[-2/3*(3**(1/2))],[-2/3*(3**(1/2))]])
        elif mtype == "2원자분자":
            points.append([[3],[0],[0]])
            points.append([[1],[0],[0]])
        elif mtype == "직선형":
            points.append([[2],[0],[0]])
            points.append([[4],[0],[0]])
            points.append([[0],[0],[0]])
        elif mtype == "평면 삼각형":
            points.append([[2],[0],[0]])
            points.append([[4],[0],[0]])
            points.append([[1],[3**(1/2)],[0]])
            points.append([[1],[-3**(1/2)],[0]])
        elif mtype == "굽은형":
            points.append([[2],[0],[0]])
            points.append([[2],[2*math.sin(math.radians(52.25))],[-2*math.cos(math.radians(52.25))]])
            points.append([[2],[-2*math.sin(math.radians(52.25))],[-2*math.cos(math.radians(52.25))]])
        elif mtype == "삼각뿔형":
            a = 2
            b = 2/3*(3**(1/2))*a*math.cos(math.radians(36.5))
            h = (a**2 - b**2)**(1/2)
            points.append([[2],[0],[0]])
            points.append([[2 + b],[0],[-h]])
            points.append([[2 + -math.sin(math.radians(30))*b],[math.cos(math.radians(30))*b],[-h]])
            points.append([[2 + -math.sin(math.radians(30))*b],[-math.cos(math.radians(30))*b],[-h]])
        run = True
        pygame.init()
        pygame.display.set_caption("3D molecule projection")
        screen = pygame.display.set_mode((width, height))
        clock = pygame.time.Clock()
        fps = 60
        angle = 0
        cube_position = [width//2, height//2]
        scale = 600
        speed = 0.01

    while run:
        clock.tick(fps)
        screen.fill(white)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False

        index = 0
        projected_points = [j for j in range(len(points))]

        rotation_x = [[1, 0, 0],
                      [0, math.cos(angle), -math.sin(angle)],
                      [0, math.sin(angle), math.cos(angle)]]

        rotation_y = [[math.cos(angle), 0, -math.sin(angle)],
                      [0, 1, 0],
                      [math.sin(angle), 0, math.cos(angle)]]

        rotation_z = [[math.cos(angle), -math.sin(angle), 0],
                      [math.sin(angle), math.cos(angle), 0],
                      [0, 0 ,1]]

        for point in points:
            rotated_2d = matrix_multiplication(rotation_y, point)
            rotated_2d = matrix_multiplication(rotation_x, rotated_2d)
            rotated_2d = matrix_multiplication(rotation_z, rotated_2d)
            distance = 5
            z = 1/(distance - rotated_2d[2][0])
            projection_matrix = [[z, 0, 0],
                                [0, z, 0]]
            projected_2d = matrix_multiplication(projection_matrix, rotated_2d)

            x = int(projected_2d[0][0] * scale) + cube_position[0]
            y = int(projected_2d[1][0] * scale) + cube_position[1]
            projected_points[index] = [x, y]

            if index == 0:
                pygame.draw.circle(screen, red, (x, y), 50)
            else :
                pygame.draw.circle(screen, blue, (x, y), 50)
            index += 1

        #draw edges
        for i in range(len(points)-1):
            connect_point(0, i+1, projected_points, screen)

        angle += speed
        pygame.display.update()

    pygame.quit()

def find_type(molecule):
    non_genetic_pair = {"H":0, "Li":0, "Be":0, "B":0, "C":0,
                        "N":1, "O":2, "F":3}
    valence_emission_electromyography = {"H":1, "Li":1, "Be":2, "B":3, "C":4,
                        "N":5, "O":6, "F":7}
    
    element = molecule.split(',')
    mtype = "wrong molecule"
    if len(element) == 2 and element[0] == element[1] and (non_genetic_pair[element[0]] != 0 or element[0] == "H"):
        mtype = "2원자분자"
    elif len(element) == 2:
        mtype = "직선형"
    elif len(element) == 3:
        if non_genetic_pair[element[0]] == 0:
            mtype = "직선형"
        elif non_genetic_pair[element[0]] == 2:
            mtype = "굽은형"
    elif len(element) == 4:
        if non_genetic_pair[element[0]] == 1:
            mtype = "삼각뿔형"
        elif non_genetic_pair[element[0]] == 0:
            mtype = "평면 삼각형"
    elif len(element) == 5:
        mtype = "사면체"

    sum_elec = valence_emission_electromyography[element[0]]
    for i in range(1,len(element)):
        sum_elec += valence_emission_electromyography[element[i]]-2*non_genetic_pair[element[i]]
    if sum_elec != 8 and (element[0] not in ["H", "Li", "Be", "B"]):
        mtype = "wrong molecule"
        print("triggered")

    print(sum_elec, element[0], mtype)
    return mtype

def connect_point(i, j, k, screen, black=(20,20,20)):
    a = k[i]
    b = k[j]
    pygame.draw.line(screen, black, (a[0], a[1]), (b[0], b[1]), 2)

if __name__ == "__main__":
    main()
