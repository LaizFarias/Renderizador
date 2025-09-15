#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: <SEU NOME AQUI>
Disciplina: Computação Gráfica
Data: <DATA DE INÍCIO DA IMPLEMENTAÇÃO>
"""

import time         # Para operações com tempo
import gpu          # Simula os recursos de uma GPU
import math         # Funções matemáticas
import numpy as np  # Biblioteca do Numpy

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante
    # matrizes globais
    view_matrix = np.identity(4)
    projection_matrix = np.identity(4)
    transform_stack = []
    z_buffer = None

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far

    @staticmethod
    def polypoint2D(point, colors):
        """Função usada para renderizar Polypoint2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polypoint2D
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
        # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
        # você pode assumir inicialmente o desenho dos pontos com a cor emissiva (emissiveColor).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Polypoint2D : pontos = {0}".format(point)) # imprime no terminal pontos
        print("Polypoint2D : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo:
        coluna = colors.get("emissiveColor", [1.0, 1.0, 1.0])

        # ranges de RGB
        r = int(max(0, min(255, round(coluna[0] * 255))))
        g = int(max(0, min(255, round(coluna[1] * 255))))
        b = int(max(0, min(255, round(coluna[2] * 255))))

        i = 0 

        while i < (len(point)//2):
            x = int(round(point[2*i + 0]))
            y = int(round(point[2*i + 1]))
            if 0 <= x < GL.width and 0 <= y < GL.height:
                gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, [r, g, b])
            i += 1 
        # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)
        
    @staticmethod
    def polyline2D(lineSegments, colors):
        """Função usada para renderizar Polyline2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polyline2D
        # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
        # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
        # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
        # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
        # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
        # vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        print("Polyline2D : lineSegments = {0}".format(lineSegments)) # imprime no terminal
        print("Polyline2D : colors = {0}".format(colors)) # imprime no terminal as cores
        
        # Exemplo:
        coluna = colors.get("emissiveColor", [1.0, 1.0, 1.0])
        r = int(max(0, min(255, round(coluna[0] * 255))))
        g = int(max(0, min(255, round(coluna[1] * 255))))
        b = int(max(0, min(255, round(coluna[2] * 255))))

        if (len(lineSegments)//2) < 2:
            return

        i = 0 

        while i < (len(lineSegments)//2) - 1:
            x0 = int(round(lineSegments[2*i + 0]))
            y0 = int(round(lineSegments[2*i + 1]))
            x1 = int(round(lineSegments[2*(i+1) + 0]))
            y1 = int(round(lineSegments[2*(i+1) + 1]))

            dx = abs(x1 - x0)
            dy = abs(y1 - y0)

            if x0 < x1:
                sx = 1
            else:
                sx = -1

            if y0 < y1:
                sy = 1
            else:
                sy = -1

            err = dx - dy

            while True:
                if 0 <= x0 < GL.width and 0 <= y0 < GL.height:
                    gpu.GPU.draw_pixel([x0, y0], gpu.GPU.RGB8, [r, g, b])

                if x0 == x1 and y0 == y1:
                    break

                e2 = 2 * err

                if e2 > -dy:
                    err -= dy
                    x0 += sx

                if e2 < dx:
                    err += dx
                    y0 += sy

            i += 1
        # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)

    @staticmethod
    def circle2D(radius, colors):
        """Função usada para renderizar Circle2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Circle2D
        # Nessa função você receberá um valor de raio e deverá desenhar o contorno de
        # um círculo.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Circle2D
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

        print("Circle2D : radius = {0}".format(radius)) # imprime no terminal
        print("Circle2D : colors = {0}".format(colors)) # imprime no terminal as cores
        
        # Exemplo:
        pos_x = GL.width//2
        pos_y = GL.height//2
        gpu.GPU.draw_pixel([pos_x, pos_y], gpu.GPU.RGB8, [255, 0, 255])  # altera pixel (u, v, tipo, r, g, b)
        # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)


    @staticmethod
    def triangleSet2D(vertices, colors):
        """Função usada para renderizar TriangleSet2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#TriangleSet2D
        # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
        # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
        # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
        # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
        # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).
        print("TriangleSet2D : vertices = {0}".format(vertices)) # imprime no terminal
        print("TriangleSet2D : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo:
        coluna = colors.get("emissiveColor", [1.0, 1.0, 1.0])
        R = int(max(0, min(255, round(coluna[0] * 255))))
        G = int(max(0, min(255, round(coluna[1] * 255))))
        B = int(max(0, min(255, round(coluna[2] * 255))))

        if len(vertices) % 6 != 0:
            return

        for k in range(0, len(vertices), 6):
            ax, ay = vertices[k + 0], vertices[k + 1]
            bx, by = vertices[k + 2], vertices[k + 3]
            cx, cy = vertices[k + 4], vertices[k + 5]

            #box
            x_min_calc = math.floor(min(ax, bx, cx))
            x_max_calc = math.ceil (max(ax, bx, cx))
            y_min_calc = math.floor(min(ay, by, cy))
            y_max_calc = math.ceil (max(ay, by, cy))

            if x_min_calc < 0:
                x_min = 0
            else:
                x_min = x_min_calc
            if x_max_calc > GL.width - 1:
                x_max = GL.width - 1
            else:
                x_max = x_max_calc

            if y_min_calc < 0:
                y_min = 0
            else:
                y_min = y_min_calc
            if y_max_calc > GL.height - 1:
                y_max = GL.height - 1
            else:
                y_max = y_max_calc

            # area orientada
            area2 = (cx - ax) * (by - ay) - (cy - ay) * (bx - ax)
            if area2 == 0:
                continue
            elif area2 > 0:
                orient_positive = True
            else:
                orient_positive = False

            y = y_min
            while y <= y_max:
                x = x_min
                while x <= x_max:
                    w0 = (x - bx) * (cy - by) - (y - by) * (cx - bx)
                    w1 = (x - cx) * (ay - cy) - (y - cy) * (ax - cx)
                    w2 = (x - ax) * (by - ay) - (y - ay) * (bx - ax)

                    inside = False
                    if orient_positive:
                        if w0 >= 0 and w1 >= 0 and w2 >= 0:
                            inside = True
                    else:
                        if w0 <= 0 and w1 <= 0 and w2 <= 0:
                            inside = True

                    if inside:
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, [R, G, B])
                    x += 1
                y += 1

    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleSet
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, você pode assumir
        # inicialmente, para o TriangleSet, o desenho das linhas com a cor emissiva
        # (emissiveColor), conforme implementar novos materias você deverá suportar outros
        # tipos de cores.


        emissive = colors.get("emissiveColor", [1.0, 1.0, 1.0])
        R = int(emissive[0] * 255)
        G = int(emissive[1] * 255)
        B = int(emissive[2] * 255)
        rgb = [R, G, B]

        w, h = GL.width, GL.height
        if GL.z_buffer is None:
            GL.z_buffer = np.full((h, w), GL.far, dtype=float)

        i = 0
        while i < len(point):
            # --- Projeta os 3 vértices ---
            verts = []
            j = 0
            while j < 3:
                vx, vy, vz = point[i+3*j], point[i+3*j+1], point[i+3*j+2]
                v = np.array([vx, vy, vz, 1.0])
                if GL.transform_stack:
                    v = GL.transform_stack[-1] @ v
                v = GL.view_matrix @ v
                v = GL.projection_matrix @ v
                if v[3] != 0:
                    v /= v[3]
                sx = int((v[0] + 1) * w / 2)
                sy = int((1 - v[1]) * h / 2)
                sz = (v[2] + 1) / 2
                verts.append([sx, sy, sz])
                j += 1

            # ordena vértices pelo y (de cima para baixo)
            verts.sort(key=lambda v: v[1])
            x0, y0, z0 = verts[0]
            x1, y1, z1 = verts[1]
            x2, y2, z2 = verts[2]

            # --- varredura por linhas (scanline) ---
            y = y0
            while y <= y2:
                if y < y1:  # parte superior
                    if y1 != y0:
                        xa = x0 + (x1-x0)*(y-y0)/(y1-y0)
                        za = z0 + (z1-z0)*(y-y0)/(y1-y0)
                    else:
                        xa, za = x0, z0
                    if y2 != y0:
                        xb = x0 + (x2-x0)*(y-y0)/(y2-y0)
                        zb = z0 + (z2-z0)*(y-y0)/(y2-y0)
                    else:
                        xb, zb = x0, z0
                else:       # parte inferior
                    if y2 != y1:
                        xa = x1 + (x2-x1)*(y-y1)/(y2-y1)
                        za = z1 + (z2-z1)*(y-y1)/(y2-y1)
                    else:
                        xa, za = x1, z1
                    if y2 != y0:
                        xb = x0 + (x2-x0)*(y-y0)/(y2-y0)
                        zb = z0 + (z2-z0)*(y-y0)/(y2-y0)
                    else:
                        xb, zb = x0, z0

                if xa > xb:  # garante que xa é o da esquerda
                    xa, xb = xb, xa
                    za, zb = zb, za

                x = int(xa)
                while x <= int(xb):
                    z = za + (zb-za)*(x-xa)/(xb-xa) if xb != xa else za
                    if 0 <= x < w and 0 <= y < h:
                        if z < GL.z_buffer[y, x]:
                            gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, rgb)
                            GL.z_buffer[y, x] = z
                    x += 1
                y += 1

            i += 9

    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.


        #posição da câmera
        px, py, pz = position

        #aqui vamos assumir que a câmera sempre olha para -Z e "cima" é Y
        GL.view_matrix = np.identity(4)
        GL.view_matrix[0, 3] = -px
        GL.view_matrix[1, 3] = -py
        GL.view_matrix[2, 3] = -pz

        #matriz de projeção perspectiva 
        aspecto = GL.width / GL.height
        f = 1.0 / math.tan(fieldOfView / 2.0)

        P = np.zeros((4, 4))
        P[0, 0] = f / aspecto
        P[1, 1] = f
        P[2, 2] = (GL.far + GL.near) / (GL.near - GL.far)
        P[2, 3] = (2 * GL.far * GL.near) / (GL.near - GL.far)
        P[3, 2] = -1

        GL.projection_matrix = P

    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # ESSES NÃO SÃO OS VALORES DE QUATÉRNIOS AS CONTAS AINDA PRECISAM SER FEITAS.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo para depois potencialmente usar em outras chamadas. 
        # Quando começar a usar Transforms dentre de outros Transforms, mais a frente no curso
        # Você precisará usar alguma estrutura de dados pilha para organizar as matrizes.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        matriz = np.identity(4)

        if scale:
            sx, sy, sz = scale
            matriz[0, 0] *= sx
            matriz[1, 1] *= sy
            matriz[2, 2] *= sz

        # aplica rotação
        if rotation:
            eixo_x, eixo_y, eixo_z, ang = rotation
            c = math.cos(ang)
            s = math.sin(ang)

            if eixo_x == 1 and eixo_y == 0 and eixo_z == 0:
                rot = np.array([[1, 0, 0, 0],
                                [0, c,-s, 0],
                                [0, s, c, 0],
                                [0, 0, 0, 1]])
            elif eixo_x == 0 and eixo_y == 1 and eixo_z == 0:
                rot = np.array([[ c, 0, s, 0],
                                [ 0, 1, 0, 0],
                                [-s, 0, c, 0],
                                [ 0, 0, 0, 1]])
            elif eixo_x == 0 and eixo_y == 0 and eixo_z == 1:
                rot = np.array([[c,-s, 0, 0],
                                [s, c, 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]])
            else:
                # caso não seja eixo puro
                rot = np.identity(4) 

            matriz = rot @ matriz

        # aplica translação
        if translation:
            tx, ty, tz = translation
            trans = np.identity(4)
            trans[0, 3] = tx
            trans[1, 3] = ty
            trans[2, 3] = tz
            matriz = trans @ matriz

        # se já tem transformações anteriores acumula
        if GL.transform_stack:
            matriz = GL.transform_stack[-1] @ matriz

        # empilha resultado
        GL.transform_stack.append(matriz)

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        if len(GL.transform_stack) > 0:
            GL.transform_stack.pop()


    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleStripSet
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        base = 0
        s = 0
        while s < len(stripCount):
            total = stripCount[s]
            i = 0
            while i < total - 2:
                if i % 2 == 0:
                    tri = point[3*(base+i):3*(base+i)+3] + \
                          point[3*(base+i+1):3*(base+i+1)+3] + \
                          point[3*(base+i+2):3*(base+i+2)+3]
                else:
                    tri = point[3*(base+i):3*(base+i)+3] + \
                          point[3*(base+i+2):3*(base+i+2)+3] + \
                          point[3*(base+i+1):3*(base+i+1)+3]
                GL.triangleSet(tri, colors)
                i += 1
            base += total
            s += 1

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#IndexedTriangleStripSet
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        atual = []
        p = 0
        while p < len(index):
            idx = index[p]
            if idx == -1:
                j = 0
                while j < len(atual) - 2:
                    if j % 2 == 0:
                        tri = point[3*atual[j]:3*atual[j]+3] + \
                              point[3*atual[j+1]:3*atual[j+1]+3] + \
                              point[3*atual[j+2]:3*atual[j+2]+3]
                    else:
                        tri = point[3*atual[j]:3*atual[j]+3] + \
                              point[3*atual[j+2]:3*atual[j+2]+3] + \
                              point[3*atual[j+1]:3*atual[j+1]+3]
                    GL.triangleSet(tri, colors)
                    j += 1
                atual = []
            else:
                atual.append(idx)
            p += 1

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#IndexedFaceSet
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão não possui uma ordem oficial, mas em geral se o primeiro ponto com os dois
        # seguintes e depois este mesmo primeiro ponto com o terçeiro e quarto ponto. Por exemplo: numa
        # sequencia 0, 1, 2, 3, 4, -1 o primeiro triângulo será com os vértices 0, 1 e 2, depois serão
        # os vértices 0, 2 e 3, e depois 0, 3 e 4, e assim por diante, até chegar no final da lista.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        face = []
        k = 0
        while k < len(coordIndex):
            idx = coordIndex[k]
            if idx == -1:
                if len(face) >= 3:
                    v0 = face[0]
                    n = 1
                    while n < len(face) - 1:
                        v1, v2 = face[n], face[n+1]
                        tri = coord[3*v0:3*v0+3] + coord[3*v1:3*v1+3] + coord[3*v2:3*v2+3]
                        GL.triangleSet(tri, colors)
                        n += 1
                face = []
            else:
                face.append(idx)
            k += 1

    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone
        # A função cone é usada para desenhar cones na cena. O cone é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento bottomRadius especifica o
        # raio da base do cone e o argumento height especifica a altura do cone.
        # O cone é alinhado com o eixo Y local. O cone é fechado por padrão na base.
        # Para desenha esse cone você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cone : bottomRadius = {0}".format(bottomRadius)) # imprime no terminal o raio da base do cone
        print("Cone : height = {0}".format(height)) # imprime no terminal a altura do cone
        print("Cone : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder
        # A função cylinder é usada para desenhar cilindros na cena. O cilindro é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da base do cilindro e o argumento height especifica a altura do cilindro.
        # O cilindro é alinhado com o eixo Y local. O cilindro é fechado por padrão em ambas as extremidades.
        # Para desenha esse cilindro você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cylinder : radius = {0}".format(radius)) # imprime no terminal o raio do cilindro
        print("Cylinder : height = {0}".format(height)) # imprime no terminal a altura do cilindro
        print("Cylinder : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#DirectionalLight
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""
