from sympy import *

def diferencias_divididas (entrada, n):

  matriz_diferencias = []

  diferencias_orden_0 = []
  for i in range(n):
    diferencia = entrada[i][1]
    diferencias_orden_0.append(diferencia)

  matriz_diferencias.append(diferencias_orden_0)

  for j in range(n-1):
    vector_diferencias = []
    for i in range(n-j-1):
      diferencia = (matriz_diferencias[j][i+1]-matriz_diferencias[j][i])/(entrada[i+1+j][0]-entrada[i][0])
      vector_diferencias.append(diferencia)
    matriz_diferencias.append(vector_diferencias)

  return(matriz_diferencias)

def polinomio_newton (entrada, n):
  x = Symbol('x')
  matriz_diferencias_divididas = diferencias_divididas(entrada,n)

  Pn = 0
  for i in range(n):
    termino = matriz_diferencias_divididas[i][0]
    for j in range(i):
      termino = termino * (x-entrada[j][0])
    Pn = Pn + termino
  return (Pn)

def main ():
  
    ''' 
    Tabla de entrada segun el enunciado:
    ------------------------
    | x    | 1 | 2 | 3 | 5 |
    ------------------------
    | f(x) | 4 | 3 | 4 | 8 |
    ------------------------
    '''

    tabla_entrada = [[1,4],[2,3],[3,4],[5,8]]
    n = 4

    polinomio_newton = polinomio_newton(tabla_entrada, n)

    print(polinomio_newton)


main()