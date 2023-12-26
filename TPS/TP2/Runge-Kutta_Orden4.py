import matplotlib.pyplot as plt

def f(x, y):
    return 1.2 * x - 0.6 * x * y

def g(x, y):
    return 0.3 * x * y - 0.8 * y


def Runge_Kutta_4(xi, yi, h):

    m1 = f(xi, yi)
    k1 = g(xi, yi)

    m2 = f(xi + (h/2) * m1, yi + (h/2) * k1)
    k2 = g(xi + (h/2) * m1, yi + (h/2) * k1)

    m3 = f(xi + (h/2) * m2, yi + (h/2) * k2)
    k3 = g(xi + (h/2) * m2, yi + (h/2) * k2)

    m4 = f(xi + h * m3, yi + h * k3)
    k4 = g(xi + h * m3, yi + h * k3)


    xi1 = xi + (h/6) * (m1 + 2 * m2 + 2 * m3 + m4)
    yi1 = yi + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return xi1, yi1

def print_rk4(puntos):
    print(f"t_i \t\t | \t x_i \t\t | \t y_i")

    for punto in puntos:
        (t, x, y) = punto
        t = '%.5f' % t
        x = '%.5f' % x
        y = '%.5f' % y

        print(f"{t} \t | \t {x} \t | \t {y}")

def take_t(puntos):
    return puntos[0]

def take_x(puntos):
    return puntos[1]

def take_y(puntos):
    return puntos[2]

def graph_rk4(puntos_t, puntos_x, puntos_y):
    plt.plot(puntos_t, puntos_x, label = "x → presa", color='deepskyblue') 
    
    plt.plot(puntos_t, puntos_y, label = "y → depredador", color='palevioletred') 
    
    plt.xlabel('tiempo') 
    plt.ylabel('población') 
    plt.title('Modelo Depredador-Presa') 
    
    plt.legend(loc='upper right') 

    plt.grid(linewidth = 0.5)
    
    plt.show() 

def main():
    h = 0.1
    t = 0
    x = 2
    y = 1

    puntos = [[t,x,y]]

    while t <= 30:
        t = t + h
        x, y = Runge_Kutta_4(x, y, h)
        puntos.append([t,x,y])

    print(f"Simulacion a traves del Metodo de Runge-Kutta de Orden Cuatro: ")
    print_rk4(puntos)

    graph_rk4(list(map(take_t, puntos)), list(map(take_x, puntos)), list(map(take_y, puntos)))

main()