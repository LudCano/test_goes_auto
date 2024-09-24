import numpy as np
import matplotlib.pyplot as plt
import os
print('ESTO ESTA CORRIENDO AHORA WAAAAAA')

print("Archivos en el directorio:", os.listdir())

# Datos para la gráfica
x = np.linspace(-10, 10, 100)  # 100 puntos entre -10 y 10
y = x ** 2  # y = x^2


# Crear la gráfica
plt.figure(figsize=(10, 6))
plt.plot(x, y, label='y = x^2', color='blue')
plt.title('Gráfica de y = x^2')
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(0, color='black', lw=0.5, ls='--')
plt.axvline(0, color='black', lw=0.5, ls='--')
plt.grid()
plt.legend()

# Guardar la gráfica
plt.savefig('output/graph.png')  # Guarda la gráfica en un archivo PNG
plt.close()  # Cierra la figura

print("Archivos en el directorio:", os.listdir())
