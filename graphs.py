import matplotlib.pyplot as plt

valores = [3219.73, 146.173, 85.87]
labels = ["Antes do ajuste", "Depois do ajuste", "Tratamento"]

plt.bar(labels, valores)
plt.ylabel("Valores")
plt.title("Gráfico de Barras")

plt.show()
