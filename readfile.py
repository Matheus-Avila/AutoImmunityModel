def ordenar_valores_arquivo(arquivo):
    lista_valores = []
    with open(arquivo, 'r') as f:
        for linha in f:
            alpha, beta, quality = linha.strip().split(',')
            dicionario = {"quality": float(quality), 'alpha': float(alpha), 'beta': float(beta)}
            lista_valores.append(dicionario)

    lista_valores.sort(key=lambda x: abs(x["quality"] - 70.63))

    return lista_valores
caminho_arquivo = "./cp.txt"
resultado = ordenar_valores_arquivo(caminho_arquivo)

for item in resultado:
    print(item)
