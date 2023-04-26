import os

# Cria o diretório 'modelmpi' na pasta atual
os.mkdir('modelmpi')

# Cria o diretório 'result' dentro do diretório 'modelmpi'
os.mkdir(os.path.join('modelmpi', 'result'))

# Cria os diretórios 'ant', 'da', 'dc', 'matrix', 'mic', 'odc' e 'tke' dentro do diretório 'result'
for subdir in ['ant', 'da', 'dc', 'matrix', 'mic', 'odc', 'tke']:
    os.mkdir(os.path.join('modelmpi', 'result', subdir))
