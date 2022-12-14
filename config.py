import os

def init():
    modelDirectoryExists = os.path.isdir('./modelmpi')
    if(modelDirectoryExists == True) : 
        return './modelmpi/result/'
    else:
        return './results/'