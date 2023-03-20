import subprocess
import numpy as np


kp = open('times.txt', 'w')

commands = [
   2
]

def processResult(result: str):
    replaced1 = result.replace("\nComputation Done in ", "")
    replaced2 = replaced1.replace(" seconds!!\n", "")
    return replaced2

def exec(command):
    output = subprocess.run(command.split(" "), stdout=subprocess.PIPE)
    result = output.stdout

    return processResult(result.decode('utf-8'))


def main():
    exec("mpicc -O3 mainMPI.c modelMPI.c -o main")
    for command in commands:
        comm = "mpiexec -n {} ./main".format(command)
        results = []
        sum = 0
        for i in range(1):
            finalResult = exec(comm)
            parsed = float(finalResult)
            results.append(parsed)
        print(results)
        for k in results:
            sum += k
        mean = sum / 1
        stdd = 0
        for k in results:
            stdd += pow((k - mean), 2)
        stdd = stdd / 1
        to_write = "Mean for MPI with {} processes is: {} and stdd: {}\n".format(command, mean, stdd)
        kp.write(to_write)
main()