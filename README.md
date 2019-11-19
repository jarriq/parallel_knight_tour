# Parallel Knight's Tour
Problema do passeio do cavalo usando OpenMPI e OpenMP em C++

# Instruções

Compilar:
mpic++ tour.cpp -fopenmp -o tour

Distribuir .o, onde m.txt é uma lista de IP's na rede local:
for ip in `cat m.txt`; do scp t10 $ip:;done

Rodar, onde N deve ser 5 <= N <= 8:
mpirun --machinefile m.txt -pernode tour <N>
