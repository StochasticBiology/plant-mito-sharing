gcc -o3 damage-loop.c -lm -o damage-loop.ce
./damage-loop.ce 0 > out-0.txt &
./damage-loop.ce 1 > out-1.txt &
./damage-loop.ce 2 > out-2.txt &
