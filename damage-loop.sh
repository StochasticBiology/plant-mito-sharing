gcc -o0 damage-loop.c -lm -o damage-loop.ce -g
./damage-loop.ce 0 0 > out-0-0.txt &
./damage-loop.ce 1 0 > out-1-0.txt &
./damage-loop.ce 2 0 > out-2-0.txt &
./damage-loop.ce 3 0 > out-3-0.txt &
./damage-loop.ce 4 0 > out-4-0.txt &
./damage-loop.ce 5 0 > out-5-0.txt &

./damage-loop.ce 0 1 > out-0-1.txt &
./damage-loop.ce 1 1 > out-1-1.txt &
./damage-loop.ce 2 1 > out-2-1.txt &
./damage-loop.ce 3 1 > out-3-1.txt &
./damage-loop.ce 4 1 > out-4-1.txt &
./damage-loop.ce 5 1 > out-5-1.txt &
