#!/bin/sh
for i in 5 15 25 35 45 55 65 75 85
do
    for j in 1 2 3 4 5 6 7 8 9
    do
        ./pbrt test$j.pbrt -theta_i $i -alpha $j -numrays 1000000
    done
done
