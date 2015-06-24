#!/bin/bash

#rm pgo.gp
#egrep "pgo|pgn0|pgn1|pgn2" timing.log.$1 > pgo$1.dat
#touch pgo$1.dat
if [ ! -f "timing.log.$1" ] || ! grep -lq Done timing.log.$1 #grep -lq "Ctrl-C" timing.log.$1
then
    echo "'timing.log.$1' not found! Running timing test."
    nekbmpi timing $1
    until grep -l Done logfile
    do
        sleep 2
    done
fi

for p in pgo #pgn0 pgn1 pgn2
do
    echo "Creating $p$1.dat"
    >$p$1.dat
    for ((n=1; n<$1; n++))
    do
        #grep $p timing.log.$1 | awk -v n="$n" '$1 ~ n' > $p$1-$n.dat
        grep $p timing.log.$1 | awk -v n="$n" '$1 ~ n' >> $p$1.dat
        echo "" >> $p$1.dat
    done
    
    gnuplot -e "filename='$p$1.dat'" pgo.gp
done
#head -n -1 pgo$1.dat > temp.txt ; mv temp.txt pgo$1.dat
