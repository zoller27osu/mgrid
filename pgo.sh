#!/bin/bash

#rm pgo.gp
#egrep "pgo|pgn0|pgn1|pgn2" $1.log.$2 > pgo$2.dat
#touch pgo$2.dat
if [ -z "${2}" ]; then
    echo "Usage: $0 case_name number_of_nodes"
    exit
fi

if [ ! -f "$1.log.$2" ] || ! grep -lq Done $1.log.$2 #grep -lq "Ctrl-C" $1.log.$2
then
    echo "'$1.log.$2' not found! Running $1."
    if [ -n "${PE_ENV}" ]; then
        #printf %0.f\\n 15.4
        hours=$(echo $2 | awk '{print int((log($1)/log(2))+0.5);}')
	echo "$0: Calculated time needed: $hours hours."
        nekq $1 $2 $hours
    else
        nekbmpi $1 $2
    fi
    until grep -l Done logfile
    do
        sleep 10
    done
fi

for p in pgo #pgn0 pgn1 pgn2
do
    echo "Creating $p$2.dat"
    >$p$2.dat
    for ((n=1; n<$2; n++))
    do
        #grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n' > $p$2-$n.dat
        grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n' >> $p$2.dat
        echo "" >> $p$2.dat
    done
    
    gnuplot -e "filename='$p$2.dat'" pgo.gp
done
#head -n -1 pgo$2.dat > temp.txt ; mv temp.txt pgo$2.dat
