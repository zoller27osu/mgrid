#!/bin/bash

GRAPHS=(pgo) #pgn0 pgn1 pgn2)

#egrep "pgo|pgn0|pgn1|pgn2" $1.log.$2 > pgo$2.dat

if [ -z "${2}" ]; then
    echo "Usage: $0 case_name number_of_nodes" #[graph_output_extension]"
    exit
fi
#output=${3:-png}

if [ ! -f "$1.log.$2" ] || grep -lq Done $1.log.$2 #grep -lq "Ctrl-C" $1.log.$2
then
    echo "'$1.log.$2' not found (or incomplete)! Running $1."
    if [ -n "${PE_ENV}" ]; then
        #printf %0.f\\n 15.4
        hours=$(echo $2 | awk '{print int((log(log($1)/log(2))/log(2))+0.5);}')
        echo "$0: Calculated time needed: $hours hours."

        nekq $1 $2 $hours #> .temp
        #cat .temp
        #PBS_JOBID=$(tail -n1 <.temp)
        #rm .temp
    else
        nekbmpi $1 $2
    fi
    until grep -l Done logfile
    do
        sleep 30s
        echo "Waiting for case to finish."
    done
fi

#length=${#GRAPHS[@]}
#width=1
#if [ length > 1 ]; then
#    width=2
#fi
#height=$length/$width

for p in "${GRAPHS[@]}"
do
    datafile=$p$2.dat
    echo "Creating $datafile"
    > $datafile
    for ((n=1; n<$2; n++))
    do
        #grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n' > $p$2-$n.dat
#        grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n'              >> $datafile
        echo                                                    >> $datafile
    done
    echo "Done creating $datafile."

    TERMINALS=("set term postscript eps enh solid lw 2 size 19.20,10.80 48; set output '$1_$2.eps'" "set terminal x11 enh font 'Arial,20' persist")
    #PLOTFILES=($p$2.gp ${p}pw$2.gp)

    plotfile=$p$2.gp
    echo "Creating $plotfile"
    > $plotfile
    echo "# Written by pgo.sh"                                              >> $plotfile
    echo                                                                    >> $plotfile
    echo "set style data lines"                                         >> $plotfile
    echo "set logscale xy 10"                                           >> $plotfile
    for t in "${TERMINALS[@]}"; do
        #echo "set term postscript eps enh solid lw 2 size 19.20,10.80 48"   >> $plotfile
        #echo "set output '$1_$2.eps'"                                       >> $plotfile
        echo "$t"                                                           >> $plotfile
        echo                                                                >> $plotfile
        #echo "set nokey"                                                    >> $plotfile
        #echo "set size 1,1"                                                 >> $plotfile
        #echo "set origin 0,0"                                               >> $plotfile
        #echo                                                                >> $plotfile
        echo "set multiplot layout 1,2"                                     >> $plotfile
        echo                                                                >> $plotfile
        echo "set title '1/2 RT Ping-Pong Time ($p)'"                       >> $plotfile
        echo "plot '$datafile' using 4:5 title '$datafile'"                 >> $plotfile
        echo                                                                >> $plotfile
        echo "set title '1/2 RT Ping-Pong Time Per Word ($p)'"              >> $plotfile
        echo "plot '$datafile' using 4:6 title '$datafile'"                 >> $plotfile
        echo                                                                >> $plotfile
        echo "unset multiplot"                                              >> $plotfile
        echo                                                                >> $plotfile
        #echo "set terminal x11 enhanced font 'Arial,20' persist"            >> $plotfile
    done
    echo
    echo "~~~ start $plotfile: ~~~"
    cat $plotfile
    echo "~~~~~ end $plotfile ~~~~"
    echo
    gnuplot $plotfile
done
#head -n -1 pgo$2.dat > temp.txt ; mv temp.txt pgo$2.dat
