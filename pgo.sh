#!/bin/bash

MULTIPLOT="" #"true"

GRAPHS=(pgo) #pgn0 pgn1 pgn2)

#egrep "pgo|pgn0|pgn1|pgn2" $1.log.$2 > pgo$2.dat

if [ -z "${2}" ]; then
    echo "Usage: $0 case_name number_of_nodes " #[graph_output_extension]"
    exit
fi
#output=${3:-png}

if [ ! -f "$1.log.$2" ] || ! grep -lq Done $1.log.$2 #grep -lq "Ctrl-C" $1.log.$2
then
    echo "'$1.log.$2' not found (or incomplete)! Running $1."
    if [ -n "${PE_ENV}" ]; then
        #printf %0.f\\n 15.4
        hours=$(echo $2 | awk '{print int((log(log($1)/log(2))/log(2))+0.5);}')
hours=1
        echo "$0: Calculated time needed: $hours hours."

        nekq $1 $2 $hours #> .temp
        #cat .temp
        #PBS_JOBID=$(tail -n1 <.temp)
        #rm .temp
        PBS_JOBID="FILL_THIS_IN"
        until [ -f "$1.log.$2" ]
        do
            sleep 30s
            echo "Waiting for execution of job $PBS_JOBID to begin."
        done
    else
        nekbmpi $1 $2
    fi
    until grep -l Done logfile
    do
        sleep 30s
        echo -n "Waiting for case $1 to finish: "
        awk 'NF{s=$0}END{print s}' logfile
    done
fi

for p in "${GRAPHS[@]}"
do
    datafile=$p$2.dat
    echo "Creating $datafile"

    OLD_METHOD=""
    grep $p $1.log.$2 > .tmp
  if [ "$OLD_METHOD"]; then

    > $datafile
    for ((n=1; n<$2; n++))
    do
        #grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n' > $p$2-$n.dat
        #grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n'              >> $datafile
        line=$(awk -v n="$n" '$1 ~ "^"n"$"' .tmp)
        if [ -z "$line" ]; then
            echo "Data for node $n missing, ending $datafile now."
            break
        fi
        echo "$line"                                            >> $datafile
        echo                                                    >> $datafile
    done
  
  else
    #TODO: look into consolidating these into 1 awk command
    awk '{if (NR>1 && save!=$1) print "";} {save=$1; print;}' .tmp > $datafile
  fi

    rm .tmp
    echo "Done creating $datafile."

    if [ "$MULTIPLOT" ]; then
        PLOTFILES=($p$2.gp)
    else
        PLOTFILES=($p$2.gp ${p}pw$2.gp)
    fi
    
    if [ "$MULTIPLOT" ]; then
#        length=2/"${#PLOTFILES[@]}"
#        width=1
#        if [ $length -gt 1 ]; then
#            width=2
#        fi
#        height=$length/$width
#        echo "lwh = $length, $width, $height"
        width=1
        height=2
    fi

    for i in "${!PLOTFILES[@]}"; do
        TERMINALS=("set term postscript eps enh color solid lw 2 size 19.20,10.80 48; set output '$1_$2_$i.eps'"
"set terminal x11 enh font 'Arial,20' persist")
    
        plotfile=${PLOTFILES[$i]}
        echo "Creating $plotfile"
        > $plotfile
        echo "# Written by pgo.sh"                                      >> $plotfile
        echo                                                            >> $plotfile
        echo "set style data lines"                                     >> $plotfile
        echo "set logscale xy 10"                                       >> $plotfile
        echo "set xlabel 'Message size in words'"                       >> $plotfile
        echo "set ylabel 'Time (s)'"                                    >> $plotfile
        echo                                                            >> $plotfile
        for t in "${TERMINALS[@]}"; do
            echo "$t"                                                   >> $plotfile
            #echo                                                        >> $plotfile
        #echo "set nokey"                                           >> $plotfile
        #echo "set size 1,1"                                        >> $plotfile
        #echo "set origin 0,0"                                      >> $plotfile
        #echo                                                       >> $plotfile
            if [ "$MULTIPLOT" ]; then
                echo "set multiplot layout $width,$height"              >> $plotfile
            fi
            #echo                                                       >> $plotfile
            #if [ "$MULTIPLOT" ] || [ $i -eq 0 ]; then
            if [ $i -eq 0 ]; then
                echo "set title '1/2 RT Ping-Pong Time ($p)'"           >> $plotfile
                echo "plot '$datafile' using 4:5 title '$datafile'"     >> $plotfile
            fi
            if [ "$MULTIPLOT" ]; then
                echo                                                    >> $plotfile
            fi
            if [ "$MULTIPLOT" ] || [ $i -eq 1 ]; then
                echo "set title '1/2 RT Ping-Pong Time Per Word ($p)'"  >> $plotfile
                echo "plot '$datafile' using 4:6 title '$datafile'"     >> $plotfile
            fi
            #echo                                                       >> $plotfile
            if [ "$MULTIPLOT" ]; then
                echo "unset multiplot"                                  >> $plotfile
            fi
            echo                                                        >> $plotfile
        done
        
        echo
        echo "~~~ start $plotfile: ~~~"
        cat $plotfile
        echo "~~~~~ end $plotfile ~~~~"
        echo
        gnuplot $plotfile

        #gs -dNOPAUSE -dEPSCrop -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile=$1_$2_$i.png $1_$2_$i.eps
    done
done
#head -n -1 pgo$2.dat > temp.txt ; mv temp.txt pgo$2.dat1
