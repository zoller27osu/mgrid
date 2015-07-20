#!/bin/bash

ALTS=(J F)
LEGENDS=(Original Optimized)
COLORS=(CC0000 00CC00 0000CC CC00CC 00CCCC CCCC00 000000 CC8000 808080)
#echo "set linetype 2 lc rgb '#00C000' lw 1"

#egrep "pgo|pgn0|pgn1|pgn2" $1.log.$2 > pgo$2.dat
if [ -z "${1}" ]; then
    echo "Usage: $0 pbs_jobid" #[output_extension]"
    exit
fi
#output=${2:-postscript}

if [ ! -f "$1.out" ]; then
    echo "$1.out not found! Aborting."
    exit
fi

mkdir -p x2p_out
> $1.gp
echo "# Written by x2p.sh"                                              >> $1.gp
#echo "set terminal $output size 1920,1080 enhanced font 'Arial,20'"     >> $1.gp #pngcairo
#echo "set terminal $output size 1280,1024 font giant"                   >> $1.gp #png
echo "set terminal postscript eps enh solid lw 2 size 19.20,10.80 48"   >> $1.gp
echo "set output '$1.eps'"                                              >> $1.gp
#echo "set multiplot layout 1,1"                                         >> $1.gp
echo                                                                    >> $1.gp
echo "set style data lines"                                             >> $1.gp
#echo "set style line lw 3"                                             >> $1.gp
echo "set title 'Performance Improvment of Jacobi smoothing function'"  >> $1.gp
echo "set xlabel 'Number of gridpoints'"                                >> $1.gp
echo "set ylabel 'Execution time (s)'"                                  >> $1.gp
echo "set logscale xy 10"                                               >> $1.gp
for i in "${!COLORS[@]}"; do
    echo "set style line ($i+1) lc rgb '#${COLORS[$i]}' lw 2"           >> $1.gp
done
echo                                                                    >> $1.gp

#echo "Creating $1.dat"
first=true
for i in "${!ALTS[@]}"
do
    p=${ALTS[$i]}
    echo "Creating $1-$p.dat"
    > $1-$p.dat
    #echo "Adding case $p."
    #grep $p $1.log.$2 | awk -v n="$n" '$2 ~ n' > $p$2-$n.dat
    grep error$p $1.out                                         >> $1-$p.dat
    #echo ""                                                    >> $1.dat
    mv $1-$p.dat x2p_out/.    

    if [ $first ]
    then
        echo -n "plot '$1-$p.dat' using 10:4 title '${LEGENDS[$i]}' ls ($i+1)"  >> $1.gp
        first=""
    else
        echo ", \\"                                                             >> $1.gp
        echo -n "     '$1-$p.dat' using 10:4 title '${LEGENDS[$i]}' ls ($i+1)"  >> $1.gp
    fi
    #gnuplot -e "filename='$p$2.dat'" pgo.gp
done

# Remove the last (extraneous & error-causing) comma
#sed -i '$s/,$//' $1.gp
# 2x echo: 1 to end that last line, and 1 for a blank line.
echo                                                                    >> $1.gp
echo                                                                    >> $1.gp
echo "set terminal x11 enhanced font 'Arial,20' persist"                >> $1.gp
echo "replot"                                                           >> $1.gp
#echo "unset multiplot"                                                  >> $1.gp

mv $1.gp x2p_out/.
cd x2p_out

#gnuplot -e "filename='$1-$p.dat'" $1.gp
echo
echo "~~~ start $1.gp: ~~~"
cat $1.gp
echo "~~~~~ end $1.gp ~~~~"
echo
gnuplot $1.gp
#head -n -1 pgo$2.dat > temp.txt ; mv temp.txt pgo$2.dat
