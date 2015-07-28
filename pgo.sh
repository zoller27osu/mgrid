#!/bin/bash
shopt -s extglob # enables extended parameter pattern matching

MULTIPLOT="true" #"true"

GRAPHS=(pgo) #pgn0 pgn1 pgn2)

# Usage info
show_help() {
#Usage: ${0##*/} [-hm] [-c CASE_NAME] [-n NUM_NODES]...
cat << EOF

Usage: ${0##*/} [-hm] [-c CASE_NAME] [-n NUM_NODES] [CASE_NAME] [NUM_NODES]...
Runs CASE_NAME, any derivative of examples/timing, if previous output is not found.
Processes and graphs the output to x11 and eps (importable PostScript).

    -h[elp]             display this help and exit
    -m[ultiplot]        plot the data as a Multiplot. (otherwise, multiple separate plots will be created)
    -c CASE_NAME        use output from the case with the given name (i.e. what you would supplied to makenek).
    -n NUM_NODES        use output from the named case, run with NUM_NODES processors.

CASE_NAME and NUM_NODES are parsed (in that order) at the end only if they were not previously specified with -c or -n respectively.

EOF
}
#if [ -z "${2}" ]; then
#    echo "Usage: $0 case_name number_of_nodes " #[graph_output_extension]"
#    exit 1
#fi
#output=${3:-png}

if [ $# -le 1 ]; then
    show_help
    echo
    read -p "$# argument(s) detected. Press [Enter] to continue anyway."
fi

# Initialize our own variables:
cn=""
np=""
MULTIPLOT=""

OPTIND=1 # Reset in case getopts is ever used previously
#while getopts "hmc:n:" opt; do
while [ $# -gt 0 ]; do # && [[ ."$1" = .--* ]]; do
    # Sets opt to the first non-hyphen character in $1
    opt="$1"
#    opt="${opt,,}" # should convert $opt to all lowercase...
    opt="${opt##+(-)}" #${1##+(-) #${1##-} #"$1"
    if [ ${#opt} -eq ${#1} ] && [[ "$1" != "?" ]]; then # first bit tests if $opt is shorter than $1
        # Unhyphenated, non-? argument: can't be handled by case statement.
#  echo "Unhyphenated argument \"${1}\" found, assuming hyphenated arguments are done..."
        break
    fi

    if [[ "$opt" =~ ^.*=.* ]]; then # opt contains a "="
        # the option is using the "-option=value" syntax
        OPTARG="${opt#*=}" # OPTARG gets opt minus the "=" and everything before it
#  echo "Detected \"-option=argument\" syntax, using argument \"$OPTARG\""
    else
        shift #expose next argument
        OPTARG=$1
    fi
    opt=${opt:0:1} # sets opts to its length 1 substring starting at 0 (i.e. its first character)
#  echo "option starts with \"$opt\""
    if [[ "$OPTARG" =~ ^-{1,2}.* ]]; then # OPTARG is not an argument but another hyphenated option
        #echo "WARNING: You may have left an argument blank. Double check your command."
        OPTARG=""
    fi

    case "$opt" in
        "m") #|"-multiplot")
            MULTIPLOT="true";;
        "?"|"h") #|"-help" )
            show_help;
            exit 1;; # exit 0;;
        "c") #|"-case_name")
            cn="$OPTARG"
            shift;;
        "n") #|"-num_nodes")
            np="$OPTARG"
            shift;;
        *)
            echo "ERROR: Invalid option: \""$opt"\"" >&2
            show_help >&2
            exit 1;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --. [$1 is now the first unhyphenated argument]

VAR_NAMES=("cn"        "np")
ARG_NAMES=("CASE_NAME" "NUM_NODES")
VAR_DEFAULTS=("timing" "8")
for i in "${!VAR_NAMES[@]}"; do
#    if [ -z "$1" ]; then
#        echo "All unhyphenated arguments have been parsed."
#        break
#    fi
    var=${VAR_NAMES[$i]}
    if [ -z "${!var}" ]; then
        arg=${ARG_NAMES[$i]}
        if [ -n "$1" ]; then
#  echo "Parsing unhyphenated argument \"$1\" as \"$arg\" (var name \"$var\")" #since it was not supplied as a hyphenated argument."
            eval "$var=$1"
            shift
        else
            default=${VAR_DEFAULTS[$i]}
            if [ -z "$default" ]; then
#  echo "ERROR: Mandatory argument \"$arg\" was not specified and has no default!" >&2
                exit 1
            else
#  echo "Set \"$arg\" to its default, \"$default\", since it was not specified by the user."
                eval "$var=$default"
            fi
        fi
    fi
done
if [ $# -gt 0 ]; then
    echo "WARNING: $# extra argument(s) were not parsed: \"$*\""
    read -p "    Press [Enter] to continue anyways."
fi

logfile="$cn.log.$np"

#  echo "logfile  =$logfile"
#  echo "MULTIPLOT=$MULTIPLOT"
#  exit 0

#
# END OF ARGUMENT PARSING
#

if [ ! -f "$logfile" ] || ! grep -lq Done $logfile #grep -lq "Ctrl-C" $logfile
then
    echo "'$logfile' not found (or incomplete)! Running $cn."
    if [ -n "${PE_ENV}" ]; then
        #printf %0.f\\n 15.4
#        hours=$(echo $np | awk '{print int((log(log($1)/log(2))/log(2))+0.5);}')
#hours=1
#        echo "$0: Calculated time needed: $hours hours."
#        nekq $cn $np $hours #> .tmp

        total_secs=$(echo $np | awk '{print int((log($1)/log(2))+0.5)*256}')
#  echo "total_secs = $total_secs"
          hours=$(echo $total_secs | awk '{print int($1 / 3600)}')
        minutes=$(echo $total_secs | awk '{print int(($1 % 3600)/60)}')
        seconds=$(echo $total_secs | awk '{print $1 % 60}')
        walltime=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)
        echo "Requested walltime = \"$walltime\""

        # cannot do via "x=$(nekq...); echo $x" or the output gets messed up
        nekq $cn $np $walltime > $0.tmp
        PBS_JOBID=$(tail -n1 <$0.tmp)
        cat $0.tmp
        echo "!!!!!!! End of nekq output !!!!!!!"
        echo
        rm $0.tmp
        until [ -f "$logfile" ]
        do
            echo "Waiting for execution of job \"$PBS_JOBID\" to begin."
            sleep 30s
        done
    else
        nekbmpi $cn $np
    fi
    until grep -l Done logfile
    do
        sleep 30s
        echo -n "Waiting for case $cn to finish: "
        awk 'NF{s=$0}END{print s}' logfile
    done
fi

for p in "${GRAPHS[@]}"
do
    datafile=$p$np.dat
    echo "Creating $datafile"

    grep $p $logfile > .tmp

#    OLD_METHOD="true"
#  if [ "$OLD_METHOD" ]; then
#    > $datafile
#    for ((n=1; n<$np; n++))
#    do
#        #data=$(awk -v n="$n" '$1 ~ "^"n"$"' .tmp)
#        data=$(awk -v n="$n" '$1 == n' .tmp)
#        if [ -z "$data" ]; then
#            echo "Data for node $n missing, ending $datafile now."
#            break
#        fi
##  echo "$data"
#        echo "$data"                                            >> $datafile
#        echo                                                    >> $datafile
#    done
#    head -n -1 $datafile > .tmp; mv .tmp $datafile
#  else
    #TODO: look into consolidating the grep line above and this one 1 awk command
    awk '{if (NR>1 && save!=$1) print "";} {save=$1; print;}' .tmp > $datafile
    last_node=$(awk 'END{print $1}' $datafile)
    max_node=`expr $np - 1`
    if [ $last_node -lt $max_node ]; then
        echo
        echo "    WARNING: Data only goes up to node $last_node of $max_node!"
        echo
    fi
    rm .tmp
#  fi

    echo "Done creating $datafile."

    if [ "$MULTIPLOT" ]; then
        PLOTFILES=($p$np.gp)
    else
        PLOTFILES=($p$np.gp ${p}pw$np.gp)
    fi

    if [ "$MULTIPLOT" ]; then
        width=1
        length=`expr 2 / ${#PLOTFILES[@]}`
        if [ $length -gt 2 ]; then
            width=2
        fi
        height=`expr $length / $width`
#  echo "lwh = $length, $width, $height"
    fi

    for i in "${!PLOTFILES[@]}"; do
        OUTPUT_FILE="${cn}_$np"
        if ! [ "$MULTIPLOT" ]; then
            OUTPUT_FILE="${OUTPUT_FILE}_$i"
        fi
        TERMINALS=("set term postscript eps enh color solid lw 2 size 19.20,10.80 48; set output '${OUTPUT_FILE}.eps'"
"set terminal x11 enh font 'Arial,20' persist")

        plotfile=${PLOTFILES[$i]}
        echo
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
            if [ "$MULTIPLOT" ]; then
                #echo                                                    >> $plotfile
                #echo "set nokey"                                        >> $plotfile
                #echo "set size 1,1"                                     >> $plotfile
                #echo "set origin 0,0"                                   >> $plotfile
                echo "set multiplot layout $width,$height"              >> $plotfile
                #echo                                                    >> $plotfile
            fi
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
            #echo                                                        >> $plotfile
            if [ "$MULTIPLOT" ]; then
                echo "unset multiplot"                                  >> $plotfile
            fi
            echo                                                        >> $plotfile
        done
        
        if [ $i -le 0 ]; then
            #echo
            echo "~~~ start $plotfile: ~~~"
            cat $plotfile
            echo "~~~~~ end $plotfile ~~~~"
            #echo
        fi
        gnuplot $plotfile

        echo "gs -dNOPAUSE -dEPSCrop -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile=${OUTPUT_FILE}.png ${OUTPUT_FILE}.eps"
    done
done
