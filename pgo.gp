#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.2 patchlevel 6 
#    	last modified Sep 2009
#    	System: Linux 2.6.32-504.3.3.el6.x86_64
#    
#    	Copyright (C) 1986 - 1993, 1998, 2004, 2007 - 2009
#    	Thomas Williams, Colin Kelley and many others
#    
#    	Type `help` to access the on-line reference manual.
#    	The gnuplot FAQ is available from http://www.gnuplot.info/faq/
#    
#    	Send bug reports and suggestions to <http://sourceforge.net/projects/gnuplot>
#    
set terminal x11 0 persist
# set output

if (!exists("filename")) filename='plot.dat'
#n = 8
#from=system('tail -n +2 '.filename. '| cut -f 1 -d " " | sort | uniq')
#select_source(w) = sprintf('< awk ''{if ($1 == "%s") print }'' %s', w, filename)

set style data lines
set multiplot layout 1,2

set title '1/2 RT Ping-Pong Time'
set logscale xy 10
#plot for [f in from] select_source(f) using 4:5 #title f
plot filename using 4:5 title filename

set title '1/2 RT Ping-Pong Time Per Word'
#unset logscale y
#plot for [f in from] select_source(f) using 4:6 #title f
plot filename using 4:6 title filename

unset multiplot
