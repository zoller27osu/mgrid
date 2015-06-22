#bw_ping_pong.gp
set title "Blue Waters: 1/2 pingpong"
set key autotitle columnhead
set logscale xy 10
set term postscript enh color
set terminal x11 persist
#set output "bw_ping_pong.ps"
set xlabel "Message Size in Words"
set ylabel "Time (seconds)"
set xrange [1:100000]
plot "bw_ping.out" u ($3):($4/1000000) w l
