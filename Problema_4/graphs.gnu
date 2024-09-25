reset session

# PARÁMETROS ------------------------------------------------------------------
arch_e = "data/energy.txt"
arch_m = "data/magnetization.txt"
arc_m_counts = "data/magnetization_counts.txt"

$Data << EOD
0.5	-1.99997	0.0072887
0.6	-1.99995	0.0060367
0.7	-1.99989	0.00334602
0.8	-1.99959	0.00857384
0.9	-1.99882	0.0138695
1	-1.99698	0.0276457
1.1	-1.99378	0.0448503
1.2	-1.98844	0.0688029
1.3	-1.98015	0.103068
1.4	-1.96773	0.147041
1.5	-1.95089	0.197971
1.6	-1.92842	0.262276
1.7	-1.89764	0.345321
1.8	-1.85807	0.437286
1.9	-1.80842	0.560837
2	-1.74933	0.675597
2.1	-1.66967	0.82356
2.2	-1.58482	0.927544
2.3	-1.48565	1.01652
2.4	-1.38304	1.03167
2.5	-1.28383	1.00891
2.6	-1.18384	0.934822
2.7	-1.09578	0.841062
2.8	-1.01754	0.734236
2.9	-0.943536	0.641038
3	-0.8897	0.560278
3.1	-0.836257	0.481965
3.2	-0.789097	0.426539
3.3	-0.750292	0.373179
3.4	-0.714015	0.328365
3.5	-0.683342	0.294327
3.6	-0.654859	0.268608
3.7	-0.632904	0.243761
3.8	-0.608298	0.220822
3.9	-0.585954	0.204537
4	-0.563543	0.185907
4.1	-0.546561	0.173067
4.2	-0.532914	0.161877
4.3	-0.514345	0.151965
4.4	-0.498991	0.143341
4.5	-0.482937	0.132519
4.6	-0.472306	0.125234
4.7	-0.459871	0.119347
4.8	-0.448674	0.110765
4.9	-0.441371	0.106963
5	-0.429456	0.0996862
EOD

# AJUSTES Y GRÁFICOS ----------------------------------------------------------
set grid
set size 0.9
set origin 0.07,0.05
set border 31

unset label
unset arrow
unset logscale
set xtics auto
set ytics auto

set key font "Times New Roman, 20"
set tics font "Times New Roman, 22"
set xtics offset 0,-1
set ytics offset -1,0

numero = system("read number; echo $number")

# energy ----------------------------------------------------------
if (numero == 0) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "Mean energy per spin e" font "Times New Roman, 22" offset -5,0
    plot [][-2.05:-0.4] arch_e u 1:2 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1 notitle
}

# Cv --------------------------------------------------------------
if (numero == 1) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "c_V per spin" font "Times New Roman, 22" offset -5,0
    plot [][-0.05:1.2] arch_e u 1:3 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1 notitle
}

# Magnetización ------------------------------------------------
if (numero == 2) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "Mean absolute magnetization m" font "Times New Roman, 22" offset -5,0

    set ytics add ("0.1" 0.1 0, "0.3" 0.3 0, "0.5" 0.5 0, "0.7" 0.7 0, "0.9" 0.9 0)

    plot [0.5:5][0:1.05] arch_m u 1:2 every ::3:0::0 w lp lw 2 lc "blue" dashtype 2 pt 7 ps 1.2 t "6x6" at graph 0.85,0.92, \
        arch_m u 1:2 every ::3:1::1 w lp lw 1 lc "red" dashtype 2 pt 3 ps 1.2 t "16x16" at graph 0.85,0.85, \
        arch_m u 1:2 every ::3:2::2 w lp lw 1 lc "orange" dashtype 2 pt 5 ps 1.2 t "32x32" at graph 0.85,0.78, \
        arch_m u 1:2 every ::3:3::3 w lp lw 1 lc "#1bcc23" dashtype 2 pt 2 ps 1.2 t "64x64" at graph 0.85,0.71
}

# Binder Cumulant ------------------------------------------------
if (numero == 3) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "Binder Cumulant B(T)" font "Times New Roman, 22" offset -5,0

    set ytics add ("0.1" 0.1 0, "0.3" 0.3 0, "0.5" 0.5 0, "0.7" 0.7 0, "0.9" 0.9 0, "-0.1" -0.1 0)

    plot [0.5:5][-0.05:1.05] arch_m u 1:3 every ::3:0::0 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1.2 t "6x6" at graph 0.85,0.92, \
        arch_m u 1:3 every ::3:1::1 w lp lw 1 lc "red" dashtype 2 pt 3 ps 1.2 t "16x16" at graph 0.85,0.85, \
        arch_m u 1:3 every ::3:2::2 w lp lw 1 lc "orange" dashtype 2 pt 5 ps 1.2 t "32x32" at graph 0.85,0.78, \
        arch_m u 1:3 every ::3:3::3 w lp lw 1 lc "#1bcc23" dashtype 2 pt 2 ps 1.2 t "64x64" at graph 0.85,0.71
}

if (numero == 4) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "Binder Cumulant B(T)" font "Times New Roman, 22" offset -5,0

    set ytics add ("0.1" 0.1 0, "0.3" 0.3 0, "0.5" 0.5 0, "0.7" 0.7 0, "0.9" 0.9 0, "-0.1" -0.1 0)
    set arrow from 2 / (log(1+sqrt(2))), 0.5 to 2 / (log(1+sqrt(2))), 1.05 dashtype 2 nohead
    set label "T_c ≈ 2.269" font "Times New Roman, 18" at graph 0.4, 0.22

    plot [2.15:2.35][0.5:1.05] arch_m u 1:3 every ::3:0::0 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1.2 t "6x6" at graph 0.2,0.69, \
        arch_m u 1:3 every ::3:1::1 w lp lw 1 lc "red" dashtype 2 pt 3 ps 1.2 t "16x16" at graph 0.2,0.62, \
        arch_m u 1:3 every ::3:2::2 w lp lw 1 lc "orange" dashtype 2 pt 5 ps 1.2 t "32x32" at graph 0.2,0.55, \
        arch_m u 1:3 every ::3:3::3 w lp lw 1 lc "#1bcc23" dashtype 2 pt 2 ps 1.2 t "64x64" at graph 0.2,0.48
}

if (numero == 5) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "Binder cumulant B(T)" font "Times New Roman, 22" offset -5,0
    set object 1 rect from 2.15,0.8 to 2.35,0.98 dashtype 2 fs empty front
    unset grid

    set multiplot
    plot [0.5:6.5][-0.05:1.05] arch_m u 1:3 every ::3:0::0 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1.2 t "6x6" at graph 0.15,0.29, \
        arch_m u 1:3 every ::3:1::1 w lp lw 1 lc "red" dashtype 2 pt 3 ps 1.2 t "16x16" at graph 0.15,0.22, \
        arch_m u 1:3 every ::3:2::2 w lp lw 1 lc "orange" dashtype 2 pt 5 ps 1.2 t "32x32" at graph 0.15,0.15, \
        arch_m u 1:3 every ::3:3::3 w lp lw 1 lc "#1bcc23" dashtype 2 pt 2 ps 1.2 t "64x64" at graph 0.15,0.08

    set size 0.4
    set origin 0.51, 0.48
    unset object
    set grid
    set xtics 2.15, 0.05 ,2.35
    set xlabel "T" font "Times New Roman, 22" offset 0,17
    set ylabel "B(T)" font "Times New Roman, 22" offset 41,0
    unset arrow
    set arrow from 2 / (log(1+sqrt(2))), 0.8 to 2 / (log(1+sqrt(2))), 1.01 dashtype 2 nohead
    set label "T_c ≈ 2.269" font "Times New Roman, 18" at graph 0.15, 0.22

    plot [2.15:2.35][0.8:1.01] arch_m u 1:3 every ::3:0::0 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1.2 notitle, \
        arch_m u 1:3 every ::3:1::1 w lp lw 1 lc "red" dashtype 2 pt 3 ps 1.2 notitle, \
        arch_m u 1:3 every ::3:2::2 w lp lw 1 lc "orange" dashtype 2 pt 5 ps 1.2 notitle, \
        arch_m u 1:3 every ::3:3::3 w lp lw 1 lc "#1bcc23" dashtype 2 pt 2 ps 1.2 notitle

    unset multiplot
}

# Histograma 6x6
if (numero == 6) {
    set xlabel "Magnetization per spin m" font "Times New Roman, 22" offset 0,-2
    set ylabel "Counts/n_samples" font "Times New Roman, 22" offset -5,0
    n=36 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    #plot[-1.5:1.5][0:0.5] arc_m_counts u ($1/n):($2/10**5) every :::0::0 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.5:1.5][0:0.2] arc_m_counts u ($1/n):($3/10**5) every :::0::0 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    plot[-1.5:1.5][0:0.1] arc_m_counts u ($1/n):($4/10**5) every :::0::0 smooth freq w boxes lw 2 lc "#62a1db" notitle
}

if (numero == 7) {
    set xlabel "Magnetization per spin m" font "Times New Roman, 22" offset 0,-2
    set ylabel "Counts / nº samples" font "Times New Roman, 22" offset -5,0
    n=16*16 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    #plot[-1.5:1.5][0:0.5] arc_m_counts u ($1/n):($2/10**5) every :::1::1 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.5:1.5][0:0.02] arc_m_counts u ($1/n):($3/10**5) every :::1::1 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    #plot[-1.5:1.5] arc_m_counts u ($1/n):($4/10**5) every :::1::1 smooth freq w boxes lw 2 lc "#62a1db" notitle
}

if (numero == 8) {
    set xlabel "Magnetization per spin m" font "Times New Roman, 22" offset 0,-2
    set ylabel "Counts / nº samples" font "Times New Roman, 22" offset -5,0
    n=32*32 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    #plot[-1.5:1.5][0:0.5] arc_m_counts u ($1/n):($2/10**5) every :::2::2 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.5:1.5][0:0.005] arc_m_counts u ($1/n):($3/10**5) every :::2::2 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    plot[-1.5:1.5] arc_m_counts u ($1/n):($4/10**5) every :::2::2 smooth freq w boxes lw 2 lc "#62a1db" notitle
}

if (numero == 9) {
    set xlabel "Magnetization per spin m" font "Times New Roman, 22" offset 0,-2
    set ylabel "Counts" font "Times New Roman, 22" offset -5,0
    n=64*64 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    #plot[-1.5:1.5] arc_m_counts u ($1/n):($2/10**5) every :::3::3 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.5:1.5] arc_m_counts u ($1/n):$3 every :::3::3 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    plot[-1.5:1.5] arc_m_counts u ($1/n):4 every :::3::3 smooth freq w boxes lw 2 lc "#62a1db" notitle
}