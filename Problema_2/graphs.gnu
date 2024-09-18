reset session

# PARÁMETROS ------------------------------------------------------------------
arch = "data/probabilities.txt"
arch_binder = "data/binder.txt"
arch_states = "data/densidades_em_2.txt"

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

set xlabel "Magnetization per spin m" font "Times New Roman, 22" offset 0,-2
set ylabel "Probability {/Symbol p}_m" font "Times New Roman, 22" offset -5,0

numero = system("read number; echo $number")

# 2x2 ----------------------------------------------------------
if (numero == 0) {
    n=4 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    plot[-1.5:1.5][0:0.5] arch u ($1/4):2 every :::0::0 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.5:1.5][0:0.5] arch u ($1/4):3 every :::0::0 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    #plot[-1.5:1.5][0:0.5] arch u ($1/4):4 every :::0::0 smooth freq w boxes lw 2 lc "#62a1db" notitle
}

# 4x4 ----------------------------------------------------------
if (numero == 1) {
    n=16 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    plot[-1.3:1.3][0:0.2] arch u ($1/16):2 every :::1::1 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.3:1.3][0:0.2] arch u ($1/16):3 every :::1::1 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    #plot[-1.3:1.3][0:0.2] arch u ($1/16):4 every :::1::1 smooth freq w boxes lw 2 lc "#62a1db" notitle
}

# 6x6 ----------------------------------------------------------
if (numero == 2) {
    n=36 #number of intervals
    max=1. #max value
    min=-1. #min value
    width=(max-min)/n #interval width
    set boxwidth width*0.9
    set style fill solid 0.7

    plot[-1.1:1.1][0:0.1] arch u ($1/36):2 every :::2::2 smooth freq w boxes lw 2 lc "#b01111" notitle
    #plot[-1.1:1.1][0:0.1] arch u ($1/36):3 every :::2::2 smooth freq w boxes lw 2 lc "#dd9f40" notitle
    #plot[-1.1:1.1][0:0.1] arch u ($1/36):4 every :::2::2 smooth freq w boxes lw 2 lc "#62a1db" notitle
}

# Binder ----------------------------------------------------------
if (numero == 3) {
    set xlabel "Temperature" font "Times New Roman, 22" offset 0,-2
    set ylabel "Binder cumulant" font "Times New Roman, 22" offset -5,0
    unset grid
    set object 1 rect from 1.95,0.9 to 2.45,0.98 dashtype 2 fs empty front

    set multiplot
    p [0:5][0.1:1.05] arch_binder u 1:2 every :::0::0 w lp lw 2 pt 7 lc "blue" t "2x2" at graph 0.84,0.92, \
        arch_binder u 1:2 every :::1::1 w lp lw 2 pt 7 lc "red" t "4x4" at graph 0.84,0.85, \
        arch_binder u 1:2 every :::2::2 w lp lw 2 pt 7 lc "#1bcc23" t "6x6" at graph 0.84,0.78
    
    set size 0.4
    set origin 0.17, 0.24
    unset object
    set grid
    set ytics 0.91, 0.03, 0.97
    set xlabel "Temperature" font "Times New Roman, 22" offset 0,17
    set ylabel "Binder cumulant" font "Times New Roman, 22" offset 41,0
    unset arrow
    set arrow from 2 / (log(1+sqrt(2))), 0.9 to 2 / (log(1+sqrt(2))), 0.98 dashtype 2 nohead

    p [1.95: 2.45][0.9:0.98] arch_binder u 1:2 every :::0::0 w lp lw 2 pt 7 lc "blue" notitle, \
        arch_binder u 1:2 every :::1::1 w lp lw 2 pt 7 lc "red" notitle, \
        arch_binder u 1:2 every :::2::2 w lp lw 2 pt 7 lc "#1bcc23" notitle
    unset multiplot
}

# Binder Límites ----------------------------------------------------------
if (numero == 4) {
    set xlabel "Temperature" font "Times New Roman, 22" offset 0,-2
    set ylabel "Binder cumulant" font "Times New Roman, 22" offset -5,0
    set ytics add ("0.1" 0.1 1, "0.3" 0.3 1, "0.5" 0.5 1, "0.7" 0.7 1)

    # Ajustes
    f1(x) = a1 + b1 / (x ** c1)
    a1 = 0.227746
    b1 = 1.79563
    c1 = 0.872942
    #fit f1(x) arch_binder u 1:2 every ::30:0::0 via a1,b1,c1

    f2(x) = a2 + b2 / (x ** c2)
    a2 = 0.0738571
    b2 = 5.69328
    c2 = 1.95846
    #fit f2(x) arch_binder u 1:2 every ::30:1::1 via a2,b2,c2

    f3(x) = a3 + b3 / (x ** c3)
    a3 = 0.0327684
    b3 = 6.08997
    c3 = 2.3255
    #fit f3(x) arch_binder u 1:2 every ::20:2::2 via a3,b3,c3
 
    # Gráfico
    p [][0:1.01] arch_binder u 1:2 every :::0::0 w l lw 2 lc "blue" t "2x2" at graph 0.84,0.92, \
        arch_binder u 1:2 every :::1::1 w l lw 2 lc "red" t "4x4" at graph 0.84,0.85, \
        arch_binder u 1:2 every :::2::2 w l lw 2 lc "#1bcc23" t "6x6" at graph 0.84,0.78
}

if (numero == 5) {
    # Set labels and title

    set xlabel "Energía (E)"
    set ylabel "Magnetización (M)"

    # Set palette and colorbox for color representation
    set palette rgbformulae 22,13,-31
    set colorbox
    #set xtics -10, 2, 10

    # Plot points, coloring them by N(E, M)
    #plot [-10:10][-6:6] arch_states using 1:2:3 every :::0::0 with points pointtype 5 pointsize 3 palette notitle
    #plot arch_states using 1:2:3 every :::1::1 with points pointtype 5 pointsize 3 palette notitle
    plot arch_states using 1:2:($3/1e8) every :::2::2 with points pointtype 5 pointsize 2 palette notitle
}