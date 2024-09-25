reset session

# PARÁMETROS ------------------------------------------------------------------
arch_e = "data/energy.txt"
arch_m = "data/magnetization.txt"

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
    plot [][-2.05:-0.4] arch_e u 1:2 every :::1::1 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1 notitle
}

# Cv --------------------------------------------------------------
if (numero == 1) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "c_V per spin" font "Times New Roman, 22" offset -5,0

    plot [][-0.05:1.2] arch_e u 1:3 every :::1::1 w lp lw 1 pt 7 ps 1 lc "blue" dt 3 notitle
}

# Magnetización ---------------------------------------------------
if (numero == 2) {
    set xlabel "Temperature T" font "Times New Roman, 22" offset 0,-2
    set ylabel "Mean absolute magnetization m" font "Times New Roman, 22" offset -5,0
    plot [0:5][0:1.05] arch_m u 1:2 every ::3:0::0 w lp lw 1 lc "blue" dashtype 2 pt 7 ps 1 t "4x4" at graph 0.85,0.92, \
        arch_m u 1:2 every ::3:1::1 w lp lw 1 lc "red" dashtype 2 pt 7 ps 1 t "8x8" at graph 0.85,0.85, \
        arch_m u 1:2 every ::3:2::2 w lp lw 1 lc "orange" dashtype 2 pt 7 ps 1 t "16x16" at graph 0.85,0.78, \
        arch_m u 1:2 every ::3:3::3 w lp lw 1 lc "#1bcc23" dashtype 2 pt 7 ps 1 t "32x32" at graph 0.85,0.71
}