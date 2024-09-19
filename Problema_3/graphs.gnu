reset session

# PARÁMETROS ------------------------------------------------------------------
arch_e = "data/energy.txt"
arch_m = "data/magnetization.txt"

$Data <<EOD
0.5 0
0.75    0.00291608
1	    0.0232274
1.25	0.0837681
1.5	    0.19675
1.75	0.394346
2	    0.687991
2.25	0.983136
2.5	    1.00659
2.75	0.786446
3	    0.553771
3.25	0.398411
3.5	    0.294892
3.75	0.23055
4	    0.187003
4.25    0.155309
4.5	    0.132553
4.75	0.114256
5	    0.100026
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
    set xlabel "Temperature" font "Times New Roman, 22" offset 0,-2
    set ylabel "Mean energy per spin" font "Times New Roman, 22" offset -5,0
    plot [][-2.05:-0.4] arch_e u 1:2 every :::1::1 w lp lw 2 lc "blue" dashtype 2 pt 7 ps 1.5 notitle
}

# Cv --------------------------------------------------------------
if (numero == 1) {
    set angle degrees
    set samples 200

    colX     = 1
    colY     = 2
    j        = {0,1}                                 # imaginary unit
    a(dx,dy) = dx==0 && dy==0 ? NaN : atan2(dy,dx)   # angle of segment between two points
    L(dx,dy) = sqrt(dx**2 + dy**2)                   # length of segment
    r        = 0.3                                 # relative distance of ctrl points

    stats $Data u 0 nooutput   # get number of points+1
    N = STATS_records+1
    array P0[N]
    array PA[N]
    array PB[N]
    array P1[N]

    x1=x2=y1=y2=ap1=NaN
    stats $Data u (x0=x1, x1=x2, x2=column(colX), i=int($0)+1, \
                y0=y1, y1=y2, y2=column(colY), P0[i]=x0+j*y0, \
                dx1=x1-x0, dy1=y1-y0, d1=L(dx1,dy1), dx1n=dx1/d1, dy1n=dy1/d1, \
                dx2=x2-x1, dy2=y2-y1, d2=L(dx2,dy2), dx2n=dx2/d2, dy2n=dy2/d2, \
                a1=a(dx1,dy1), a2=a(dx2,dy2), a1=a1!=a1?a2:a1, \
                ap0=ap1, ap1=a(cos(a1)+cos(a2),sin(a1)+sin(a2)), \
                PA[i]=x0+d1*r*cos(ap0) + j*(y0+d1*r*sin(ap0)), \
                PB[i]=x1-d1*r*cos(ap1) + j*(y1-d1*r*sin(ap1)), P1[i]=x1+j*y1, 0) nooutput
    # add last segment
    P0[i+1] = x1+j*y1
    PA[i+1] = x1+d1*r*cos(ap1)+j*(y1+d1*r*sin(ap1))
    PB[i+1] = x2-d2*r*cos(a2) +j*(y2-d2*r*sin(a2))
    P1[i+1] = x2+j*y2

    # Cubic Bézier function with t[0:1] as parameter between two points
    # p0: start point, pa: 1st ctrl point, pb: 2nd ctrl point, p1: endpoint
    p(i,t) = t**3 * (  -P0[i] + 3*PA[i] - 3*PB[i] + P1[i]) + \
            t**2 * ( 3*P0[i] - 6*PA[i] + 3*PB[i]        ) + \
            t    * (-3*P0[i] + 3*PA[i]                  ) + P0[i]

    set xlabel "Temperature" font "Times New Roman, 22" offset 0,-2
    set ylabel "c_V per spin" font "Times New Roman, 22" offset -5,0

    plot $Data u 1:2 w p pt 7 ps 1.5 lc "blue" dt 3 notitle, \
        for [i=2:|P0|] [0:1] '+' u (real(p(i,$1))):(imag(p(i,$1))) w l lw 2 lc "blue" dashtype 2 notitle
}

# Magnetización ---------------------------------------------------
if (numero == 2) {
    set xlabel "Temperature" font "Times New Roman, 22" offset 0,-2
    set ylabel "Mean absolute magnetization" font "Times New Roman, 22" offset -5,0
    plot [0:5][0:1.05] arch_m u 1:2 every ::3:0::0 w lp lw 2 lc "blue" dashtype 2 pt 7 ps 1.2 t "4x4" at graph 0.85,0.92, \
        arch_m u 1:2 every ::3:1::1 w lp lw 2 lc "red" dashtype 2 pt 7 ps 1.2 t "8x8" at graph 0.85,0.85, \
        arch_m u 1:2 every ::3:2::2 w lp lw 2 lc "orange" dashtype 2 pt 7 ps 1.2 t "16x16" at graph 0.85,0.78, \
        arch_m u 1:2 every ::3:3::3 w lp lw 2 lc "#1bcc23" dashtype 2 pt 7 ps 1.2 t "32x32" at graph 0.85,0.71
}