reset session

# PARÁMETROS ------------------------------------------------------------------
arch_direct = "data/data_direct_pi.txt"
arch_markov = "data/data_markov_pi.txt"

$Data <<EOD
10	0.78
20	0.789
30	0.78555
40	0.78765
50	0.785664
60	0.785275
70	0.785393
80  0.785406
100  0.78539816
EOD

$Data2 <<EOD
10 0.175
20	0.668
30	0.7883
40	0.783765
50	0.786039
60	0.785366
70	0.785447
80	0.785382
100 0.78539816
EOD

$Data3 <<EOD
0.15    4.24034e-01
0.2 3.49927e-01
0.25    2.20159e-01
0.3 1.92788e-01
0.35    1.9444e-01
0.4 1.51777e-01
0.45    1.05393e-01
0.5 7.10008e-02
0.55    5.39632e-02
0.6 8.84932e-02
0.65    7.64331e-02
0.7 6.49112e-02
0.75    4.42453e-02
0.8 4.61698e-02
0.85    1.07814e-01
0.9 9.08965e-02
0.95    1.25463e-01
1   5.32677e-02
1.05    8.89829e-02
1.1 6.36194e-02
1.15    4.95963e-02
1.2 5.51433e-02
1.25    5.53506e-02
1.3 5.23233e-02
1.35    9.71535e-02
1.4 1.05013e-01
1.45    1.04093e-01
1.5 7.50798e-02
1.55    6.09132e-02
1.6 8.93785e-02
1.65    1.15946e-01
1.7 9.90933e-02
1.75    7.54935e-02
1.8 6.85275e-02
1.85    1.60019e-01
1.9 8.39404e-02
1.95    5.76277e-02
2   8.0835e-02
2.05    1.90668e-01
2.1 1.53266e-01
2.15    1.86127e-01
2.2 1.25385e-01
2.25    1.64841e-01
2.3 1.0827e-01
2.35    3.49409e-01
2.4 1.12411e-01
2.45    1.23974e-01
2.5 2.22641e-01
2.55    2.52029e-01
2.6 2.39054e-01
2.65    2.66982e-01
2.7 2.74407e-01
2.75    2.95755e-01
2.8 2.53366e-01
2.85    2.94503e-01
2.9 2.80078e-01
2.95    3.07941e-01
3   2.09899e-01
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

set xlabel "Number of trials" font "Times New Roman, 22" offset 0,-2
set ylabel "N_{hits} / N" font "Times New Roman, 22" offset -5,0

numero = system("read number; echo $number")

# ESTIMACIÓN PI/4 SAMPLING DIRECTO: --------------------------------------------
if (numero == 0){
    set logscale x 10
    set arrow from 1,pi/4 to 10e8,pi/4 nohead dashtype 2 lw 2 lc "red"
    set xtics add ("10^1" 1e1 0, "10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "10^9" 1e9 0)

    p [1:1e9][0.775:0.795] arch_direct u 1:2 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle
}

# DESVIACIÓN ESTÁNDAR SAMPLING DIRECTO: ----------------------------------------
if (numero == 1){
    set xlabel "Number of trials" font "Times New Roman, 22" offset 0,-2
    set ylabel "MSE" font "Times New Roman, 22" offset -1,0

    # FIT
    f(x) = a*x**b
    a = 0.311693
    b = -1.05225
    #fit f(x) arch_direct u 1:4 every ::6:1::1 via a, b

    set logscale x 10
    set logscale y 10
    set xtics add ("10^1" 1e1 0, "10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "10^9" 1e9 0)
    set ytics add ("10^{-1}" 1e-1 0, "10^{-2}" 1e-2 0, "10^{-3}" 1e-3 0, "10^{-4}" 1e-4 0, "10^{-5}" 1e-5 0, "10^{-6}" 1e-6 0, "10^{-7}" 1e-7 0, "10^{-8}" 1e-8 0, "10^{-9}" 1e-9 0, "10^{-10}" 1e-10 0)

    p [1:1e9][1e-10:1e-1] arch_direct u 1:4 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle, f(x) lc "red" notitle
}

# SMOOTH PATH SAMPLING DIRECTO: ------------------------------------------------
if (numero == 2) {
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

    set multiplot 

    set logscale x 10
    set xtics add ("10^1" 1e1 0, "10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "10^9" 1e9 0)

    p [1:1e9][0.779:0.79] arch_direct u 1:2 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle

    unset logscale
    unset border
    set xtics format " "
    set label "{/Symbol p}/4 ≈ 0.7854" font "Times New Roman, 12" textcolor "red" at graph 0.025, 0.55
    set arrow from 0,pi/4 to 90,pi/4 nohead dashtype 2 lw 2 lc "red" back
    plot [0:90][0.779:0.79] $Data u 1:2 w p pt 7 lc "blue" dt 3 notitle, \
        for [i=2:|P0|] [0:1] '+' u (real(p(i,$1))):(imag(p(i,$1))) w l lc "blue" dashtype 2 notitle

    unset multiplot
    ### end of script
}

# ESTIMACIÓN PI/4 SAMPLING MARKOV: ---------------------------------------------
if (numero == 3){
    set multiplot
    unset grid
    set logscale x 10
    set arrow from 1,pi/4 to 10e8,pi/4 nohead dashtype 2 lw 2 lc "red"
    set xtics add ("10^1" 1e1 0, "10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "10^9" 1e9 0)

    p [1:1e9][0.1:0.85] arch_markov u 1:2 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle

    set size 0.55
    set origin 0.38, 0.22
    unset arrow
    set grid
    set xtics auto
    set xtics add ("10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "" 1e9 0)
    set ytics 0.782, 0.002, 0.79
    set arrow from 100,pi/4 to 10e8,pi/4 nohead dashtype 2 lw 2 lc "red"

    p [100:1e9][0.782:0.79] arch_markov u 1:2 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle

    unset multiplot
}

# DESVIACIÓN ESTÁNDAR SAMPLING MARKOV: -----------------------------------------
if (numero == 4){
    set xlabel "Number of trials" font "Times New Roman, 22" offset 0,-2
    set ylabel "MSE" font "Times New Roman, 22" offset -1,0

    # FIT
    f(x) = a*x**b
    a = 1.9985 
    b = -1.02532
    #fit f(x) arch_markov u 1:4 every ::6:1::1 via a, b

    set logscale x 10
    set logscale y 10
    set xtics add ("10^1" 1e1 0, "10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "10^9" 1e9 0)
    set ytics add ("10^{-1}" 1e-1 0, "10^{-2}" 1e-2 0, "10^{-3}" 1e-3 0, "10^{-4}" 1e-4 0, "10^{-5}" 1e-5 0, "10^{-6}" 1e-6 0, "10^{-7}" 1e-7 0, "10^{-8}" 1e-8 0, "10^{-9}" 1e-9 0)

    p [1:1e9][1e-9:1] arch_markov u 1:4 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle, f(x) lc "red" notitle
}

# SMOOTH PATH SAMPLING MARKOV: ---------------------------------------------
if (numero == 5){
    set angle degrees
    set samples 200

    colX     = 1
    colY     = 2
    j        = {0,1}                                 # imaginary unit
    a(dx,dy) = dx==0 && dy==0 ? NaN : atan2(dy,dx)   # angle of segment between two points
    L(dx,dy) = sqrt(dx**2 + dy**2)                   # length of segment
    r        = 0.3                                 # relative distance of ctrl points

    stats $Data2 u 0 nooutput   # get number of points+1
    N = STATS_records+1
    array P0[N]
    array PA[N]
    array PB[N]
    array P1[N]

    x1=x2=y1=y2=ap1=NaN
    stats $Data2 u (x0=x1, x1=x2, x2=column(colX), i=int($0)+1, \
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

    set multiplot 

    # Long range data
    unset grid
    set logscale x 10
    set xtics add ("10^1" 1e1 0, "10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "10^9" 1e9 0)

    p [1:1e9][0.1:0.85] arch_markov u 1:2 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle

    # Cubic Bézier Function Long
    unset logscale
    unset border
    set xtics format " "
    set label "{/Symbol p}/4 ≈ 0.7854" font "Times New Roman, 12" textcolor "red" at graph 0.05, 0.88
    set arrow from 0,pi/4 to 90,pi/4 nohead dashtype 2 lw 2 lc "red" back
    plot [0:90][0.1:0.85] $Data2 u 1:2 w p pt 7 lc "blue" dt 3 notitle, \
        for [i=2:|P0|] [0:1] '+' u (real(p(i,$1))):(imag(p(i,$1))) w l lc "blue" dashtype 2 notitle

    # Short range data
    set size 0.55
    set origin 0.38, 0.22
    unset arrow
    set grid
    set logscale
    set border 31
    set xtics auto
    set xtics add ("10^2" 1e2 0, "10^3" 1e3 0, "10^4" 1e4 0, "10^5" 1e5 0, "10^6" 1e6 0, "10^7" 1e7 0, "10^8" 1e8 0, "" 1e9 0)
    set ytics 0.782, 0.002, 0.79
    set ytics add ("0.784" 0.784 0, "0.786" 0.786 0, "0.788" 0.788 0, "0.790" 0.790 0)
    unset label

    p [100:1e9][0.782:0.79] arch_markov u 1:2 every :::1::1 w p ps 1.5 pt 7 lc "blue" notitle

    # Cubic Bézier Function Short
    set origin 0.38, 0.2205
    set ylabel " " font "Times New Roman, 22" offset -1,0
    unset logscale
    unset border
    set xtics format " "
    set arrow from 20,pi/4 to 90,pi/4 nohead dashtype 2 lw 2 lc "red" back
    plot [20:90][0.782:0.79] $Data2 u 1:2 w p pt 7 lc "blue" dt 3 notitle, \
        for [i=2:|P0|] [0:1] '+' u (real(p(i,$1))):(imag(p(i,$1))) w l lc "blue" dashtype 2 notitle

    unset multiplot
    ### end of script
}

# DESVIACIÓN ESTÁNDAR SAMPLING MARKOV DELTA: -----------------------------------
if (numero == 6){
    set angle degrees
    set samples 200

    colX     = 1
    colY     = 2
    j        = {0,1}                                 # imaginary unit
    a(dx,dy) = dx==0 && dy==0 ? NaN : atan2(dy,dx)   # angle of segment between two points
    L(dx,dy) = sqrt(dx**2 + dy**2)                   # length of segment
    r        = 0.2                                 # relative distance of ctrl points

    stats $Data3 u 0 nooutput   # get number of points+1
    N = STATS_records+1
    array P0[N]
    array PA[N]
    array PB[N]
    array P1[N]

    x1=x2=y1=y2=ap1=NaN
    stats $Data3 u (x0=x1, x1=x2, x2=column(colX), i=int($0)+1, \
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

    set multiplot
    set xlabel "δ" font "Times New Roman, 22" offset 0,-2
    set ylabel "MSE (x 10^{-5})" font "Times New Roman, 22" offset -4,0
    set y2label "Rejection rate" font "Times New Roman, 22" offset 4,0

    #set arrow from 1.225, 0 to 1.225, 0.27 dashtype 2 nohead
    #set arrow from 0.75, 0 to 0.75, 0.2 dashtype 2 nohead

    set y2tics 0, 0.1, 0.9
    set ytics nomirror
    set y2range [0:0.92]

    p [-0.1:3.2] arch_markov u 1:($4*10**5) every ::3:3::3 w p ps 1 pt 7 lc "blue" notitle, arch_markov u 1:2 every ::1:5::5 w lp ps 1 pt 7 lc "red" axes x1y2 notitle, arch_markov u 1:2 every ::15:5:15:5 w lp ps 1 pt 7 lc "#1bcc23" axes x1y2 notitle, arch_markov u 1:2 every ::23:5:26:5 w lp ps 1 pt 7 lc "#1bcc23" axes x1y2 notitle, $Data3 u 1:2 w p pt 7 lc "blue" dt 3 notitle,  \
        for [i=2:|P0|] [0:1] '+' u (real(p(i,$1))):(imag(p(i,$1))) w l lc "blue" dashtype 2 notitle, arch_markov u 1:($4*10**5) every ::15:3:15:3 w lp ps 1 pt 7 lc "#1bcc23" notitle, arch_markov u 1:($4*10**5) every ::23:3:26:3 w lp ps 1 pt 7 lc  "#1bcc23" notitle

    set size 0.854, 0.862
    set origin 0.107,0.088
    set border 2 linecolor "blue" lw 2
    set ylabel " "
    set y2label " "
    set xlabel " "
    unset y2tics
    unset xtics
    unset grid
    set ytics format " "
    unset arrow
    p [-0.1:3.2][0:0.45] -1 notitle

    set size 0.854, 0.862
    set origin 0.088,0.088
    set border 8 linecolor "red"
    set ylabel " "
    set y2label " "
    set xlabel " "
    unset ytics
    unset xtics
    set y2tics format " "
    unset arrow
    p [-0.1:3.2][0:0.45] -1 notitle

    unset multiplot
}