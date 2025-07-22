!rm Movie/Conf1*
!rm Movie/Conf2*

#load nice palette of colors
load "colors.plt"

#ALPHA
a=1.2

set term pngcairo enh solid color
set size ratio 1
set xr [0:60]
set yr [0:60]

##ANIMATION WITH CELL SPINS
do for[i=0:500:10]{
file=sprintf("../Alpha%.1f/out_t%.1i",a,i)
filecom=sprintf("../Alpha%.1f/com_out",a)
print file
#file='../Alpha'.a.'/out_t'.i.''
#filecom='../Alpha'.a.'/com_out'
j=i/10
fileout=sprintf("Movie/Conf1_Alpha%.1f_t%05i.png",a,i)
set output fileout
p file u 1:2:3 palette pt 5 ps 1 not, filecom i j u 2:3 pt 7 ps 0.5 lt -1 not
}
!convert -delay 10 -loop 0 Movie/Conf1*.png animation1.gif

##ANIMATION WITH CELL ALPHAs
set palette defined (1"grey",2"red");
set cbrange [1:2]
do for[i=0:500:10]{
file=sprintf("../Alpha%.1f/out_t%.1i",a,i)
filecom=sprintf("../Alpha%.1f/com_out",a)
j=i/10
fileout=sprintf("Movie/Conf2_Alpha%05i_t%05i.png",a,i)
set output fileout
p file u 1:2:4 palette pt 5 ps 1 not, filecom i j u 2:3 pt 7 ps 0.5 lt -1 not
}
!convert -delay 10 -loop 0 Movie/Conf2*.png animation2.gif


#ANIMATION WITH CELL TYPES
set palette defined (0"grey",1"red",2"green");
set cbrange [0:2]
do for[i=0:500:10]{
file=sprintf("../Alpha%.1f/out_t%.1i",a,i)
filecom=sprintf("../Alpha%.1f/com_out",a)
j=i/10
fileout=sprintf("Movie/Conf3_Alpha%05i_t%05i.png",a,i)
set output fileout
p file u 1:2:5 palette pt 5 ps 1 not, filecom i j u 2:3 pt 7 ps 0.5 lt -1 not
}
!convert -delay 10 -loop 0 Movie/Conf3*.png animation3.gif

