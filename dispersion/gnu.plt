reset
set terminal eps font "Helvetica,22" enhanced size 8,8
#set size square

ps="0.440994 0.781548"
From="1 3"
To="2 4"
phasesF="'molecular I' 'molecular II'"
phasesT="'molecular II' 'quasiatomic'"
lvlOPT="0 1"
set terminal eps font "Helvetica,22" enhanced size 8,8
set output "dispersion.eps"

#
# HELP WITH LABELING MULTIPLOT!
#

offX=0.07	# left/right  margin
offY=0.07	# top /bottom margin
loffX=0.018	# move label along X
loffY=0.015	# move label along Y
dx=0.08		# space between plots in X
dy=0.08		# space between plots in Y

# names of plots
set label 1001 "(a)" at screen 0.0+offX+loffX,1.0-offY-loffY center front
set label 1002 "(b)" at screen 0.5+dx*0.5+loffX,1.0-offY-loffY center front
set label 1003 "(c)" at screen 0.0+offX+loffX,0.5-dy*0.5-loffY center front
set label 1004 "(d)" at screen 0.5+dx*0.5+loffX,0.5-dy*0.5-loffY center front

# draw grid
#do for [i=0:10]{
#set arrow from screen 0,0.1*i to screen 1,0.1*i nohead lt 3 lc -1
#set arrow from screen 0.1*i,0 to screen 0.1*i,1 nohead lt 3 lc -1
#}
#
## draw all margins
#do for [i=0:10]{
#set arrow from screen 0.1*i,0.5 to screen 0.1*i,0.5-dy*0.5 heads filled
#set arrow from screen 0.1*i,0.5 to screen 0.1*i,0.5+dy*0.5 heads filled
#set arrow from screen 0.1*i,0.0 to screen 0.1*i,offY heads filled
#set arrow from screen 0.1*i,1.0 to screen 0.1*i,1.0-offY heads filled
#set arrow from screen 0.5,0.1*i to screen 0.5-dx*0.5,0.1*i heads filled
#set arrow from screen 0.5,0.1*i to screen 0.5+dx*0.5,0.1*i heads filled
#set arrow from screen 0.0,0.1*i to screen offX,0.1*i heads filled
#set arrow from screen 1.0,0.1*i to screen 1.0-offX,0.1*i heads filled
#}

set multiplot layout 2,2 \
              margins screen 0.0+offX,1.0-offX,0.0+offY,1.0-offY \
              spacing screen dx,dy

#
# END WITH unset multiplot
#


set style arrow 1 nohead lt 1 lc -1 lw 3 dt '.'
set style arrow 2 nohead lt 1 lc 14 lw 3 dt '-'
set style arrow 3 nohead lt 1 lc  7 lw 3 dt '-'
set style arrow 4 heads filled size screen 0.01,15,45 fixed lt 1 lw 5 lc 14
set style arrow 5 heads filled size screen 0.01,15,45 fixed lt 1 lw 5 lc 7

ps="0.440994 0.781548"
From="1 3"
To="2 4"
Ysmin="-2 -2"
Ysmax="2 2"

FrE1Min=" -1.47520021402171	-1.10192402620007"
FrE1Max=" -0.842447178566	-0.742010590791112"
FrE2Min=" 0.479251073680972	0.831524197980498"
FrE2Max=" 1.31417583311681	1.04277959294546"
ToE1Min=" -1.34844880539873	-1.74022391403419"
ToE1Max=" -1.09714846376527	0.950128102528751"
ToE2Min=" 1.18776717308646	-1.60243608009776"
ToE2Max=" 1.42390538279933	0.986134855439932"

phasesF="'molecular I' 'molecular II'"
phasesT="'molecular II' 'quasiatomic'"
lvlOPT="0 1"
do for [i=1:2]{

p=word(ps,i)
f=word(From,i)
t=word(To,i)
ymin=word(Ysmin,i)
ymax=word(Ysmax,i)
nameofphaseF=word(phasesF,i)
nameofphaseT=word(phasesT,i)
lvl=word(lvlOPT,i)

ph1e10=word(FrE1Min,i)
ph1e11=word(FrE1Max,i)
ph1e20=word(FrE2Min,i)
ph1e21=word(FrE2Max,i)
ph2e10=word(ToE1Min,i)
ph2e11=word(ToE1Max,i)
ph2e20=word(ToE2Min,i)
ph2e21=word(ToE2Max,i)

set yrange[ymin:ymax]
set y2range[ymin:ymax]
unset y2tics
set ytics 5 offset 0,0 font "Helvetica,16"
set ytics add ("{/Symbol e}_a" 0)
set ytics add ("+1 Ry" 1)
set ytics add ("+2 Ry" 2)
set ytics add ("+3 Ry" 3)
set ytics add ("-1 Ry" -1)
set ytics add ("-2 Ry" -2)
set ytics add ("-3 Ry" -3)
set xrange[0:3000]
set xtics 1000 offset 0,0.5 font "Helvetica,16"
set xtics add ("{/Symbol G}" 0)
set xtics add ("X" 1000)
set xtics add ("M" 2000)
set xtics add ("{/Symbol G}" 3000)
unset xlabel
unset ylabel

set arrow 1 from 1000,ymin to 1000,ymax as 1
set arrow 2 from 2000,ymin to 2000,ymax as 1


set arrow 3 from   0,ph1e10 to 3000,ph1e10 as 2
set arrow 4 from   0,ph1e11 to 3000,ph1e11 as 2
set arrow 5 from   0,ph1e20 to 3000,ph1e20 as 3
set arrow 6 from   0,ph1e21 to 3000,ph1e21 as 3
set arrow 7 from 2000,ph1e10 to 2000,ph1e11 as 4
set arrow 8 from 1000,ph1e20 to 1000,ph1e21 as 5

set label 10 sprintf("W =%1.2f Ry",ph1e11-ph1e10) at  1500,(ph1e10-0.2) center font "Helvetica,16" tc lt 14
set label 11 sprintf("W'=%1.2f Ry",ph1e21-ph1e20) at  1500,(ph1e21+0.2) center font "Helvetica,16" tc lt 7

	mu=0.0 #system(sprintf("awk '{if($1==%s){print $2}}' phase%s_mu.dat",p,f))+.0
	set arrow 9 from 0,mu to 3000,mu nohead lt 1 lc 0 lw 5
	set label 5 "{/Symbol m}={/Symbol e}_a" at 1500,-0.2 center font "Helvetica,16"

	#set title offset 0,-0.8 sprintf("From: %s",nameofphaseF)
plot	sprintf("band%d.dat",i-1) u 0:4 w l lw 5 lc 14 notitle,\
	sprintf("band%d.dat",i-1) u 0:5 w l lw 5 lc 7 notitle
#plot	"band0.dat" u 0:4 w l lw 5 lc 14 notitle,\
#	"band0.dat" u 0:5 w l lw 5 lc 7 notitle

unset arrow 9
unset label 1001
unset label 1002
unset label 1003
unset label 1004


set yrange[ymin:ymax]
set y2range[ymin:ymax]
unset ytics
set y2tics 5 offset 0,0 font "Helvetica,16"
set y2tics add ("{/Symbol e}_a" 0)
set y2tics add ("+1 Ry" 1)
set y2tics add ("+2 Ry" 2)
set y2tics add ("+3 Ry" 3)
set y2tics add ("-1 Ry" -1)
set y2tics add ("-2 Ry" -2)
set y2tics add ("-3 Ry" -3)
set ytics add ("" 0)
set ytics add ("" 1)
set ytics add ("" 2)
set ytics add ("" 3)
set ytics add ("" -1)
set ytics add ("" -2)
set ytics add ("" -3)

set arrow 3 from   0,ph2e10 to 3000,ph2e10 as 2
set arrow 4 from   0,ph2e11 to 3000,ph2e11 as 2
set arrow 5 from   0,ph2e20 to 3000,ph2e20 as 3
set arrow 6 from   0,ph2e21 to 3000,ph2e21 as 3
set arrow 7 from 2000,ph2e10 to 2000,ph2e11 as 4
set arrow 8 from 1000,ph2e20 to 1000,ph2e21 as 5

set label 10 sprintf("W =%1.2f Ry",ph2e11-ph2e10) at  1500,(ph2e10-0.1) center font "Helvetica,16" tc lt 14
set label 11 sprintf("W'=%1.2f Ry",ph2e21-ph2e20) at  1500,(ph2e21+0.2) center font "Helvetica,16" tc lt 7

	if(i==2){mu = -0.5}

	set arrow 9 from 0,mu to 3000,mu nohead lt 1 lc 0 lw 5
	set label 5 "{/Symbol m}={/Symbol e}_a" at 1500,-0.2 center font "Helvetica,16"
	if(mu!=0){
          set label 5 sprintf("{/Symbol m}={/Symbol e}_a - %.1f Ry",-mu) at 1500,(mu-0.2) center font "Helvetica,16"
	  set arrow 10 from 0,0 to 3000,0 nohead lt -1 lc 0 lw 1 
          set label 6 "{/Symbol e}_a" at 1500,(0.2) center font "Helvetica,16"
        }

	#set title sprintf("To: %s",nameofphaseT)
# plot	sprintf("phase%s/phase%s_eps1_p=%s.dat",t,t,p) u 0:3 w l lw 5 lc 14 notitle,\
# 	sprintf("phase%s/phase%s_eps2_p=%s.dat",t,t,p) u 0:3 w l lw 5 lc 7 notitle
plot	sprintf("band%d.dat",i+1) u 0:4 w l lw 5 lc 14 notitle,\
	sprintf("band%d.dat",i+1) u 0:5 w l lw 5 lc 7 notitle

unset arrow 9
}


unset multiplot

reset


