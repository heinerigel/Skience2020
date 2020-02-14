#!/bin/csh
set set R1 = "-R0/40000/-1.e-10/1.e-10"
set proj=-JX15c/10c
#
set list=`ls snapshot???`
set i=0
foreach f($list)
@ i++
echo $i", "$f
psxy $R1 $proj -W2pt -BS10000W1 -Y8c    $f > tmp.ps
if (-e tmp.eps) \rm tmp.eps
ps2eps  tmp.ps >& /tmp/null 
mv tmp.eps tmp"$i".eps
#convert -density 150x150 -rotate 90 tmp.eps snapshot$i.png
end
convert -dispose Background -delay 40 -rotate 90 tmp?.eps tmp??.eps tmp???.eps   -delay 1000 movie.gif
eog movie.gif
\rm tmp*.eps
