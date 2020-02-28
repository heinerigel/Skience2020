#!/bin/csh
set list='E_av_test  E_harmonic_av_test  periodic_test  periodic_test2  velocity_av_test'
foreach d($list)
cd $d
\rm snapshot* velocity mu out  density  source.gnu velocity trace* >& /dev/null
cd ..
end
pushd ../sem1d/src
make clean
popd
\rm ../sem1d/bin/sem1d  >& /dev/null

