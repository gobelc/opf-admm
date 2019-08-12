#/bin/sh

cd $3/$1/
/usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nojvm -nosplash -r "update_x('/home/gonzalo/workspace/opf-admm/$1/',$2);quit"
