#/bin/sh

cd $2/$1/
/usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nojvm -nosplash -r "update_y('/home/gonzalo/workspace/opf-admm/$1/');quit"
