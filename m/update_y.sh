#/bin/sh

cd $2/$1/
/usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nojvm -nosplash -r "update_y('$2/$1/');quit"
