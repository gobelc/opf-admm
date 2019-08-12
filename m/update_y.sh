#/bin/sh

cd /home/olaznog/workspace/opf-admm/$1/
/usr/local/MATLAB/R2017a/bin/matlab -nodisplay -nojvm -nosplash -r "update_y('/home/olaznog/workspace/opf-admm/$1/');quit"
