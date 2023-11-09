for FILE in output/*
do 
   if grep 'BrdU' ${FILE}; then
    ln -s ../${FILE} final/
   fi
done
