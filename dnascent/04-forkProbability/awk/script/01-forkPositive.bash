for FILE in output/*
  do
    RETURN_VALUE=`awk -F"\t" '$2>0.7 || $3>0.7' ${FILE}`
      if [[ ! -z "${RETURN_VALUE}" ]]
          then
            ln -s ../${FILE} final/
      fi
  done
