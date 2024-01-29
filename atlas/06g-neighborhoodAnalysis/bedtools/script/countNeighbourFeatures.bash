IFS='
'

OLD_ORIGIN_ID="-1"

CDS="0"
CDS_BASES="0"
PSEUDOGENE="0"
PSEUDOGENE_BASES="0"
SUM="0"
SUM_BASES="0"

for LINE in `cat output/mfa-intersect-genomeDisjointFeatures.tsv`
do
   ORIGIN_ID=`echo ${LINE} | cut -f4`
   if [ "${ORIGIN_ID}" -gt "${OLD_ORIGIN_ID}" ]
   then
      FIRST_WINDOW_START=`echo ${LINE} | cut -f2`
      FIRST_WINDOW_END=`echo ${LINE} | cut -f3`
      FEATURE=`echo ${LINE} | cut -f7`
      if [ "${FEATURE}" -eq "CDS" ]
      then
         let "CDS = CDS + 1"
         let "SUM = SUM + 1"
         FEATURE_BASES=`echo ${LINE} | cut -f14`
         let "CDS_BASES = CDS_BASES + ${FEATURE_BASES}"
         let "SUM_BASES = SUM_BASES + ${FEATURE_BASES}"
      elif [ "${FEATURE}" -eq "pseudogene" ]
      then
         let "PSEUDOGENE = PSEUDOGENE +1"
         let "SUM = SUM +1"
         FEATURE_BASES=`echo ${LINE} | cut -f14`
         let "PSEUDOGENE_BASES = PSEUDOGENE_BASES + ${FEATURE_BASES}"
         let "SUM_BASES = SUM_BASES + ${FEATURE_BASES}"
      fi
      elif
      then
      fi
      else
         echo "Feature not implemented yet: ${FEATURE}"
      fi
      FEATURE_ID=`echo ${LINE} | cut -f13 | cut -d';' -f1`
      OLD_ORIGIN_ID=${ORIGIN_ID}
   else
      SECOND_WINDOW_START=`echo ${LINE} | cut -f2`
      SECOND_WINDOW_END=`echo ${LINE} | cut -f3`

   fi
done
