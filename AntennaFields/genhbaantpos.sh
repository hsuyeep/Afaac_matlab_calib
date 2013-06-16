cat AntennaFieldCS002.conf.HBA | awk -F ' ' '{if (NR == 2) {xbase=$3; ybase=$4; zbase=$5} if (NR > 3 && NR < 51) printf "%.9f %.9f %.9f\n", $1+xbase, $2+ybase, $3+zbase}'
