
OLDIFS=$IFS
IFS=,
while read ids channel campaign; do echo python script_parallel.py $ids $channel $campaign ngc2168 >> joblist.sh ; done < ngc2168_altai.csv 
while read ids channel campaign; do echo python script_parallel.py $ids $channel $campaign ngc2682 >> joblist.sh ; done < ngc2682_altai.csv 
while read ids channel campaign; do echo python script_parallel.py $ids $channel $campaign pleiades >> joblist.sh ; done < pleiades_altai.csv 
while read ids channel campaign; do echo python script_parallel.py $ids $channel $campaign hyades >> joblist.sh ; done < hyades_altai.csv 
while read ids channel campaign; do echo python script_parallel.py $ids $channel $campaign praesepe >> joblist.sh ; done < praesepe_altai.csv 
while read ids channel campaign; do echo python script_parallel.py $ids $channel $campaign ngc6774 >> joblist.sh ; done < ngc6774_altai.csv 
IFS=$OLDIFS

grep -v "EPIC" joblist.sh > temp && mv temp joblist.sh
