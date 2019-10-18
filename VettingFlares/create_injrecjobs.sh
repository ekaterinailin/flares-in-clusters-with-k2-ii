 

OLDIFS=$IFS
IFS=,
END=100
while read ids campaign
    
    do  
        count=0
        while [ "$count" -le 10 ]
            do 
                
                echo python /home/eilin/InjectionRecoveryK2/src/injrec.py $ids $campaign $count >> injrec.sh #/home/eilin/src/
                count=$(($count+1))
            done 
    done < targets.csv 
    

IFS=$OLDIFS

#grep -v "EPIC" joblist.sh > temp && mv temp joblist.sh
