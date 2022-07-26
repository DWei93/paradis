odi=`pwd`
for DIR in set1 set2 set3; do
dir=0.8um/${DIR}/point-
cat data | while read line; do
    i=`echo ${line} | cut -d " " -f 1`
    if [ ! -d ../${dir}${i} ]; then
         mkdir ../${dir}${i} -p
    fi
    cp bin ../${dir}${i}/ -r
    sed -i "79d" input/Copper.ctrl
    polar=`echo ${line} | cut -d " " -f 7` 
    sed -i "79i    ${polar}" input/Copper.ctrl
    sed -i "80d" input/Copper.ctrl
    am=`echo ${line} | cut -d " " -f 8` 
    sed -i "80i    ${am}" input/Copper.ctrl
    echo "$i $polar $am"
    cp input ../${dir}${i}/ -r
    cp gen ../${dir}${i}/
    cp DWei-woody ../${dir}${i}/DWei-woody-${i}
    cd ../${dir}${i}; qsub DWei-woody-${i}; cd $odi
done
done
