#/bin/sh
datestr=`date "+%m%d%Y-%H-%M"`
for i in {0,1,2,3,4,5,6,8,9,10}
do 
    echo $i
    #./a.out $i | tee ../Benchmarks/benchmark-$i-$datestr.txt
done
./a.out 7 seidel-skip-list.txt | tee ../Benchmarks/benchmark-7-$datestr.txt

