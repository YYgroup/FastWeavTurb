echo "#!/bin/sh">submit.sh
echo "yhrun -N 1 -n 64 -p mars ./tube >out 2>&1">>submit.sh
chmod 700 submit.sh 
yhbatch -N 1 -n 64 -p mars -J weav.hzs  ./submit.sh
rm submit.sh
