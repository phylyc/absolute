mkdir ABSOLUTEv1.3

cp ABSOLUTE_stub.R ABSOLUTEv1.3
cp batch_exec_ABSv1.3.R ABSOLUTEv1.3
cp README ABSOLUTEv1.3
cp -rL ~/ABS_lib ABSOLUTEv1.3/library

tar -zcvf ABSOLUTEv1.3.tar.gz ABSOLUTEv1.3

rm -rf ABSOLUTEv1.3
