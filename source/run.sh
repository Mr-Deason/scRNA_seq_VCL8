t0=`date +'%Y-%m-%d %H:%M:%S'`
s0=$(date --date="$t0" +%s);

./get_barcode $1

t1=`date +'%Y-%m-%d %H:%M:%S'`
s1=$(date --date="$t1" +%s);

./new_error_correct_and_split $1

t2=`date +'%Y-%m-%d %H:%M:%S'`
s2=$(date --date="$t2" +%s);

./compute_TCCs $1

t3=`date +'%Y-%m-%d %H:%M:%S'`
s3=$(date --date="$t3" +%s);

echo "kallisto running time: "$((s3-s2))"s"

./prep_TCC_matrix_simd $1

t4=`date +'%Y-%m-%d %H:%M:%S'`
s4=$(date --date="$t4" +%s);
echo "running time: "$((s4-s0))"s"

echo ""$((s1-s0))", "$((s2-s1))", "$((s3-s2))", "$((s4-s3))", "$((s4-s0))""
