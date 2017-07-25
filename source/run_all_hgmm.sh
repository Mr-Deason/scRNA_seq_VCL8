starttime=`date +'%Y-%m-%d %H:%M:%S'`

./get_barcode ~/downloads/hgmm_fastqs/config_cp.json

./error_correct_and_split ~/downloads/hgmm_fastqs/config_cp.json

kstart=`date +'%Y-%m-%d %H:%M:%S'`
./compute_TCCs ~/downloads/hgmm_fastqs/config_cp.json
kend=`date +'%Y-%m-%d %H:%M:%S'`
kstart_s=$(date --date="$kstart" +%s);
kend_s=$(date --date="$kend" +%s);
echo "kallisto running time: "$((kend_s-kstart_s))"s"

./prep_TCC_matrix ~/downloads/hgmm_fastqs/config_cp.json

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
echo "running time: "$((end_seconds-start_seconds))"s"
