starttime=`date +'%Y-%m-%d %H:%M:%S'`

./get_barcode ../example_dataset/config.json
./error_correct_and_split ../example_dataset/config.json
./compute_TCCs ../example_dataset/config.json
./prep_TCC_matrix ../example_dataset/config.json

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
echo "running time: "$((end_seconds-start_seconds))"s"