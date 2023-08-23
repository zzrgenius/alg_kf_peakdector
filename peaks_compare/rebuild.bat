make clean
rm ./pd.exe
rm ./*.txt
cp ./data/*.txt ./
make
.\pd.exe
cp *.txt E:\work\work_project\develop\alg_four_in_one\data_check\filter_compare