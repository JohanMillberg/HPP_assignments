# Script for checking format of A2.tar.gz file for Assignment 2 subimssion

# Remove old tmpdir_for_checking if it exists
rm -rf tmpdir_for_checking

echo Creating temporary directory tmpdir_for_checking

# Create new tmpdir_for_checking directory
mkdir tmpdir_for_checking || exit 1

echo Copying A2.tar.gz file into tmpdir_for_checking directory

cp A2.tar.gz tmpdir_for_checking || exit 1

# cd into the tmpdir_for_checking directory
cd tmpdir_for_checking || exit 1

echo Unpacking the A2.tar.gz file
tar -xzf A2.tar.gz || exit 1

echo Doing cd into A2 directory
cd A2 || exit 1


if [ "$#" -ne 1 ];
then
    check_ass=0;
    echo "Checking all parts"
else
    check_ass=$1
    echo "Checking part " $check_ass
fi


if [ $check_ass -eq 0 ] || [ $check_ass -eq 1 ]
then

# Part 1

echo Doing cd into part1 directory
cd part1 || exit 1

echo Checking that Makefile is available in part1 directory
ls *akefile || exit 1

echo Checking that at least one .c file is available in part1 directory
ls *.c || exit 1

echo Checking that there are no .o files in the part1 directory
ls *.o 2> /dev/null && exit 1

echo Checking that there is no triang file in the part1 directory
ls triang 2> /dev/null && exit 1

echo Checking that doing make generates triang file in the part1 directory
make || exit 1
ls triang || exit 1

echo Check that running triang with input 11 really gives 11 lines of output
timeout 2 ./triang 11 | wc -l | grep 11 || exit 1

echo Check that running triang with input 12 gives OK output
./triang 12 | tail -1 | grep 462 || exit 1
./triang 12 | tail -1 | grep 330 || exit 1
./triang 12 | tail -1 | grep 165 || exit 1
./triang 12 | tail -2 | head -1 | grep 252 || exit 1

echo Checking that doing make clean really removes triang file in the part1 directory
make clean || exit 1
ls triang 2> /dev/null && exit 1

cd .. || exit 1

fi


if [ $check_ass -eq 0 ] || [ $check_ass -eq 2 ]
then

# Part 2

echo Doing cd into part2 directory
cd part2 || exit 1

echo Checking that Makefile is available in part2 directory
ls *akefile || exit 1

echo Checking that at least one .c file is available in part2 directory
ls *.c || exit 1

echo Checking that the little_bin_file file is available in part2 directory
ls little_bin_file || exit 1

echo Checking that there are no .o files in the part2 directory
ls *.o 2> /dev/null && exit 1

echo Checking that there is no readfile file in the part2 directory
ls readfile 2> /dev/null && exit 1

echo Checking that doing make generates readfile file in the part2 directory
make || exit 1
ls readfile || exit 1

echo Checking that running ./readfile in the part2 directory works and gives the expected output
./readfile || exit 1
./readfile | head -1 | grep [-]3477 > /dev/null || exit 1
./readfile | head -2 | tail -1 | grep 6892[.]367295 > /dev/null || exit 1
./readfile | head -3 | tail -1 | grep D > /dev/null || exit 1
./readfile | head -4 | tail -1 | grep 24[.]1 > /dev/null || exit 1

echo Checking that doing make clean really removes readfile file in the part2 directory
make clean || exit 1
ls readfile 2> /dev/null && exit 1

cd .. || exit 1


fi


if [ $check_ass -eq 0 ] || [ $check_ass -eq 3 ]
then

# Part 3

echo Doing cd into part3 directory
cd part3 || exit 1

echo Checking that Makefile is available in part3 directory
ls *akefile || exit 1

echo Checking that at least one .c file is available in part3 directory
ls *.c || exit 1

echo Checking that there are no .o files in the part3 directory
ls *.o 2> /dev/null && exit 1

echo Checking that there is no january file in the part3 directory
ls january 2> /dev/null && exit 1

echo Checking that doing make generates january file in the part3 directory
make || exit 1
ls january || exit 1


echo A 1 -15.2 -5.1 > tmp_test_input.txt
echo A 5 -1 1 >> tmp_test_input.txt
echo A 3 -11 -2 >> tmp_test_input.txt
echo D 2 >> tmp_test_input.txt
echo A 6 -4 -2 >> tmp_test_input.txt
echo D 3 >> tmp_test_input.txt
echo A 11 -8 -5 >> tmp_test_input.txt
echo A 6 -1 0 >> tmp_test_input.txt
echo A 12 -5 -2.3 >> tmp_test_input.txt
echo P >> tmp_test_input.txt
echo Q >> tmp_test_input.txt


echo Testing ./january program for the following usage example:
cat tmp_test_input.txt
cat tmp_test_input.txt | timeout 2 ./january > tmp_test_output.txt || exit 1
grep 1 tmp_test_output.txt | grep "[-]15[.]2" | grep "[-]5[.]1" || exit 1
grep 5 tmp_test_output.txt | grep "[-]1[.]0" | grep "1[.]0" || exit 1
grep 6 tmp_test_output.txt | grep "[-]1[.]0" | grep "0[.]0" || exit 1
grep 11 tmp_test_output.txt | grep "[-]8[.]0" | grep "[-]5[.]0" || exit 1
grep 12 tmp_test_output.txt | grep "[-]5[.]0" | grep "[-]2[.]3" || exit 1
grep ^2 tmp_test_output.txt 2> /dev/null && exit 1
grep ^3 tmp_test_output.txt 2> /dev/null && exit 1

rm -f tmp_test_output.txt


# Prepare test input file tmp_test_input_2.txt
echo A 15 -3.1 -1.5 > tmp_test_input_2.txt
echo A 17 0 3 >> tmp_test_input_2.txt
echo A 14 -1.3 4.8 >> tmp_test_input_2.txt
echo A 16 1.1 8.5 >> tmp_test_input_2.txt
echo A 29 0 6 >> tmp_test_input_2.txt
echo A 12 0 9 >> tmp_test_input_2.txt
echo D 17 >> tmp_test_input_2.txt
echo D 12 >> tmp_test_input_2.txt
echo D 29 >> tmp_test_input_2.txt
echo D 15 >> tmp_test_input_2.txt
echo P >> tmp_test_input_2.txt
echo Q >> tmp_test_input_2.txt

echo Testing ./january program for the following usage example:
cat tmp_test_input_2.txt
cat tmp_test_input_2.txt | timeout 2 ./january > tmp_test_output.txt || exit 1
grep 14 tmp_test_output.txt | grep "[-]1[.]3" | grep "4[.]8" || exit 1
grep 16 tmp_test_output.txt | grep "1[.]1" | grep "8[.]5" || exit 1
grep 17 tmp_test_output.txt 2> /dev/null && exit 1
grep 12 tmp_test_output.txt 2> /dev/null && exit 1
grep 29 tmp_test_output.txt 2> /dev/null && exit 1
grep 15 tmp_test_output.txt 2> /dev/null && exit 1

rm -f tmp_test_output.txt



# Prepare test input file tmp_test_input_3.txt
echo D 3 > tmp_test_input_3.txt
echo A 16 1.1 8.5 >> tmp_test_input_3.txt
echo D 3 >> tmp_test_input_3.txt
echo P >> tmp_test_input_3.txt
echo Q >> tmp_test_input_3.txt

echo Testing ./january program for the following usage example:
cat tmp_test_input_3.txt
cat tmp_test_input_3.txt | timeout 2 ./january > tmp_test_output.txt || exit 1
grep 16 tmp_test_output.txt | grep "1[.]1" | grep "8[.]5" || exit 1
grep 3 tmp_test_output.txt 2> /dev/null && exit 1

rm -f tmp_test_output.txt



# Prepare test input file tmp_test_input_4.txt
echo A 14 -1.3 4.8 > tmp_test_input_4.txt
echo A 16 1.1 8.5 >> tmp_test_input_4.txt
echo P >> tmp_test_input_4.txt
echo Q >> tmp_test_input_4.txt

echo Testing ./january program for the following usage example:
cat tmp_test_input_4.txt
cat tmp_test_input_4.txt | timeout 2 ./january > tmp_test_output.txt || exit 1
grep 14 tmp_test_output.txt | grep "[-]1[.]3" | grep "4[.]8" || exit 1
grep 16 tmp_test_output.txt | grep "1[.]1" | grep "8[.]5" || exit 1

rm -f tmp_test_output.txt



echo Checking that doing make clean really removes january file in the part3 directory
make clean || exit 1
ls january 2> /dev/null && exit 1

cd .. || exit 1

fi

# If we get to this point, then all the different tests above have passed.
echo
echo Congratulations, your A2.tar.gz file seems OK!
echo
