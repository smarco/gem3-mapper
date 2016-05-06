#!/bin/bash

counter=../bin/gt.gtfcount
function test_gtfcount(){
    out="$(mktemp -t gtfcount_result)"
    # test default paired end run
    target=$1; shift;
    echo -n "Runing with:: target:$target params:'$@' "
    $counter -i testdata/counts.map -a testdata/counts.gtf -g $out $@ >/dev/null
    diff $out testdata/$target || { echo "Default paramter test failed, check $out" >&2; exit 1; }
    rm $out
    echo ": Done"
}

test_gtfcount result_counts_p.txt -p
test_gtfcount result_counts_p_s.txt -p -s
test_gtfcount result_counts_p_m.txt -p -m
test_gtfcount result_counts_p_m_s.txt -p -m -s
test_gtfcount result_counts_p_e_1.0.txt -p -e 1.0
test_gtfcount result_counts_p_e_0.5.txt -p -e 0.5
test_gtfcount result_counts_p_e_0.5_s.txt -p -e 0.5 -s
test_gtfcount result_counts_p_e_1.0_m.txt -p -e 1.0 -m
test_gtfcount result_counts_p_m_w.txt -p -m -w
test_gtfcount result_counts_p_m_w_s.txt -p -m -w -s

