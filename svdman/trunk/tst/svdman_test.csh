#!/bin/csh -f
#
#
# If this doesn't work, try adding a "-minus" directive to the command line
#
../src/svdman -i data.tab -o test -evaluate -writesvd
diff data_B.dat test_B.dat > diff_B.txt
diff data_VT.dat test_VT.dat > diff_VT.txt
diff data_U.dat test_U.dat > diff_U.txt
diff data_groups.dat test_groups.dat > diff_groups.txt
diff data_evaluate.dat test_evaluate.dat > diff_evaluate.txt

