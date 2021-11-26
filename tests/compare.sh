# @param $1: known tutorial output
# @param $2: workflow output
diff -bws <(sort$1) <(sort$2)
# -bw ignores spacing differences
# -s reports if the two are same