BEGIN {
    OFS="\t"
}

{
    print ".", prefix""$3, $4, "+", prefix""$3, $8, "+"
}
