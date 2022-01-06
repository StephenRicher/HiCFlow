function header() {
	return  $0 ~ /^@/
}

BEGIN {
    if (! flag) {
        print "Flag modifier not set (e.g -v flag=0x40)"
        exit 1
    }
    OFS="\t"
}

{
    if (! header()) {
        $2 = $2 + strtonum(flag)
		$1 = NR
    }
    print $0
}
