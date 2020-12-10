#!/usr/bin/awk -f

function inHeader() {
	return $0 ~ /^@/
}

function first_in_pair(count) {
	return count % 2
}

BEGIN {
    FS="\t"
    OFS="\t"
}

{
if (endHeader || inHeader()) {
	print $0
} else {
	endHeader=0
	N++
	if (first_in_pair(N)) {
		read1qname = $1
		read1 = $0
		next
	} else {
		read2qname = $1
		read2 = $0
		if (read1qname == read2qname) {
			print read1
			print read2
		} else {
			read1 = read2
			read1qname = read2qname
			N++
		}
	}
}
}
