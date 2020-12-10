#!/usr/bin/awk -f

function header() {
	return $0 ~ /^@/
}

function first_in_pair(count) {
	return count % 2
}

function qname(read) {
	return read[1]
}

function print_record(read) {
	for (field in read) {
		printf "%s\t", read[field]
	}
	printf "\n"
}

BEGIN {
    FS="\t"
    OFS="\t"
}

{
if (header()) {
	print $0
} else {
	N++
	if (first_in_pair(N)) {
		split($0, read1, "\t")
		next
	} else {
		split($0, read2, "\t")
		if (qname(read1) == qname(read2)) {
			print_record(read1)
			print_record(read2)
		} else {
			split($0, read1, "\t")
			N++
		}
	}
}
}
