function header(line) {
	return  $line ~ /^@/
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
sam_record = $0
if (header(sam_record)) {
		print sam_record
} else {
	N++
	if (first_in_pair(N)) {
		split(sam_record, read1, "\t")
		next
	} else {
		split(sam_record, read2, "\t")
	}

	if (qname(read1) == qname(read2)) {
		print_record(read1)
		print_record(read2)
	} else {
		split(sam_record, read1, "\t")
		N++
	}
}
}
