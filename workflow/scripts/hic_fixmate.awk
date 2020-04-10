function header(line) {
	return  $line ~ /^@/
}

function unmapped(flag) {
	return and(flag, 4)
}

function low_quality(quality) {
	return quality < min_quality
}

function first_in_pair(count) {
	return count % 2
}

function print_record(read) {
	for (field in read) {
		printf "%s\t", read[field] > out
	}
	printf "\n" > out
}

BEGIN {
    FS="\t"
    OFS="\t"
		if (min_quality == "") {
			min_quality = 15
		}
		if (default_out == "") {
			default_out = "/dev/stdout"
		}
		if (stats_out == "") {
			stats_out = "/dev/stderr"
		}
		if (unmapped_out == "") {
			unmapped_out = "/dev/null"
		}
		if (low_quality_out == "") {
			low_quality_out = "/dev/null"
		}
}


{
sam_record = $0
if (header(sam_record)) {
	if (sam_record != previous_header) {
		print sam_record > default_out
		print sam_record > unmapped_out
		print sam_record > low_quality_out
	}
	previous_header = sam_record
} else {
	N++
	if (first_in_pair(N)) {
		$2 = $2 + 0x1 + 0x40
		r1_qual = $5
		r1_flag = $2
		read1 = $0
		next
	} else {
		$2 = $2 + 0x1 + 0x80
		r2_qual = $5
		r2_flag = $2
		read2 = $0
	}

	if (unmapped(r1_flag) || unmapped(r2_flag)) {
		unmap++
		out = unmapped_out
	} else if (low_quality(r1_qual) || low_quality(r2_qual)) {
		low_qual++
		out = low_quality_out
	} else {
		out = default_out
	}

	print read1 > out
	print read2 > out
}
}

END {
	printf "Total Reads: %s\n", N/2 > stats_out
	printf "Unmapped: %s\n", unmap > stats_out
	printf "Low quality: %s\n", low_qual > stats_out
}
