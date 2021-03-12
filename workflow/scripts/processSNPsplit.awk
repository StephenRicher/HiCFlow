function header() {
	return  $0 ~ /^@/
}

BEGIN {
	if (! prefix) {
		prefix = ""
	}
    OFS="\t"
}

{
    if (header()) {
        next
    } else if ($7 != "=") { # Skip inter-chromosomal pairs
		next
	}
    sequence++
	# Get allelic assignment for each read
	split($NF, arr, ":")
	assign = arr[3]
	# Remap to default allele tags
	if (assign == "G1") {
		assign = "ref"
	} else if (assign == "G2") {
		assign = "alt"
	} else {
		assign = "both-ref"
	}
	# If first in pair, save assignment and move to 2nd
	if (sequence % 2 != 0) {
		assign1 = assign
		next
	} else {
		assign2 = assign
		# Define output file suffix
		if (assign1 == "ref") {
			if (assign2 == "ref") {
				suffix = "ref_ref"
			} else if (assign2 == "alt") {
				suffix = "ref_alt"
			} else {
				suffix = "ref_both-ref"
			}
		} else if (assign1 == "alt") {
			if (assign2 == "ref") {
				suffix = "ref_alt"
			} else if (assign2 == "alt") {
				suffix = "alt_alt"
			} else {
				suffix = "alt_both-ref"
			}
		} else {
			if (assign2 == "ref") {
				suffix = "ref_both-ref"
			} else if (assign2 == "alt") {
				suffix = "alt_both-ref"
			} else {
				suffix = "both-ref_both-ref"
			}
		}
		chr = $3
		out = prefix""chr"_"suffix
		# Fix column order, both-ref takes lowest priority
		if (assign1 != assign2 && assign1 != "both-ref") {
			if (assign2 == "both-ref" || assign2 == "alt") {
				print($3, $8, assign1, $3, $4, assign2, "") > out
			} else {
				print($3, $4, assign2, $3, $8, assign1, "") > out
			}
		} else {
			print($3, $4, assign2, $3, $8, assign1, "") > out
		}
	}
}
