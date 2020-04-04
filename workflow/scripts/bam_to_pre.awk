function header(line) {
	return  $line ~ /^@/
}

function reverse_strand(flag) {
	return and(flag, 0x10)
}

function first_in_pair(count) {
	return count % 2
}

{
  if (header($0)) {
		next
  } else {
    N++

    if (reverse_strand($2)) {
      s = 1
    } else {
        s = 0
    }

    if (first_in_pair(N)) {
      m = $2
      printf "%s\t%s\t%s\t%s\t0\t", $1, s, $3, $4
    } else {
      printf "%s\t%s\t%s\t1\t%s\t%s\n", s, $3, $4, m, $5
    }
  }
}
