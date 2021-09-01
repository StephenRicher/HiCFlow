function header() {
    return  $0 ~ /^@/
}

function unpaired() {
    return (and($2, 4) != 0) || (and($2, 8) != 0)
}

BEGIN {
    OFS="\t"
    prevHeader = ""
}

{
    if (header()) {
        print $0
    } else if (unpaired()) {
        pass
    } else if (prevHeader == $1) {
        print prevLine
        print $0
    } else {
        prevHeader = $1
        prevLine = $0
    }
}
