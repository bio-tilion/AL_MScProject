BEGIN {
    FS="\t"; OFS="\t"
}
NR==FNR {
    newcol[NR]=$4
    next
}
{
    $3=$3"\t"newcol[FNR]
    print $0
}