BEGIN {
} {
  if (k != $1) {
    if (k) {
      print k" "v;
    }
    k = $1;
    v = "";
  }
  if (v) {
    v = v" "$2;
  } else {
    v = $2;
  }
} END {
  if (k) {
    print k" "v;
  }
}
