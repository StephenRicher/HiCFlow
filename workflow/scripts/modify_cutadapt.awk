{
  if (sample != "") {
    if ($0 ~ /^Command line parameters/) {
		    $0 = $0" # "sample
      }
  }
  print $0
}
