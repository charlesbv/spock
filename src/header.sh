for f in *.c; do
  cat header $f > $f.new
  mv $f.new $f
done
