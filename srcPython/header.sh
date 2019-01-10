for f in *.py; do
  cat header $f > $f.new
  mv $f.new $f
done
