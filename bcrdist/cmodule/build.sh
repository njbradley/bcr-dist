python setup.py build
for name in build/lib*; do
	mv $name/* build/
done
