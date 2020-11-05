compile:
	mkdir -p out
	mkdir -p build
	mkdir -p bin
	javac src/*/*/*.java -cp .:lib/colt.jar:lib/concurrent.jar:lib/htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar:lib/py4j0.10.9.1.jar -d out
	cd out; jar cmf manifests/count.mf ../build/isocount.jar *
	cd out; jar cmf manifests/cellpile.mf ../build/cellpile.jar *
	cd out; jar cmf manifests/python.mf ../build/pycellpile.jar *
	rm -Rf build/lib
	cp -r lib build/
	cp py/* build/  #to help testing

loc:
	wc  -l src/*/*/*java

gitaddall:
	git add src/*.java



############ For building the python binding

getpytools:
	python3 -m pip install --upgrade setuptools wheel

pypkg:
	python3 setup.py sdist bdist_wheel

testpypkg:
	python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps cellpile

uploadpy:
	#https://packaging.python.org/tutorials/packaging-projects/
	#python3 -m pip install --upgrade twine
	python3 -m twine upload --repository testpypi dist/*

topypi:
	twine upload dist/*

#-Xmx6G

#~/miniconda3/share/py4j/py4j0.10.9.1.jar
