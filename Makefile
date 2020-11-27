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

docs:
	cd docs; make html


############ For building the python binding

getpytools:
	#to build python packages
	python3 -m pip install --upgrade setuptools wheel

pypkg:
	#build the python package
	cd py; python3 setup.py sdist bdist_wheel

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


test_isocount_build:
	java -jar build/isocount.jar build /home/mahogny/Dropbox/applyPI/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3.gz test/featurefile_test.ff

test_isocount_count:
	java -jar build/isocount.jar count test/featurefile_test.ff test/SRR11827037_out test/out.SRR11827037_out




clean:
	rm -Rf py/scherlock_pkg_USER.egg-info py/build py/dist
