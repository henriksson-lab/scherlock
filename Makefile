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

gitaddall:
	git add src/*.java


#-Xmx6G

#~/miniconda3/share/py4j/py4j0.10.9.1.jar
