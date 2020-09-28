compile:
	mkdir -p out
	mkdir -p build
	javac src/*/*/*.java -cp .:lib/colt.jar:lib/concurrent.jar:lib/htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar -d out
	cd out; jar cmf manifests/count.mf ../build/isocount.jar *
	cd out; jar cmf manifests/cellpile.mf ../build/cellpile.jar *

gitaddall:
	git add src/*.java


#-Xmx6G
