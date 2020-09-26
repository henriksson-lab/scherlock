compile:
	mkdir -p out
	javac src/* -cp htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar:. -d out

gitaddall:
	git add src/*.java


