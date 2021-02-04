

# Compile
javac -cp /corgi/anbjork/scherlock/build/cellpile.jar CellPileTestWrite.java
javac -cp /corgi/anbjork/scherlock/build/cellpile.jar CellPileTestReadJohan.java
javac -cp /corgi/anbjork/scherlock/build/cellpile.jar CellPileTestReadAnton.java

# Run write on test data
# See paths to data in .java file
java -cp /corgi/anbjork/scherlock/build/cellpile.jar:. CellPileTestWrite

# The two read ones are both outdated, so I didn't update them, since
# they seem to not be in use.
# The compilation and running should be similar to write test

