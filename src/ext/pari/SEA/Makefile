distrib:
	cd ..; tar czf ellsea.tar.gz --exclude CVS SEA
test::
	gp -q -f --test < test/seain > gp.out
	diff -c test/sea32 gp.out
