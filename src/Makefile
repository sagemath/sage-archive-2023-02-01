clean:
	@echo "Deleting Sage library build artifacts..."
	rm -f c_lib/.sconsign.dblite
	find c_lib -name '*.os' | xargs rm -f
	find c_lib -name '*.so' | xargs rm -f
	find c_lib -name '*.dylib' | xargs rm -f
	rm -rf build
	find . -name '*.pyc' | xargs rm -f
