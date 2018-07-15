#makefile
test:
	gcc pi3.c lhep.c bn_ext.c  test.c error_hdl.c -I /usr/local/include/flint  -lgmp -lflint -lrelic -lrelic_s -o test
.PHONY : clean
 clean:
	rm -f test
