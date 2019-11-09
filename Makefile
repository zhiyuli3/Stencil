stencil: stencil.c 
	#gcc  -O0 -std=c99 -Wall $^ -o $@
	#gcc  -O1 -std=c99 -Wall $^ -o $@
	#gcc  -O2 -std=c99 -Wall $^ -o $@
	#gcc  -O1 -ftree-vectorize -std=c99 -Wall $^ -o $@
	#gcc  -O2 -ftree-vectorize -std=c99 -Wall $^ -o $@
	#gcc  -O3 -std=c99 -Wall $^ -o $@
	# icc -O0 -std=c99  stencil.c -o stencil
	# icc -O1 -std=c99  stencil.c -o stencil
	# icc -O2 -std=c99  stencil.c -o stencil
	# icc -O3 -std=c99  stencil.c -o stencil
	# icc -O2 -no-vec -std=c99  stencil.c -o stencil
	# icc -O3 -no-vec -std=c99  stencil.c -o stencil
	# branch mpi
	icc -Ofast -std=c99  stencil.c -o stencil