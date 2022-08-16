py_ext: venv libalign_reads.a align_reads_gen.cpp 
	. align_reads_venv/bin/activate; python3 setup.py build_ext -f; python3 setup.py install

venv: align_reads_venv/bin/activate

align_reads_venv/bin/activate:
	python3 -m venv align_reads_venv
	. align_reads_venv/bin/activate; pip install pip --upgrade; pip install -r requirements.txt

libalign_reads.a: build/libalign_reads.a venv bioparser/CMakeLists.txt

bioparser/CMakeLists.txt: 
	git submodule init && git submodule update

build/libalign_reads.a: build/Makefile src/aligner.cpp tests/Inputs_tests.cpp
	cmake --build build/

build/Makefile: bioparser/CMakeLists.txt CMakeLists.txt
	. align_reads_venv/bin/activate; cmake -S . -B build/

rebuild:  
	cmake --build build/
	. align_reads_venv/bin/activate; python3 setup.py build_ext -f; python3 setup.py install
