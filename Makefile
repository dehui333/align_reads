py_ext: venv libalign_reads.a align_reads_gen.cpp
	. align_reads_venv/bin/activate; python3 setup.py build_ext -f; python3 setup.py install

venv: align_reads_venv/bin/activate

align_reads_venv/bin/activate:
	python3 -m venv align_reads_venv
	. align_reads_venv/bin/activate; pip install pip --upgrade; pip install -r requirements.txt

libalign_reads.a: build/libalign_reads.a venv

build/libalign_reads.a: build/Makefile src/aligner.cpp
	cmake --build build/

build/Makefile:
	align_reads_venv/bin/activate; cmake -S . -B build/


	
