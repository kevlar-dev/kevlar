SHELL=bash


## #≠≠≠≠≠ build targets ≠≠≠≠≠#

## help:     print this help message and exit
help: Makefile
	@sed -n 's/^## //p' Makefile

## devenv:   install software development pre-requisites
devenv:
	pip install --upgrade pip setuptools 'pytest>=4.0,<5.0' pytest-cov pytest-xdist pycodestyle cython sphinx sphinx-argparse

## style:    check Python code style against PEP8
style:
	pycodestyle --exclude=kevlar/sandbox/*.py kevlar/*.py kevlar/*/*.py kevlar/*/*/*.py

## ext:      build C extensions
ext: kevlar/alignment.pyx kevlar/sequence.pyx kevlar/assembly.pyx src/align.c inc/align.h
	python setup.py build_ext --inplace

## test:     execute the automated test suite
test: ext
	pytest --cov=kevlar kevlar/tests/test_*.py -m 'not long and not toolong'

## test4:    execute the automated test suite with 4 parallel threads
test4: ext
	pytest -n=4 --cov=kevlar kevlar/tests/test_*.py -m 'not long and not toolong'

## testmore: execute the automated test suite, including longer-running tests
testmore: ext
	pytest -v --cov=kevlar kevlar/tests/test_*.py -m 'not toolong'

## testall:  execute the automated test suite, including all tests
testall: ext
	pytest -v --cov=kevlar kevlar/tests/test_*.py

## doc:      build the documentation locally
doc: ext
	PYTHONPATH=$$(pwd) make -C docs/ html

## loc:      compute the lines of code
loc:
	@- echo -e "\n\n===== Core kevlar ====="
	cloc --exclude-list-file=.cloc.exclude kevlar/*.py
	@- echo -e "\n\n===== kevlar C extensions ====="
	cloc kevlar/*.pyx src/*.c inc/*.h
	@- echo -e "\n\n===== kevlar CLI ====="
	cloc kevlar/cli/*.py
	@- echo -e "\n\n===== kevlar tests ====="
	cloc kevlar/tests/test_*.py
	@- echo -e "\n\n===== Sandbox scripts ====="
	cloc kevlar/sandbox/*.py

## clean:    remove development artifacts
clean:
	rm -rf __pycache__/ biokevlar.egg-info/ build/ dist/ kevlar/__pycache__/ kevlar/*.cpython*.so
