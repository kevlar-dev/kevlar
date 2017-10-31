SHELL=bash

## #≠≠≠≠≠ build targets ≠≠≠≠≠#

## help:     print this help message and exit
help: Makefile
	@sed -n 's/^## //p' Makefile

## devenv:   install software development pre-requisites
devenv:
	pip install --upgrade pip setuptools pytest pytest-cov pep8 cython sphinx sphinx-argparse

## style:    check Python code style against PEP8
style:
	pep8 kevlar/*.py kevlar/*/*.py kevlar/*/*/*.py

## ext:      build C extensions
ext: kevlar/alignment.c src/align.c inc/align.h
	python setup.py build_ext --inplace

## test:     execute the automated test suite
test: ext
	py.test --cov=kevlar kevlar/tests/*.py -m 'not long'

## testall:  execute the automated test suite, including long-running tests
testall: ext
	py.test -v --cov=kevlar kevlar/tests/*.py

## doc:      build the documentation locally
doc:
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

## clean:    remove development artifacts
clean:
	rm -rf __pycache__/ biokevlar.egg-info/ build/ dist/ kevlar/__pycache__/ kevlar/*.cpython*.so
