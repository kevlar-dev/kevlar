depend:
	pip install pysam git+https://github.com/dib-lab/khmer.git

devenv: depend
	pip install pep8 pytest

style:
	pep8 kevlar/*.py kevlar/*/*.py tests/*.py tests/data/*.py

test:
	py.test -v --cov=kevlar tests/*.py -m 'not long'

testall:
	py.test -v --cov=kevlar tests/*.py
