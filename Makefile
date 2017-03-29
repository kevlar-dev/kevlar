depend:
	pip install pysam git+https://github.com/dib-lab/khmer.git

devenv: depend
	pip install pep8 pytest-cov

style:
	pep8 kevlar/*.py kevlar/*/*.py kevlar/*/*/*.py

test:
	py.test -v --cov=kevlar kevlar/tests/*.py -m 'not long'

testall:
	py.test -v --cov=kevlar kevlar/tests/*.py
