depend:
	pip install pysam git+https://github.com/dib-lab/khmer.git

devenv: depend
	pip install pep8 pytest

style:
	pep8 kevlar/*.py bin/kevlar tests/*.py

test:
	py.test -v tests/*.py
