devenv:
	pip install --upgrade pip setuptools pytest pytest-cov pep8 cython sphinx

style:
	pep8 kevlar/*.py kevlar/*/*.py kevlar/*/*/*.py

test:
	py.test -v --cov=kevlar kevlar/tests/*.py -m 'not long'

testall:
	py.test -v --cov=kevlar kevlar/tests/*.py

doc:
	cd docs && make html

loc:
	@- echo "\n\n===== Core kevlar ====="
	cloc --exclude-list-file=<(echo kevlar/_version.py) kevlar/*.py
	@- echo "\n\n===== kevlar CLI ====="
	cloc kevlar/cli/*.py
	@- echo "\n\n===== kevlar tests ====="
	cloc kevlar/tests/test_*.py
