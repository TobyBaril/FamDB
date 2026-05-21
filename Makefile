.PHONY: check coverage

check:
	FAMDB_TEST_BLESS=1 python3 -m unittest tests.test_cli

coverage:
	FAMDB_TEST_COVERAGE=1 coverage run -m unittest
	coverage combine
	coverage html --omit='*/site-packages/*'
