.PHONY: test clean

test:
	python -m pytest -q

clean:
	python -c "import pathlib,shutil; \
[shutil.rmtree(p) for p in pathlib.Path('.').rglob('__pycache__')]; \
[shutil.rmtree(p) for p in pathlib.Path('.').rglob('.pytest_cache')]"
