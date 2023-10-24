from setuptools import setup, find_packages

setup(
    name="markerrepo",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'pyyaml',
        'tabulate',
        'ipykernel',
        'GitPython',
        'scikit-learn',
        'pybiomart',
        'scanpy',
    ],
)
