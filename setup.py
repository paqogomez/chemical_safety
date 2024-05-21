from setuptools import setup, find_packages

setup(
    name="chemical_safety",
    version="0.0.7",
    description="A package for retreiving chemical safety information",
    author="Demetrios Pagonis",
    author_email="demetriospagonis@weber.edu",

    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(include=['chemical_safety', 'chemical_safety.*']),
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy',
        'jinja2',
        'requests',
        'Levenshtein',
        'natsort',
        'flask',
        'scipy',
        'rdkit'
    ],
    entry_points={
        'console_scripts': [
            'chemical-dashboard=chemical_safety.dashboard.app:dashboard',
        ],
    },
    url='https://github.com/dpagonis/chemical_safety',
)
