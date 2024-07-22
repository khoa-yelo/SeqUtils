from setuptools import setup, find_packages

setup(
    name='SeqUtils',
    version='0.1',
    description='Utilities for working with biological sequences.',
    author='Khoa Hoang',
    author_email='your.email@example.com',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'biopython'
    ],
    entry_points={
        'console_scripts': [
            'SeqUtils=SeqUtils.seq_utils:read_fasta',
        ],
    },
    python_requires='>=3.8',
)
