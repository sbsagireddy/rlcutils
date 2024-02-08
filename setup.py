from pathlib import Path
from setuptools import find_packages, setup

# Load version number
__version__ = ''
version_file = Path(__file__).parent.absolute() / 'rlcutils' / '_version.py'

with open(version_file) as fd:
    exec(fd.read())

# Load README
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='rlcutils',
    version=__version__,
    author='Srikant Sagireddy',
    author_email='sbsagireddy@gmail.com',
    description='RLC utils',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/sbsagireddy/rlcutils',
    project_urls={
        'Source': 'https://github.com/sbsagireddy/rlcutils',
    },
    license='MIT',
    packages=find_packages(),
    package_data={'rlcutils': ['py.typed']},
    install_requires=[
        'mdanalysis',
        'numpy',
        'pandas',
        'scipy',
    ],
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    keywords=[
        'polymer physics'
    ]
)
