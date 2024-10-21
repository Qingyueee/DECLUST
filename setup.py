from setuptools import setup, find_packages

setup(
    name='DECLUST',
    version='0.1.0',
    description='A tool for cell type deconvolution in spatial transcriptomics data',
    author='Qingyue',
    author_email='qingyue.wang@umail.ucc.ie',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scanpy',
        'matplotlib',
        'scikit-learn',
    ],
    entry_points={
        'console_scripts': [
            'declust=DECLUST.main:main',
        ],
    },
    package_data={
        '': ['data/*.csv', 'data/*.h5ad'],
        'scripts': ['R/*.R'],
    },
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
