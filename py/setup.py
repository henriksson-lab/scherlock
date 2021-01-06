import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scherlock-pkg-USER", # Replace with your own username
    version="0.0.1",
    author="Johan Henriksson, Anton Björk, Henlab",
    author_email="author@example.com",
    description="A package for detailed single-cell analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/henriksson-lab/isocounter",
    packages=setuptools.find_packages(),
    install_requires=[
        'pandasql','pyjnius','IPython',
        'matplotlib','plotly',
        'numpy','pandas','validators','beautifulsoup4',
        'anndata','scanpy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False,
    python_requires='>=3.6',
    include_package_data=True,
    package_data={'': ['data/*.jar']},
)

