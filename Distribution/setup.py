from setuptools import setup, find_packages
setup(
    name = "pBasexQt",
    version = "1.0",
    packages = find_packages(),
    scripts = ['pBase.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = ['docutils>=0.3'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        'pBase': ['*.msg'],
    },

    # metadata for upload to PyPI
    author = "Nicolas Saquet",
    author_email = "nicolas.saquet@nottingham.ac.uk",
    description = "pBasex 2014 in python with a Qt font end",
    license = "LGPL",
    keywords = "pBaex, VMI, Image inversion",
    
    # could also include long_description, download_url, classifiers, etc.
)