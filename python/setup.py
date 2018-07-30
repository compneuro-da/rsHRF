import re
from setuptools import setup, find_packages

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('rsHRF/CLI.py').read(),
    re.M
).group(1)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="rsHRF",
    packages=find_packages(),
    entry_points={
        "console_scripts": ['rsHRF = rsHRF.CLI:main']
    },
    version=version,
    description="BIDs App to retrieve the haemodynamic response function from resting state fMRI data",
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Madhur Tandon",
    author_email="madhurtandon23@gmail.com",
    url="https://github.com/guorongwu/rsHRF",
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    zip_safe=False,
    python_requires=">=3.5",
    install_requires=["numpy>=1.14,<1.15", "nibabel", "matplotlib", "scipy", "pybids", "pandas", "patsy", "duecredit",
                      "joblib"],
)
