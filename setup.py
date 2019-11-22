import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="RPx", # Replace with your own username
    version="0.0.1",
    author="Benjamin Sebastian",
    author_email="bensebastian@mac.com",
    description="A discrete time rhythmicity detector",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/benjaminwsebastian/RPx",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
