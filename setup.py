import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lammps_helper-kevinwhitham",
    version="0.1.1",
    author="Kevin Whitham",
    author_email="kevin.whitham@gmail.com",
    description="Helps with LAMMPS input/output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kevinwhitham/lammps_helper",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
