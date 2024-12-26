from setuptools import setup, find_packages

setup(
    name="pyburn",  # The name of our package
    version="0.1",  # Version of our package
    description="A tool for calculating Burn Severity using Sentinel-2 imagery",  # Short description of our package
    author="Team pyburn",  # our Team name 
    author_email="praveenkumar.saminathan@mail.polimi.it",  # our email address
    url="https://github.com/Orange3456/pyburn",  # The URL to our project (hosted on GitHub)
    packages=find_packages(),  # Automatically find all packages in the current directory
    install_requires=[  # Any external packages our package depends on
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "rasterio",
        "geopandas",
        "fiona",
        "shapely",
        "gdal",
    ],
    classifiers=[  # Optional: classifiers to categorize our package on PyPI
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',  # The Python versions supported by our package
)
