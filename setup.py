from setuptools import setup, find_packages

with open("docs/README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [
        line.strip() for line in fh 
        if line.strip() and not line.startswith("#") and "setuptools" not in line.lower()
    ]

setup(
    name="spo2-algorithm",
    version="0.1.0",
    author="Mario Gómez-Barea",
    author_email="m.gomez@onalabs.com",
    description="SpO2 calculation algorithms from wrist-worn PPG sensors with IMU-based motion artifact removal",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/OnalabsInnoHub/spo2-algorithm",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.9,<3.15",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=8.4.0",
            "pytest-cov>=7.0.0",
            "black>=25.0.0",
            "flake8>=7.0.0",
            "mypy>=1.18.0",
        ],
        "docs": [
            "sphinx>=8.0.0",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)