# Copyright (c) 2015-2018. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import re
from pathlib import Path

from setuptools import find_packages, setup
from setuptools_scm import get_version

if sys.version_info < (3, 6):
    sys.exit("pygmst requires Python >= 3.6")

current_directory = os.path.dirname(__file__)

# with open("src/gtfparse/__init__.py", "r") as f:
#     version = re.search(
#         r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', f.read(), re.MULTILINE
#     ).group(1)

scm_version_options = {
    'write_to' : 'src/gtfparse/version.py',
    "root": ".",
    "relative_to": __file__,
    "local_scheme": "node-and-timestamp",
    "version_scheme": "post-release",
    }

setup(
    name="gtfparse",
    use_scm_version=scm_version_options,
    setup_requires=['setuptools_scm'],
    packages=find_packages(where="src"),
    package_dir={"gtfparse": "src/gtfparse"},
    package_data={"": ["gtfparse/tests/data/*.*"]},
    include_package_data=True,
    version=get_version(),
    description=(
        f"Parsing of General Transfer Format(GTF)/General Feature"
        f"Format 3 (GFF3) genetic annotation files"
    ),
    long_description=Path("README.rst").read_text("utf-8"),
    python_requires=">=3.6",
    url="https://github.com/milescsmith/gtfparse",
    author=["Alex Rubinsteyn", "Miles Smith"],
    license="http://www.apache.org/licenses/LICENSE-2.0.html",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
    ],
    install_requires=["numpy>=1.15", "pandas>=1.0.5", "tqdm>=4.31"],
    extras_require={"parallel": ["swifter~=0.3"]},
)
