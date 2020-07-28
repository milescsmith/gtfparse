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
import re
import sys

from setuptools import setup, find_packages

if sys.version_info < (3, 6):
    sys.exit("pygmst requires Python >= 3.6")

readme_filename = "README.md"
current_directory = os.path.dirname(__file__)
readme_path = os.path.join(current_directory, readme_filename)

readme_markdown = ""
try:
    with open(readme_path, "r") as f:
        readme_markdown = f.read()
except Exception as e:
    print(e)
    print(f"Failed to open {readme_path}")

try:
    import pypandoc

    readme_restructured = pypandoc.convert(readme_markdown, to="rst", format="md")
except Exception as e:
    readme_restructured = readme_markdown
    print(e)
    print(f"Failed to convert {readme_filename} from Markdown to reStructuredText")

with open("src/gtfparse/__init__.py", "r") as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', f.read(), re.MULTILINE
    ).group(1)


setup(
    name="gtfparse",
    packages=find_packages(where="src"),
    package_dir={"gtfparse": "src/gtfparse"},
    package_data={"": ["gtfparse/tests/data/*.*"]},
    include_package_data=True,
    version=version,
    description="GTF Parsing",
    long_description=readme_restructured,
    python_requires=">=3.6",
    url="https://github.com/openvax/gtfparse",
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
    extra_requires={"parallel": ["swifter~=0.3"]},
)
