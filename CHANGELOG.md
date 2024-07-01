# [2.2.0] [2024-07-01]

## Changed:
- Removed modin optional dep
- Updated dependencies
    - Simplified dev dependencies, replacing black and isort with ruff
    - Removed pyre
- Switched from poetry to pdm for package management

# [2.1.0] [2023-11-24]:

## Changed:
- shifted a few tuples to sets to make pyre happy
- added typechecking by pyre


# [2.0.0] [2023-11-21]:

## Changed:
- Minimum Python to 3.10
- Replaced unittest with pytest
- Replaced tox with nox
- Replaced flake8, isort, pylint, pyre-check with ruff
- Replaced relative dependencies
- Update dependencies
- Move all config files into pyproject.toml
- Replace coloredlogs with loguru
- Updated typing

## Added:
- Added .gitignore

[2.1.0]: https://github.com/milescsmith/gtfparse/compare/2.0.0...2.1.0
[2.0.0]: https://github.com/milescsmith/gtfparse/releases/tag/2.0.0