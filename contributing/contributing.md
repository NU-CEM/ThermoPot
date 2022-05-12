## Contributing to ThermoPot

Thanks for your interest in ThermoPot - we welcome your help in improving and extending this package. 

This project follows the [all-contributors](https://allcontributors.org/) specification. Contributions of any kind are welcome! A full list of possible contribution types can be found [here](https://allcontributors.org/docs/en/emoji-key). **You do not need to be an experienced programmer to help improve ChooChoo.**

All contributions will be recognised on the README.md.
You are encouraged to log your own contribution using the [all-contributors bot](https://allcontributors.org/docs/en/bot/usage). The project lead(s) will also maintain this list.

ChooChoo contributors are asked to follow the [Contributor Covenant Code of Conduct](https://github.com/NU-CEM/ThermoPot/blob/main/CODE_OF_CONDUCT.md).

## Contributions workflow

Code contributions are primarily managed through Github pull requests. For external contributions we prefer the “fork and pull” workflow:

1. open an Issue to discuss the proposed contribution. 
2. make your own project fork and implement the changes there
3. open a pull requestion (PR) to merge the changes into the main project. Further discussion might also take place at this stage.

Note: Please, where applicable and possible, write tests and documentation for your proposed code contributions (further details below). Support can be given to those who are writing tests and documentation for the first time :) 
        
## Communication

As a rule of thumb, the more public the information exchange, the better! With this in mind, we encourage you to [raise an issue](https://github.com/NU-CEM/ThermoPot/issues/new/choose) on the ThermoPot repository.

However we also understand that this public way of working is not suitable for everyone. As an alternative, you can get [get in contact](https://lucydot.github.io/about/) with the project lead Lucy. 

## Testing

Github Actions and [Pytest](https://docs.pytest.org/) is used to automatically [test all commits](https://github.com/NU-CEM/ThermoPot/actions/workflows/run-tests.yml) to the main branch of the ThermoPot repository. If possible, please write tests for your code when issuing a Pull Request; `./tests/conftest.py` contains pytest fixtures which may be of use. Test files and test data can be found in the `tests/` folder.

## Documentation

Github Actions is used to automatically [build and publish code documentation](https://github.com/NU-CEM/ThermoPot/actions/workflows/build-docs.yml) after each commit to the main branch of the ThermoPot repository. All documentation is written in Markdown using the [MkDocs framework](https://www.mkdocs.org/) and the [Material theme](https://squidfunk.github.io/mkdocs-material/). This also includes the automatic generation of API documentation using the [mkdocstrings](https://github.com/mkdocstrings/mkdocstrings) plugin. 

If possible, please update the documentation when issuing a Pull Request. This includes both the web pages written in Markdown, and in-code docstrings. The latter are particularly important, and we request that these are written [Google-style](https://gist.github.com/redlotus/3bc387c2591e3e908c9b63b97b11d24e).

All documentation source can be found in the `docs/` folder. MkDocs is configured with the `/.mkdocs.yml` folder in the project root.

## Linting

[Black](https://github.com/psf/black) is used to automatically [re-format all Python code](https://github.com/NU-CEM/ThermoPot/actions/workflows/lint-code.yml) committed to the repository. The code will be formatted according to the PEP8 specification. As this is an automatic step you do not need lint/format your code to any specification before issuing a PR.

## Publishing

Github Actions (surprise, surprise!) is used to automatically [build and publish the ThermoPot package](https://github.com/NU-CEM/ThermoPot/actions/workflows/build-n-publish.yml) after every tagged release. Before tagging the release, we need to remember to update the version number in `__init__.py`. It's so easy to forget!
