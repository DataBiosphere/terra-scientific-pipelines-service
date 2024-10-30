# Teaspoons CLI

## Python CLI structure
The CLI code is structured as follows:
```
teaspoons-cli
└── generated
│   └── [auto-generated files for thin client]
├── teaspoons
│   └── commands
│   │   └── __init__.py
│   │   └── auth_commands.py
│   │   └── pipelines_commands.py
│   └── logic
│   │   └── __init__.py
│   │   └── auth_logic.py
│   │   └── pipelines_logic.py
├── tests
│   └── commands
│   │   └── test_auth_commands.py
│   │   └── test_pipelines_commands.py
│   └── logic
│   │   └── test_auth_logic.py
│   │   └── test_pipelines_logic.py
│   └── __init__.py
│   └── auth_helper.py
│   └── cli.py
│   └── client.py
│   └── config.py
│   └── teaspoons
├── pyproject.toml
├── poetry.lock
├── README.md
```

The `generated` directory contains the auto-generated files for the thin client, including the python model classes and API calls.

In the `teaspoons` directory, we have the following files and subdirectories:
- `auth_helper.py` contains the code for authenticating with the Teaspoons service (Terra, via b2c).
- `cli.py` assembles the CLI sub-modules that are defined in `commands/`.
- `client.py` contains the code for wrapping API calls to the Teaspoons service.
- `config.py` contains the code for managing the CLI configuration via environment variables.
- `teaspoons` is the entrypoint for the CLI. It contains the main function that is called when the CLI is run.
- `utils.py` contains utility functions that are used across the CLI.
- The `commands` directory contains the CLI sub-modules. This is effectively the controller layer for the CLI.
- The `logic` directory contains the business logic for the CLI.

In the `tests` directory, we have test files that can be run with pytest.

## Using the CLI
For now, the CLI requires poetry to be installed to run. See the [Development](#development) section for instructions on how to install poetry.

To run the CLI, navigate to the `teaspoons-cli/teaspoons/` directory and run the following command:
```bash
./teaspoons COMMAND [ARGS]
```

For example, to authenticate with the Teaspoons service, run the following command:
```bash
./teaspoons auth login
```

To list the pipelines in the Teaspoons service, run the following command:
```bash
./teaspoons pipelines list
```

See WIP documentation for the CLI [here](https://docs.google.com/document/d/1ovbcHCzdyuC8RjFfkVJZiuDTQ_UAVrglSxSGaZwppoY/edit?tab=t.0#heading=h.jfsr3j3x0zjr).


## Development
You'll need to have poetry installed to manage python dependencies. Instructions for installing poetry can be found [here](https://python-poetry.org/docs/).

To run the tests, execute the following command from the teaspoons-cli directory:
```bash
poetry run pytest
```

To run the formatter, execute the following command from the teaspoons-cli directory:
```bash
poetry run black .
```

To run the linter with fixes, execute the following command from the teaspoons-cli directory:
```bash
poetry run ruff check --fix
```
To run the linter as a check without fixes, omit the `--fix` flag.


## Python thick client auto-generation
To generate the Python thick client (which will also generate the thin client), run the following command:
```bash
./gradlew :teaspoons-cli:cliBuild
```


## Python thin client auto-generation
Note: If you've already built the thick client, you don't need to generate the thin client separately.

We use the [OpenAPI Generator](https://github.com/OpenAPITools/openapi-generator) to generate the "thin" Python client,
which is then used to build the Python-based "thick" CLI tool.

To generate the Python thin client, run the following command:
```bash
./gradlew :teaspoons-cli:openApiGenerate
```

This will produce generated files at `/teaspoons-cli/teaspoons/generated/`.

Note we do not run the openApiGenerate task as part of the main Teaspoons build, as it is not necessary for the 
service itself and we don't want any potential bugs in the CLI to affect the service.

