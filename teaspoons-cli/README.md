# Teaspoons CLI

## Python CLI structure
The CLI code is structured as follows:
```
teaspoons-cli
├── teaspoons
│   └── commands
│   │   └── __init__.py
│   │   └── auth.py
│   │   └── pipelines.py
│   └── generated
│   │   └── [auto-generated files for thin client]
├── └── tests (to be created)
│   └── __init__.py
│   └── auth.py
│   └── cli.py
│   └── config.py
│   └── teaspoons
├── pyproject.toml
├── poetry.lock
├── README.md
```

Inside the `teaspoons` directory, we have the following files:
- `teaspoons` is the entrypoint for the CLI. It contains the main function that is called when the CLI is run.
- `auth.py` contains the code for authenticating with the Teaspoons service (Terra, via b2c).
- `config.py` contains the code for managing the CLI configuration via environment variables.
- `cli.py` assembles the CLI sub-modules that are defined in `commands/`. 
- A future file will be included to contain the business logic for the CLI commands.
- The `commands` directory contains the CLI sub-modules. This is effectively the controller layer for the CLI.
- The `generated` directory contains the auto-generated files for the thin client, containing the python model classes and API calls.


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

This will produce generated files at `/teaspoons-cli/generated/`.

Note we do not run the openApiGenerate task as part of the main Teaspoons build, as it is not necessary for the 
service itself and we don't want any potential bugs in the CLI to affect the service.

