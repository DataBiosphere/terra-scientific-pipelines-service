# Teaspoons CLI

## Python CLI structure
The CLI code is structured as follows:
```
teaspoons-cli
├── teaspoons
│   └── generated
│   │   └── [auto-generated files for thin client]
│   └── __init__.py
│   └── auth.py
│   └── cli.py
│   └── config.py
│   └── teaspoons
└── tests (to be created)
├── pyproject.toml
├── poetry.lock
├── README.md
```

Inside the `teaspoons` directory, we have the following files:
- `teaspoons` is the entrypoint for the CLI. It contains the main function that is called when the CLI is run.
- `auth.py` contains the code for authenticating with the Teaspoons service (Terra, via b2c).
- `config.py` contains the code for managing the CLI configuration via environment variables.
- `cli.py` contains the code for the CLI commands and options. This is like the controller layer for the CLI.
- A future file will be included to contain the business logic for the CLI commands.
- The `generated` directory contains the auto-generated files for the thin client, containing the python model classes and API calls.

## Development
You'll need to have poetry installed to manage python dependencies. Instructions for installing poetry can be found [here](https://python-poetry.org/docs/).

## Python thin client auto-generation
We use the [OpenAPI Generator](https://github.com/OpenAPITools/openapi-generator) to generate the "thin" Python client,
which is then used to build the Python-based "thick" CLI tool.

To generate the Python thin client, run the following command:
```bash
./gradlew :teaspoons-cli:openApiGenerate
```

This will produce generated files at `/teaspoons-cli/generated/`.

Note we do not run the openApiGenerate task as part of the main Teaspoons build, as it is not necessary for the 
service itself and we don't want any potential bugs in the CLI to affect the service.

## Python thick client auto-generation
To generate the Python thick client (which will also generate the thin client), run the following command:
```bash
./gradlew :teaspoons-cli:cliBuild
```
