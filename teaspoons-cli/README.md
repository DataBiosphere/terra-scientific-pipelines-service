# Teaspoons CLI

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
To generate the Python thick client, run the following command:
```bash
./gradlew :teaspoons-cli:cliBuild
```

## Development
You'll need to have poetry installed to manage python dependencies. Instructions for installing poetry can be found [here](https://python-poetry.org/docs/).
